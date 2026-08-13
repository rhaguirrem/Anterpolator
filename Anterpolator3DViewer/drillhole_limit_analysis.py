from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy import stats


MISSING_CATEGORY_LABEL = "(blank)"
BETWEEN_BOUNDARIES_LABEL = "(between_boundaries)"


@dataclass(frozen=True)
class PreparedIntervalData:
    frame: pd.DataFrame
    numeric_aliases: dict[str, str]
    categorical_aliases: dict[str, str]
    quality_summary: pd.DataFrame


ProgressCallback = Callable[[int, int, str], None]


def _coerce_numeric_series(series: pd.Series) -> pd.Series:
    numeric_series = pd.to_numeric(series, errors="coerce")
    missing_mask = numeric_series.isna()
    if not missing_mask.any():
        return numeric_series

    text_series = series.astype(str)
    normalized = text_series.str.replace("\u00a0", "", regex=False).str.replace(" ", "", regex=False)
    alternate_numeric = pd.to_numeric(normalized, errors="coerce")
    return numeric_series.where(~missing_mask, alternate_numeric)


def _normalize_category_series(series: pd.Series) -> pd.Series:
    text_series = series.where(~series.isna(), "").astype(str).str.strip()
    return text_series.mask(text_series == "", MISSING_CATEGORY_LABEL)


def _weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    if len(values) == 0:
        return float("nan")
    total_weight = float(weights.sum())
    if total_weight <= 0.0:
        return float(values.mean())
    return float(np.average(values.to_numpy(dtype=float), weights=weights.to_numpy(dtype=float)))


def _weighted_majority(values: pd.Series, weights: pd.Series) -> str:
    if len(values) == 0:
        return MISSING_CATEGORY_LABEL

    categories = _normalize_category_series(values)
    safe_weights = pd.to_numeric(weights, errors="coerce").fillna(0.0).clip(lower=0.0)

    totals: dict[str, float] = {}
    counts: dict[str, int] = {}
    for category_value, weight_value in zip(categories.tolist(), safe_weights.tolist()):
        category_text = str(category_value)
        totals[category_text] = totals.get(category_text, 0.0) + float(weight_value)
        counts[category_text] = counts.get(category_text, 0) + 1

    return sorted(totals.keys(), key=lambda item: (-totals[item], -counts[item], item))[0]


def _normalize_boundary_categories(boundary_categories: list[str] | tuple[str, ...] | None) -> set[str] | None:
    if not boundary_categories:
        return None
    normalized = _normalize_category_series(pd.Series(list(boundary_categories), dtype="object"))
    return {str(value) for value in normalized.tolist()}


def _quality_rows_from_counts(counts: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(
        [{"metric": str(key), "value": value} for key, value in counts.items()],
        columns=["metric", "value"],
    )


def _emit_progress(progress_callback: ProgressCallback | None, current: int, total: int, message: str) -> None:
    if progress_callback is None:
        return
    progress_callback(int(current), int(total), str(message))


def prepare_drillhole_interval_data(
    df: pd.DataFrame,
    hole_id_col: str,
    from_col: str,
    to_col: str,
    categorical_fields: list[str] | tuple[str, ...],
    numeric_fields: list[str] | tuple[str, ...],
) -> PreparedIntervalData:
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    required_columns = [hole_id_col, from_col, to_col, *categorical_fields, *numeric_fields]
    missing_columns = [column for column in required_columns if column not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    work = df.loc[:, required_columns].copy()
    work["__hole_id"] = work[hole_id_col].where(~work[hole_id_col].isna(), "").astype(str).str.strip()
    from_values = _coerce_numeric_series(work[from_col])
    to_values = _coerce_numeric_series(work[to_col])
    work["__from"] = np.minimum(from_values, to_values)
    work["__to"] = np.maximum(from_values, to_values)
    work["__interval_length"] = work["__to"] - work["__from"]

    numeric_aliases: dict[str, str] = {}
    for field in numeric_fields:
        alias = f"__numeric__{field}"
        work[alias] = _coerce_numeric_series(work[field])
        numeric_aliases[str(field)] = alias

    categorical_aliases: dict[str, str] = {}
    for field in categorical_fields:
        alias = f"__category__{field}"
        work[alias] = _normalize_category_series(work[field])
        categorical_aliases[str(field)] = alias

    initial_row_count = int(len(work))
    valid_mask = (
        work["__hole_id"].ne("")
        & work["__from"].notna()
        & work["__to"].notna()
        & (work["__interval_length"] > 0.0)
    )
    prepared = work.loc[valid_mask].copy()

    quality_counts = {
        "input_rows": initial_row_count,
        "valid_interval_rows": int(valid_mask.sum()),
        "dropped_rows": int((~valid_mask).sum()),
        "rows_missing_hole_id": int(work["__hole_id"].eq("").sum()),
        "rows_missing_interval_bounds": int((work["__from"].isna() | work["__to"].isna()).sum()),
        "rows_non_positive_interval": int((work["__interval_length"] <= 0.0).sum()),
    }
    for field in numeric_fields:
        alias = numeric_aliases[field]
        quality_counts[f"rows_with_numeric_{field}"] = int(prepared[alias].notna().sum())

    return PreparedIntervalData(
        frame=prepared,
        numeric_aliases=numeric_aliases,
        categorical_aliases=categorical_aliases,
        quality_summary=_quality_rows_from_counts(quality_counts),
    )


def _summarize_groups(
    subset: pd.DataFrame,
    limit_field: str,
    numeric_field: str,
    category_column: str,
    numeric_column: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for category_value, group in subset.groupby(category_column, dropna=False, sort=True):
        numeric_values = group[numeric_column]
        interval_lengths = group["__interval_length"]
        rows.append(
            {
                "limit_field": limit_field,
                "numeric_field": numeric_field,
                "category": str(category_value),
                "row_count": int(len(group)),
                "hole_count": int(group["__hole_id"].nunique()),
                "total_interval_length": float(interval_lengths.sum()),
                "mean": float(numeric_values.mean()),
                "weighted_mean": _weighted_mean(numeric_values, interval_lengths),
                "median": float(numeric_values.median()),
                "std": float(numeric_values.std(ddof=1)) if len(group) > 1 else float("nan"),
                "min": float(numeric_values.min()),
                "max": float(numeric_values.max()),
                "q25": float(numeric_values.quantile(0.25)),
                "q75": float(numeric_values.quantile(0.75)),
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(
            columns=[
                "limit_field",
                "numeric_field",
                "category",
                "row_count",
                "hole_count",
                "total_interval_length",
                "mean",
                "weighted_mean",
                "median",
                "std",
                "min",
                "max",
                "q25",
                "q75",
            ]
        )
    return result.sort_values(["limit_field", "numeric_field", "weighted_mean", "category"], ascending=[True, True, False, True]).reset_index(drop=True)


def _run_statistical_tests(
    subset: pd.DataFrame,
    limit_field: str,
    numeric_field: str,
    category_column: str,
    numeric_column: str,
    min_group_size: int,
) -> pd.DataFrame:
    result_columns = ["limit_field", "numeric_field", "test_name", "groups_tested", "min_group_size", "statistic", "p_value", "significant_0_05"]
    groups: list[tuple[str, np.ndarray]] = []
    for category_value, group in subset.groupby(category_column, dropna=False, sort=True):
        numeric_values = group[numeric_column].dropna().to_numpy(dtype=float)
        if len(numeric_values) > 0:
            groups.append((str(category_value), numeric_values))

    rows: list[dict[str, Any]] = []
    if len(groups) < 2:
        return pd.DataFrame(columns=result_columns)

    eligible_groups = [values for _name, values in groups if len(values) >= min_group_size]
    tested_group_names = " | ".join(name for name, _values in groups)
    if len(eligible_groups) >= 2:
        kruskal_stat, kruskal_p = stats.kruskal(*eligible_groups, nan_policy="omit")
        rows.append(
            {
                "limit_field": limit_field,
                "numeric_field": numeric_field,
                "test_name": "kruskal",
                "groups_tested": tested_group_names,
                "min_group_size": int(min_group_size),
                "statistic": float(kruskal_stat),
                "p_value": float(kruskal_p),
                "significant_0_05": bool(kruskal_p < 0.05),
            }
        )

    if len(groups) == 2 and all(len(values) >= 2 for _name, values in groups):
        t_stat, t_p = stats.ttest_ind(groups[0][1], groups[1][1], equal_var=False, nan_policy="omit")
        mw_stat, mw_p = stats.mannwhitneyu(groups[0][1], groups[1][1], alternative="two-sided")
        rows.extend(
            [
                {
                    "limit_field": limit_field,
                    "numeric_field": numeric_field,
                    "test_name": "welch_t",
                    "groups_tested": tested_group_names,
                    "min_group_size": 2,
                    "statistic": float(t_stat),
                    "p_value": float(t_p),
                    "significant_0_05": bool(t_p < 0.05),
                },
                {
                    "limit_field": limit_field,
                    "numeric_field": numeric_field,
                    "test_name": "mann_whitney_u",
                    "groups_tested": tested_group_names,
                    "min_group_size": 2,
                    "statistic": float(mw_stat),
                    "p_value": float(mw_p),
                    "significant_0_05": bool(mw_p < 0.05),
                },
            ]
        )

    if not rows:
        return pd.DataFrame(columns=result_columns)
    return pd.DataFrame(rows, columns=result_columns)


def _summarize_boundary_transitions(
    subset: pd.DataFrame,
    limit_field: str,
    numeric_field: str,
    category_column: str,
    numeric_column: str,
) -> pd.DataFrame:
    events: list[dict[str, Any]] = []
    ordered = subset.sort_values(["__hole_id", "__from", "__to"], ascending=[True, True, True]).reset_index(drop=True)

    for _hole_id, hole_group in ordered.groupby("__hole_id", sort=False):
        previous_row: pd.Series | None = None
        for _index, row in hole_group.iterrows():
            if previous_row is not None and previous_row[category_column] != row[category_column]:
                events.append(
                    {
                        "from_category": str(previous_row[category_column]),
                        "to_category": str(row[category_column]),
                        "hole_id": str(row["__hole_id"]),
                        "delta": float(row[numeric_column] - previous_row[numeric_column]),
                        "abs_delta": float(abs(row[numeric_column] - previous_row[numeric_column])),
                        "gap": float(row["__from"] - previous_row["__to"]),
                    }
                )
            previous_row = row

    if not events:
        return pd.DataFrame(
            columns=[
                "limit_field",
                "numeric_field",
                "from_category",
                "to_category",
                "boundary_count",
                "hole_count",
                "mean_delta",
                "median_delta",
                "mean_abs_delta",
                "std_delta",
                "min_delta",
                "max_delta",
                "mean_gap",
            ]
        )

    events_df = pd.DataFrame(events)
    rows: list[dict[str, Any]] = []
    for (from_category, to_category), group in events_df.groupby(["from_category", "to_category"], sort=True):
        rows.append(
            {
                "limit_field": limit_field,
                "numeric_field": numeric_field,
                "from_category": str(from_category),
                "to_category": str(to_category),
                "boundary_count": int(len(group)),
                "hole_count": int(group["hole_id"].nunique()),
                "mean_delta": float(group["delta"].mean()),
                "median_delta": float(group["delta"].median()),
                "mean_abs_delta": float(group["abs_delta"].mean()),
                "std_delta": float(group["delta"].std(ddof=1)) if len(group) > 1 else float("nan"),
                "min_delta": float(group["delta"].min()),
                "max_delta": float(group["delta"].max()),
                "mean_gap": float(group["gap"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(["limit_field", "numeric_field", "mean_abs_delta"], ascending=[True, True, False]).reset_index(drop=True)


def _build_merged_interval_row(
    segment_rows: pd.DataFrame,
    prepared: PreparedIntervalData,
    limit_field: str,
    numeric_fields: list[str] | tuple[str, ...],
    composite_categorical_fields: list[str],
    segment_index: int,
    segment_role: str,
    analysis_category: str,
) -> dict[str, Any]:
    segment_start = float(segment_rows["__from"].min())
    segment_end = float(segment_rows["__to"].max())
    covered_length = float(segment_rows["__interval_length"].sum())
    span_length = float(segment_end - segment_start)

    row: dict[str, Any] = {
        "limit_field": limit_field,
        "hole_id": str(segment_rows.iloc[0]["__hole_id"]),
        "segment_index": int(segment_index),
        "segment_role": segment_role,
        "analysis_category": str(analysis_category),
        "from": segment_start,
        "to": segment_end,
        "covered_length": covered_length,
        "span_length": span_length,
        "gap_length": float(max(0.0, span_length - covered_length)),
        "source_interval_count": int(len(segment_rows)),
        "__hole_id": str(segment_rows.iloc[0]["__hole_id"]),
        "__from": segment_start,
        "__to": segment_end,
        "__interval_length": covered_length,
    }

    limit_alias = prepared.categorical_aliases[limit_field]
    row["majority_limit_category"] = _weighted_majority(segment_rows[limit_alias], segment_rows["__interval_length"])

    for field in composite_categorical_fields:
        alias = prepared.categorical_aliases[field]
        row[field] = _weighted_majority(segment_rows[alias], segment_rows["__interval_length"])

    for field in numeric_fields:
        alias = prepared.numeric_aliases[field]
        valid_values = segment_rows.loc[segment_rows[alias].notna(), [alias, "__interval_length"]]
        row[field] = (
            _weighted_mean(valid_values[alias], valid_values["__interval_length"])
            if not valid_values.empty
            else float("nan")
        )

    return row


def _merge_intervals_for_limit_field(
    prepared: PreparedIntervalData,
    limit_field: str,
    numeric_fields: list[str] | tuple[str, ...],
    composite_categorical_fields: list[str],
    boundary_categories: list[str] | tuple[str, ...] | None,
    skip_holes_without_boundary: bool = True,
    hole_progress_callback: Callable[[str], None] | None = None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    limit_alias = prepared.categorical_aliases[limit_field]
    normalized_boundary_categories = _normalize_boundary_categories(boundary_categories)

    selected_columns = ["__hole_id", "__from", "__to", "__interval_length", limit_alias]
    selected_columns.extend(prepared.numeric_aliases[field] for field in numeric_fields)
    selected_columns.extend(prepared.categorical_aliases[field] for field in composite_categorical_fields)
    ordered = prepared.frame.loc[:, list(dict.fromkeys(selected_columns))].sort_values(
        ["__hole_id", "__from", "__to"], ascending=[True, True, True]
    )

    merged_rows: list[dict[str, Any]] = []
    total_holes = int(ordered["__hole_id"].nunique())
    holes_with_boundary = 0
    holes_processed = 0
    holes_skipped_without_boundary = 0

    for hole_index, (_hole_id, hole_group) in enumerate(ordered.groupby("__hole_id", sort=False), start=1):
        current_indices: list[int] = []
        current_category = MISSING_CATEGORY_LABEL
        segment_index = 0

        if normalized_boundary_categories is not None:
            hole_categories = {str(value) for value in hole_group[limit_alias].tolist()}
            has_boundary = bool(hole_categories & normalized_boundary_categories)
            if has_boundary:
                holes_with_boundary += 1
            elif skip_holes_without_boundary:
                holes_skipped_without_boundary += 1
                if hole_progress_callback is not None:
                    hole_progress_callback(f"Skipping hole {hole_index:,}/{total_holes:,}: {str(_hole_id)} has no selected boundary category")
                continue

        holes_processed += 1
        if hole_progress_callback is not None:
            hole_progress_callback(f"Processing hole {hole_index:,}/{total_holes:,}: {str(_hole_id)}")

        def flush_current(segment_role: str, analysis_category: str) -> int:
            nonlocal current_indices
            if not current_indices:
                return 0
            merged_rows.append(
                _build_merged_interval_row(
                    hole_group.loc[current_indices],
                    prepared=prepared,
                    limit_field=limit_field,
                    numeric_fields=numeric_fields,
                    composite_categorical_fields=composite_categorical_fields,
                    segment_index=segment_index,
                    segment_role=segment_role,
                    analysis_category=analysis_category,
                )
            )
            current_indices = []
            return 1

        for row_index, row in hole_group.iterrows():
            limit_value = str(row[limit_alias])
            if normalized_boundary_categories is not None:
                if limit_value in normalized_boundary_categories:
                    segment_index += flush_current("composite", BETWEEN_BOUNDARIES_LABEL)
                    merged_rows.append(
                        _build_merged_interval_row(
                            hole_group.loc[[row_index]],
                            prepared=prepared,
                            limit_field=limit_field,
                            numeric_fields=numeric_fields,
                            composite_categorical_fields=composite_categorical_fields,
                            segment_index=segment_index,
                            segment_role="boundary",
                            analysis_category=limit_value,
                        )
                    )
                    segment_index += 1
                    continue

                current_indices.append(row_index)
                continue

            if current_indices and limit_value != current_category:
                segment_index += flush_current("category_run", current_category)

            current_indices.append(row_index)
            current_category = limit_value

        if normalized_boundary_categories is not None:
            segment_index += flush_current("composite", BETWEEN_BOUNDARIES_LABEL)
        else:
            segment_index += flush_current("category_run", current_category)

    if not merged_rows:
        return pd.DataFrame(
            columns=[
                "limit_field",
                "hole_id",
                "segment_index",
                "segment_role",
                "analysis_category",
                "majority_limit_category",
                "from",
                "to",
                "covered_length",
                "span_length",
                "gap_length",
                "source_interval_count",
                "__hole_id",
                "__from",
                "__to",
                "__interval_length",
                *composite_categorical_fields,
                *numeric_fields,
            ]
        ), {
            "holes_total": total_holes,
            "holes_processed": holes_processed,
            "holes_with_boundary": holes_with_boundary,
            "holes_skipped_without_boundary": holes_skipped_without_boundary,
        }

    return pd.DataFrame(merged_rows), {
        "holes_total": total_holes,
        "holes_processed": holes_processed,
        "holes_with_boundary": holes_with_boundary,
        "holes_skipped_without_boundary": holes_skipped_without_boundary,
    }


def analyze_drillhole_limit_behavior(
    df: pd.DataFrame,
    hole_id_col: str,
    from_col: str,
    to_col: str,
    limit_fields: list[str] | tuple[str, ...],
    numeric_fields: list[str] | tuple[str, ...],
    min_group_size: int = 3,
    boundary_categories_by_field: dict[str, list[str] | tuple[str, ...]] | None = None,
    composite_categorical_fields: list[str] | tuple[str, ...] | None = None,
    skip_holes_without_boundary: bool = True,
    progress_callback: ProgressCallback | None = None,
) -> dict[str, pd.DataFrame]:
    if not limit_fields:
        raise ValueError("At least one categorical limit field is required.")
    if not numeric_fields:
        raise ValueError("At least one numeric field is required.")

    composite_field_list = list(dict.fromkeys([*limit_fields, *(composite_categorical_fields or [])]))
    boundary_categories_by_field = dict(boundary_categories_by_field or {})

    prepared = prepare_drillhole_interval_data(
        df,
        hole_id_col=hole_id_col,
        from_col=from_col,
        to_col=to_col,
        categorical_fields=composite_field_list,
        numeric_fields=list(numeric_fields),
    )

    hole_count = int(prepared.frame["__hole_id"].nunique()) if not prepared.frame.empty else 0
    total_steps = max(1, 1 + len(limit_fields) * hole_count + len(limit_fields) * len(numeric_fields) * 3)
    completed_steps = 0

    def report_progress(message: str, advance: int = 1) -> None:
        nonlocal completed_steps
        completed_steps = min(total_steps, completed_steps + max(0, int(advance)))
        _emit_progress(progress_callback, completed_steps, total_steps, message)

    report_progress(
        f"Prepared {len(prepared.frame):,} valid intervals across {hole_count:,} drillholes",
        advance=1,
    )

    summary_rows: list[dict[str, Any]] = []
    merged_frames: list[pd.DataFrame] = []
    group_frames: list[pd.DataFrame] = []
    test_frames: list[pd.DataFrame] = []
    boundary_frames: list[pd.DataFrame] = []
    quality_counts = dict(zip(prepared.quality_summary["metric"], prepared.quality_summary["value"]))

    for limit_field in limit_fields:
        merged_frame, merge_stats = _merge_intervals_for_limit_field(
            prepared,
            limit_field=limit_field,
            numeric_fields=list(numeric_fields),
            composite_categorical_fields=composite_field_list,
            boundary_categories=boundary_categories_by_field.get(limit_field),
            skip_holes_without_boundary=skip_holes_without_boundary,
            hole_progress_callback=lambda message: report_progress(message, advance=1),
        )
        merged_frames.append(merged_frame)
        quality_counts[f"holes_total_{limit_field}"] = merge_stats["holes_total"]
        quality_counts[f"holes_processed_{limit_field}"] = merge_stats["holes_processed"]
        quality_counts[f"holes_with_boundary_{limit_field}"] = merge_stats["holes_with_boundary"]
        quality_counts[f"holes_skipped_without_boundary_{limit_field}"] = merge_stats["holes_skipped_without_boundary"]

        for numeric_field in numeric_fields:
            subset = merged_frame.loc[
                merged_frame[numeric_field].notna(),
                ["__hole_id", "__from", "__to", "__interval_length", "analysis_category", numeric_field],
            ].copy()

            group_summary = _summarize_groups(subset, limit_field, numeric_field, "analysis_category", numeric_field)
            report_progress(f"Summarized groups for {limit_field} / {numeric_field}", advance=1)
            statistical_tests = _run_statistical_tests(subset, limit_field, numeric_field, "analysis_category", numeric_field, min_group_size=min_group_size)
            report_progress(f"Ran statistical tests for {limit_field} / {numeric_field}", advance=1)
            boundary_summary = _summarize_boundary_transitions(subset, limit_field, numeric_field, "analysis_category", numeric_field)
            report_progress(f"Computed boundary transitions for {limit_field} / {numeric_field}", advance=1)

            group_frames.append(group_summary)
            test_frames.append(statistical_tests)
            boundary_frames.append(boundary_summary)

            if group_summary.empty:
                summary_rows.append(
                    {
                        "limit_field": limit_field,
                        "numeric_field": numeric_field,
                        "valid_rows": 0,
                        "category_count": 0,
                        "largest_mean_gap": float("nan"),
                        "best_supported_test": "n/a",
                        "best_p_value": float("nan"),
                        "significant_0_05": False,
                        "boundary_transition_count": 0,
                    }
                )
                continue

            ordered_groups = group_summary.sort_values("weighted_mean", ascending=False).reset_index(drop=True)
            largest_gap = float(ordered_groups.loc[0, "weighted_mean"] - ordered_groups.loc[len(ordered_groups) - 1, "weighted_mean"]) if len(ordered_groups) >= 2 else 0.0
            best_test_row = statistical_tests.sort_values("p_value", ascending=True).head(1)
            boundary_count = int(boundary_summary["boundary_count"].sum()) if not boundary_summary.empty else 0
            summary_rows.append(
                {
                    "limit_field": limit_field,
                    "numeric_field": numeric_field,
                    "valid_rows": int(len(subset)),
                    "category_count": int(group_summary["category"].nunique()),
                    "largest_mean_gap": largest_gap,
                    "best_supported_test": str(best_test_row.iloc[0]["test_name"]) if not best_test_row.empty else "n/a",
                    "best_p_value": float(best_test_row.iloc[0]["p_value"]) if not best_test_row.empty else float("nan"),
                    "significant_0_05": bool(bool(best_test_row.iloc[0]["p_value"] < 0.05) if not best_test_row.empty else False),
                    "boundary_transition_count": boundary_count,
                }
            )

    quality_summary = _quality_rows_from_counts(quality_counts)
    _emit_progress(progress_callback, total_steps, total_steps, "Analysis complete")

    return {
        "quality_summary": quality_summary,
        "merged_intervals": pd.concat(merged_frames, ignore_index=True) if merged_frames else pd.DataFrame(),
        "summary": pd.DataFrame(summary_rows).sort_values(["significant_0_05", "largest_mean_gap"], ascending=[False, False]).reset_index(drop=True),
        "group_summary": pd.concat(group_frames, ignore_index=True) if group_frames else pd.DataFrame(),
        "statistical_tests": pd.concat(test_frames, ignore_index=True) if test_frames else pd.DataFrame(),
        "boundary_transitions": pd.concat(boundary_frames, ignore_index=True) if boundary_frames else pd.DataFrame(),
    }