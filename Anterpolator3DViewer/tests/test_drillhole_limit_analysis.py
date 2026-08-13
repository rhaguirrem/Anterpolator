import os
import sys

import pandas as pd
import pytest

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from drillhole_limit_analysis import analyze_drillhole_limit_behavior, prepare_drillhole_interval_data


def _build_behavior_frame() -> pd.DataFrame:
    rows = []
    for hole_index in range(12):
        hole_id = f"DH_{hole_index:02d}"
        rows.extend(
            [
                {"hole_id": hole_id, "from": 0.0, "to": 10.0, "fault_zone": "Host", "grade": 1.0 + hole_index * 0.03, "density": 2.60},
                {"hole_id": hole_id, "from": 10.0, "to": 20.0, "fault_zone": "Fault", "grade": 4.8 + hole_index * 0.05, "density": 2.45},
                {"hole_id": hole_id, "from": 20.0, "to": 30.0, "fault_zone": "Host", "grade": 1.2 + hole_index * 0.03, "density": 2.59},
            ]
        )
    return pd.DataFrame(rows)


def test_prepare_drillhole_interval_data_drops_invalid_intervals():
    df = pd.DataFrame(
        [
            {"hole_id": "A", "from": 0, "to": 5, "fault_zone": "Host", "grade": 1.1},
            {"hole_id": "", "from": 5, "to": 10, "fault_zone": "Fault", "grade": 5.0},
            {"hole_id": "A", "from": 10, "to": 10, "fault_zone": "Fault", "grade": 5.3},
            {"hole_id": "A", "from": "bad", "to": 15, "fault_zone": "Fault", "grade": 4.9},
        ]
    )

    prepared = prepare_drillhole_interval_data(
        df,
        hole_id_col="hole_id",
        from_col="from",
        to_col="to",
        categorical_fields=["fault_zone"],
        numeric_fields=["grade"],
    )

    assert len(prepared.frame) == 1
    quality = dict(zip(prepared.quality_summary["metric"], prepared.quality_summary["value"]))
    assert quality["valid_interval_rows"] == 1
    assert quality["dropped_rows"] == 3


def test_analyze_drillhole_limit_behavior_detects_group_difference_and_boundaries():
    df = _build_behavior_frame()

    result = analyze_drillhole_limit_behavior(
        df,
        hole_id_col="hole_id",
        from_col="from",
        to_col="to",
        limit_fields=["fault_zone"],
        numeric_fields=["grade"],
        min_group_size=3,
        boundary_categories_by_field={"fault_zone": ["Fault"]},
    )

    summary = result["summary"]
    assert len(summary) == 1
    assert bool(summary.iloc[0]["significant_0_05"])
    assert float(summary.iloc[0]["largest_mean_gap"]) > 3.0

    group_summary = result["group_summary"]
    host_mean = float(group_summary.loc[group_summary["category"] == "(between_boundaries)", "weighted_mean"].iloc[0])
    fault_mean = float(group_summary.loc[group_summary["category"] == "Fault", "weighted_mean"].iloc[0])
    assert fault_mean > host_mean

    tests = result["statistical_tests"]
    mann_whitney = tests.loc[tests["test_name"] == "mann_whitney_u"].iloc[0]
    assert float(mann_whitney["p_value"]) < 0.05

    boundaries = result["boundary_transitions"]
    host_to_fault = boundaries.loc[
        (boundaries["from_category"] == "(between_boundaries)") & (boundaries["to_category"] == "Fault")
    ].iloc[0]
    fault_to_host = boundaries.loc[
        (boundaries["from_category"] == "Fault") & (boundaries["to_category"] == "(between_boundaries)")
    ].iloc[0]
    assert int(host_to_fault["boundary_count"]) == 12
    assert int(fault_to_host["boundary_count"]) == 12
    assert float(host_to_fault["mean_delta"]) > 0.0
    assert float(fault_to_host["mean_delta"]) < 0.0


def test_analyze_drillhole_limit_behavior_merges_non_boundary_runs_with_weighted_averages():
    df = pd.DataFrame(
        [
            {"hole_id": "DH_01", "from": 0.0, "to": 5.0, "fault_zone": "Host", "rock": "Host", "grade": 1.0, "density": 2.0},
            {"hole_id": "DH_01", "from": 5.0, "to": 15.0, "fault_zone": "Host", "rock": "Host", "grade": 3.0, "density": 2.2},
            {"hole_id": "DH_01", "from": 15.0, "to": 17.0, "fault_zone": "Fault", "rock": "Fault", "grade": 8.0, "density": 2.5},
            {"hole_id": "DH_01", "from": 17.0, "to": 20.0, "fault_zone": "Host", "rock": "Host", "grade": 4.0, "density": 2.1},
            {"hole_id": "DH_01", "from": 20.0, "to": 26.0, "fault_zone": "Host", "rock": "Altered", "grade": 6.0, "density": 2.3},
            {"hole_id": "DH_01", "from": 26.0, "to": 31.0, "fault_zone": "Host", "rock": "Host", "grade": 5.0, "density": 2.4},
            {"hole_id": "DH_01", "from": 31.0, "to": 33.0, "fault_zone": "Fault", "rock": "Fault", "grade": 9.0, "density": 2.6},
            {"hole_id": "DH_01", "from": 33.0, "to": 40.0, "fault_zone": "Host", "rock": "Altered", "grade": 7.0, "density": 2.5},
        ]
    )

    result = analyze_drillhole_limit_behavior(
        df,
        hole_id_col="hole_id",
        from_col="from",
        to_col="to",
        limit_fields=["fault_zone"],
        numeric_fields=["grade", "density"],
        min_group_size=1,
        boundary_categories_by_field={"fault_zone": ["Fault"]},
        composite_categorical_fields=["rock"],
    )

    merged = result["merged_intervals"].reset_index(drop=True)
    assert list(merged["segment_role"]) == ["composite", "boundary", "composite", "boundary", "composite"]
    assert list(merged["source_interval_count"]) == [2, 1, 3, 1, 1]
    assert list(merged["analysis_category"]) == ["(between_boundaries)", "Fault", "(between_boundaries)", "Fault", "(between_boundaries)"]

    assert merged.loc[0, "from"] == pytest.approx(0.0)
    assert merged.loc[0, "to"] == pytest.approx(15.0)
    assert merged.loc[0, "grade"] == pytest.approx((1.0 * 5.0 + 3.0 * 10.0) / 15.0)
    assert merged.loc[0, "density"] == pytest.approx((2.0 * 5.0 + 2.2 * 10.0) / 15.0)
    assert merged.loc[0, "rock"] == "Host"
    assert merged.loc[0, "majority_limit_category"] == "Host"

    assert merged.loc[2, "from"] == pytest.approx(17.0)
    assert merged.loc[2, "to"] == pytest.approx(31.0)
    assert merged.loc[2, "grade"] == pytest.approx((4.0 * 3.0 + 6.0 * 6.0 + 5.0 * 5.0) / 14.0)
    assert merged.loc[2, "rock"] == "Host"

    assert merged.loc[4, "grade"] == pytest.approx(7.0)


def test_analyze_drillhole_limit_behavior_skips_holes_without_boundary_by_default():
    df = pd.DataFrame(
        [
            {"hole_id": "DH_01", "from": 0.0, "to": 5.0, "fault_zone": "Host", "grade": 1.0},
            {"hole_id": "DH_01", "from": 5.0, "to": 7.0, "fault_zone": "Fault", "grade": 8.0},
            {"hole_id": "DH_01", "from": 7.0, "to": 9.0, "fault_zone": "Host", "grade": 2.0},
            {"hole_id": "DH_02", "from": 0.0, "to": 4.0, "fault_zone": "Host", "grade": 3.0},
            {"hole_id": "DH_02", "from": 4.0, "to": 8.0, "fault_zone": "Host", "grade": 4.0},
        ]
    )

    default_result = analyze_drillhole_limit_behavior(
        df,
        hole_id_col="hole_id",
        from_col="from",
        to_col="to",
        limit_fields=["fault_zone"],
        numeric_fields=["grade"],
        boundary_categories_by_field={"fault_zone": ["Fault"]},
    )
    merged_default = default_result["merged_intervals"]
    assert set(merged_default["hole_id"].unique()) == {"DH_01"}
    quality_default = dict(zip(default_result["quality_summary"]["metric"], default_result["quality_summary"]["value"]))
    assert quality_default["holes_skipped_without_boundary_fault_zone"] == 1

    inclusive_result = analyze_drillhole_limit_behavior(
        df,
        hole_id_col="hole_id",
        from_col="from",
        to_col="to",
        limit_fields=["fault_zone"],
        numeric_fields=["grade"],
        boundary_categories_by_field={"fault_zone": ["Fault"]},
        skip_holes_without_boundary=False,
    )
    merged_inclusive = inclusive_result["merged_intervals"]
    assert set(merged_inclusive["hole_id"].unique()) == {"DH_01", "DH_02"}


def test_analyze_drillhole_limit_behavior_reports_progress():
    df = _build_behavior_frame()
    progress_events: list[tuple[int, int, str]] = []

    analyze_drillhole_limit_behavior(
        df,
        hole_id_col="hole_id",
        from_col="from",
        to_col="to",
        limit_fields=["fault_zone"],
        numeric_fields=["grade", "density"],
        boundary_categories_by_field={"fault_zone": ["Fault"]},
        progress_callback=lambda current, total, message: progress_events.append((current, total, message)),
    )

    assert progress_events
    assert progress_events[-1][0] == progress_events[-1][1]
    assert progress_events[-1][2] == "Analysis complete"