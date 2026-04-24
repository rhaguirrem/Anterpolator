import os
import sys
import tempfile

import pandas as pd

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from anterpolator3DViewer import (
    apply_blank_sample_domain_behavior,
    compute_domain_sensitive_assignment_mask,
    ensure_sample_domains_for_domain_operations,
    export_block_domain_sample_metrics,
    export_domained_samples_from_blocks,
    load_large_blocks_metadata,
    normalize_selected_sample_domain_column,
)


def _capture_progress_events():
    events = []

    def callback(value, maximum, message):
        events.append((int(value), int(maximum), str(message)))

    return events, callback


def test_streaming_metadata_collapses_subblocks_to_base_blocks():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        pd.DataFrame(
            [
                {"x": 2.5, "y": 2.5, "z": 2.5, "dom": "A"},
                {"x": 7.5, "y": 2.5, "z": 2.5, "dom": "A"},
                {"x": 12.5, "y": 2.5, "z": 2.5, "dom": "B"},
                {"x": 17.5, "y": 2.5, "z": 2.5, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        metadata = load_large_blocks_metadata(
            blocks_path,
            ",",
            1,
            (10, 10, 10),
            None,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            config=None,
        )

        assert metadata["domain_mapping"] == {
            (0, 0, 0): "A",
            (1, 0, 0): "B",
        }
        assert metadata["subblock_counts"] == {
            (0, 0, 0): 2,
            (1, 0, 0): 2,
        }


def test_export_domained_samples_uses_aggregated_subblock_domains():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "domained.csv")

        block_rows = []

        for x in (2.5, 7.5):
            for y in (2.5, 7.5):
                for z in (2.5, 7.5):
                    block_rows.append({"x": x, "y": y, "z": z, "dom": "X"})

        for x in (12.5, 17.5):
            for y in (2.5, 7.5):
                for z in (2.5, 7.5):
                    domain = "B" if x == 12.5 else "A"
                    block_rows.append({"x": x, "y": y, "z": z, "dom": domain})

        pd.DataFrame(block_rows).to_csv(blocks_path, index=False)
        pd.DataFrame([
            {"sx": 15.0, "sy": 5.0, "sz": 5.0},
        ]).to_csv(samples_path, index=False)

        result = export_domained_samples_from_blocks(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)

        assert result["matched_samples"] == 1
        assert result["unmatched_samples"] == 0
        assert output_df.loc[0, "dom"] == "A"


def test_export_block_domain_sample_metrics_applies_sample_filters():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "lith": "keep", "grade": 1.0},
                {"sx": 8.0, "sy": 5.0, "sz": 5.0, "lith": "drop", "grade": 5.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "lith": "keep", "grade": 3.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "lith": "keep", "grade": 2.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
            sample_filters=[
                {"field": "lith", "type": "categorical", "values": ["keep"]},
                {"field": "grade", "type": "numeric", "min": 2.0, "max": 3.0},
            ],
        )

        output_df = pd.read_csv(output_path)

        assert result["filtered_samples"] == 2
        assert result["matched_samples"] == 2
        assert len(result["filters_applied"]) == 2
        assert "sample_count_column" not in result
        assert "dom_Sample_Count" not in output_df.columns

        assert output_df.loc[0, "dom_NN_Distance"] == 13.0
        assert output_df.loc[0, "dom_Avg_Distance"] == 13.0

        assert output_df.loc[1, "dom_NN_Distance"] == 3.0
        assert output_df.loc[1, "dom_Avg_Distance"] == 3.0

        assert output_df.loc[2, "dom_NN_Distance"] == 0.0
        assert output_df.loc[2, "dom_Avg_Distance"] == 0.0


def test_export_block_domain_sample_metrics_can_emit_closest_sample_id_from_multiple_columns():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-01", "from_m": 10, "lith": "A"},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-02", "from_m": 20, "lith": "A"},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-03", "from_m": 30, "lith": "B"},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            sample_domain_col="lith",
            closest_sample_id_cols=["hole", "from_m"],
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)

        assert result["closest_sample_id_column"] == "dom_Closest_Sample_ID"
        assert result["closest_sample_id_source_columns"] == ["hole", "from_m"]
        assert output_df.loc[0, "dom_Closest_Sample_ID"] == "DDH-01 | 10"
        assert output_df.loc[1, "dom_Closest_Sample_ID"] == "DDH-02 | 20"
        assert output_df.loc[2, "dom_Closest_Sample_ID"] == "DDH-03 | 30"


def test_export_block_domain_sample_metrics_uses_explicit_sample_domains_without_block_size():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A"},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A"},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B"},
                {"sx": 30.0, "sy": 5.0, "sz": 5.0, "sample_dom": "C"},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            sample_domain_col="sample_dom",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=None,
        )

        output_df = pd.read_csv(output_path)

        assert result["matched_samples"] == 3
        assert result["unmatched_samples"] == 1
        assert output_df.loc[0, "dom_NN_Distance"] == 3.0
        assert output_df.loc[1, "dom_NN_Distance"] == 3.0
        assert output_df.loc[2, "dom_NN_Distance"] == 0.0


def test_normalize_selected_sample_domain_column_preserves_domain_when_value_column_matches():
    df = pd.DataFrame(
        [
            {"mid_x": 1.0, "mid_y": 2.0, "mid_z": 3.0, "DG": "A"},
            {"mid_x": 4.0, "mid_y": 5.0, "mid_z": 6.0, "DG": "B"},
        ]
    )

    normalized = normalize_selected_sample_domain_column(df, sample_domain_col="DG")

    assert "DG" in normalized.columns
    assert "Domain" in normalized.columns
    assert normalized["Domain"].tolist() == ["A", "B"]


def test_apply_blank_sample_domain_behavior_can_infer_from_blocks():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        df = pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "Domain": ""},
                {"x": 15.0, "y": 5.0, "z": 5.0, "Domain": "B"},
            ]
        )

        resolved_df, summary = apply_blank_sample_domain_behavior(
            df,
            blank_domain_behavior="infer_from_blocks",
            blocks_file=blocks_path,
            block_size=(10, 10, 10),
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
        )

        assert resolved_df["Domain"].tolist() == ["A", "B"]
        assert summary["blank_domains"] == 1
        assert summary["inferred_domains"] == 1
        assert summary["remaining_blank_domains"] == 0


def test_ensure_sample_domains_for_domain_operations_infers_all_when_no_domain_column_selected():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        df = pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0},
                {"x": 15.0, "y": 5.0, "z": 5.0},
            ]
        )

        resolved_df, summary = ensure_sample_domains_for_domain_operations(
            df,
            sample_domain_col=None,
            blank_domain_behavior="skip",
            blocks_file=blocks_path,
            block_size=(10, 10, 10),
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
        )

        assert resolved_df["Domain"].tolist() == ["A", "B"]
        assert summary["remaining_blank_domains"] == 0


def test_compute_domain_sensitive_assignment_mask_rejects_mismatched_domains():
    block_indices = [(0, 0, 0), (1, 0, 0), (1, 0, 0)]
    allowed_grid = {(0, 0, 0), (1, 0, 0)}
    domain_mapping = {
        (0, 0, 0): "A",
        (1, 0, 0): "B",
    }
    sample_domains = ["A", "A", "B"]

    assigned_mask = compute_domain_sensitive_assignment_mask(
        block_indices,
        allowed_grid,
        domain_mapping=domain_mapping,
        sample_domains=sample_domains,
    )

    assert assigned_mask.tolist() == [True, False, True]


def test_export_block_domain_sample_metrics_can_infer_blank_explicit_sample_domains():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A"},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "sample_dom": ""},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B"},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            sample_domain_col="sample_dom",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
            blank_sample_domain_behavior="infer_from_blocks",
        )

        assert result["matched_samples"] == 3
        assert result["unmatched_samples"] == 0


def test_export_domained_samples_reports_progress():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "domained.csv")
        events, progress_callback = _capture_progress_events()

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 5.0, "sy": 5.0, "sz": 5.0},
                {"sx": 15.0, "sy": 5.0, "sz": 5.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_domained_samples_from_blocks(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
            progress_callback=progress_callback,
        )

        assert result["matched_samples"] == 2
        assert events
        assert events[-1][0] == 100
        assert any("Reading sample file" in message for _, _, message in events)
        assert any("Assigning samples to block domains" in message for _, _, message in events)


def test_export_block_domain_sample_metrics_reports_progress():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")
        events, progress_callback = _capture_progress_events()

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
            progress_callback=progress_callback,
        )

        assert result["processed_blocks"] == 3
        assert events
        assert events[-1][0] == 100
        assert any("Reading sample file" in message for _, _, message in events)
        assert any("Reading blocks file" in message for _, _, message in events)
        assert any("Computing distance statistics" in message for _, _, message in events)