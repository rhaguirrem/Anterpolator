import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import anterpolator3DViewer as viewer_module

from anterpolator3DViewer import (
    _build_export_blocks_dataframe,
    FilterDataSource,
    apply_blank_sample_domain_behavior,
    compute_domain_sensitive_assignment_mask,
    create_blocks,
    detect_grid_rotation,
    ensure_sample_domains_for_domain_operations,
    export_block_volume_weighted_average,
    export_block_domain_sample_metrics,
    export_domain_interpolation_confidence_metrics,
    export_domained_samples_from_blocks,
    export_samples_with_block_values_from_blocks,
    load_block_domain_catalog,
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


def test_streaming_metadata_uses_single_full_pass(monkeypatch):
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

        progress_labels = []

        def _fake_iterate_csv_with_progress(path, progress_label, progress_callback=None, **read_csv_kwargs):
            progress_labels.append(progress_label)
            yield from pd.read_csv(path, **read_csv_kwargs)

        monkeypatch.setattr(viewer_module, "iterate_csv_with_progress", _fake_iterate_csv_with_progress)

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

        assert progress_labels.count("Reading grid file (rotation sample)") == 1
        assert progress_labels.count("Reading grid file (bounds + domain mapping)") == 1
        assert "Reading grid file (domain mapping)" not in progress_labels
        assert metadata["domain_mapping"] == {
            (0, 0, 0): "A",
            (1, 0, 0): "B",
        }


def test_create_blocks_rotates_samples_before_streaming_assignment(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        pd.DataFrame(
            [
                {"x": 0.0, "y": 0.0, "z": 0.0, "dom": "A"},
            ]
        ).to_csv(blocks_path, index=False)

        monkeypatch.setattr(viewer_module, "LARGE_BLOCK_FILE_THRESHOLD", 1)

        def _fake_load_large_blocks_metadata(*args, **kwargs):
            return {
                "all_min_bounds": np.array([0.0, 0.0, 0.0], dtype=float),
                "all_max_bounds": np.array([10.0, 10.0, 10.0], dtype=float),
                "unified_dims": np.array([10.0, 10.0, 10.0], dtype=float),
                "domain_mapping": {(0, 0, 0): "A"},
                "subblock_counts": {(0, 0, 0): 1},
                "mixed_domain_blocks": {},
                "rotation_matrix": np.eye(3, dtype=float),
                "rotation_center": np.array([100.0, 100.0, 100.0], dtype=float),
                "is_rotated": True,
            }

        class _FakeCell:
            def __init__(self):
                self.n_cells = 1
                self.cell_data = {}

        class _FakeMultiBlock(list):
            pass

        class _FakeInterpolator:
            def initialize_blocks(self, *args, **kwargs):
                self.initialized = True

            def create_ants(self):
                return None

        monkeypatch.setattr(viewer_module, "load_large_blocks_metadata", _fake_load_large_blocks_metadata)
        monkeypatch.setattr(
            viewer_module,
            "_require_pyvista",
            lambda: type("_FakePv", (), {"Box": staticmethod(lambda **kwargs: _FakeCell()), "MultiBlock": _FakeMultiBlock})(),
        )
        monkeypatch.setattr(viewer_module, "create_interpolator", lambda config, **kwargs: _FakeInterpolator())

        blocks = create_blocks(
            points=np.array([[100.0, 100.0, 100.0]], dtype=float),
            values=np.array([7.0], dtype=float),
            block_size=(10.0, 10.0, 10.0),
            blocks_file=blocks_path,
            blocks_delimiter=",",
            blocks_header_line=1,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            config={
                "algorithm": "ant_colony",
                "process_domains_sequentially": False,
                "expand_interpolation_exports_to_subblocks": True,
            },
            sample_domains=np.array(["A"], dtype=object),
        )

        assert blocks._sample_blocks == {(0, 0, 0): 7.0}
        assert blocks._sample_assignment_data["assigned_mask"].tolist() == [True]
        assert blocks._sample_assignment_data["block_indices"].tolist() == [[0, 0, 0]]


def test_load_block_domain_catalog_applies_filters_on_small_files():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")

        with open(blocks_path, "w", encoding="utf-8") as handle:
            handle.write("metadata line\n")
            handle.write("x,y,z,dom,lith\n")
            handle.write("5,5,5,A,keep\n")
            handle.write("15,5,5,B,drop\n")
            handle.write("25,5,5,C,keep\n")

        domains = load_block_domain_catalog(
            blocks_path,
            ",",
            2,
            "dom",
            block_filters=[{"field": "lith", "type": "categorical", "values": ["keep"]}],
        )

        assert domains == ["A", "C"]


def test_load_block_domain_catalog_uses_path_chunk_reader_for_large_filtered_files(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")

        with open(blocks_path, "w", encoding="utf-8") as handle:
            handle.write("metadata line\n")
            handle.write("x,y,z,dom,lith\n")
            handle.write("5,5,5,A,keep\n")
            handle.write("15,5,5,B,drop\n")
            handle.write("25,5,5,C,keep\n")

        monkeypatch.setattr(viewer_module, "LARGE_BLOCK_FILE_THRESHOLD", 1)

        def _fail_if_called(*args, **kwargs):
            raise AssertionError("iterate_csv_with_progress should not be used for large filtered domain catalog loads")

        monkeypatch.setattr(viewer_module, "iterate_csv_with_progress", _fail_if_called)

        progress_events = []

        def _progress_callback(value, maximum, message):
            progress_events.append((int(value), int(maximum), str(message)))

        domains = load_block_domain_catalog(
            blocks_path,
            ",",
            2,
            "dom",
            block_filters=[{"field": "lith", "type": "categorical", "values": ["keep"]}],
            progress_callback=_progress_callback,
        )

        assert domains == ["A", "C"]
        assert progress_events
        assert progress_events[-1][0] == progress_events[-1][1]
        assert progress_events[-1][2] == "Reading domain column"


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


def test_export_samples_with_block_values_from_blocks_transfers_multiple_columns():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "samples_with_blocks.csv")

        block_rows = []
        for x, grade, lith in [
            (2.5, 1.0, "ore"),
            (7.5, 3.0, "ore"),
            (12.5, 5.0, "waste"),
            (17.5, 7.0, "waste"),
        ]:
            for y in (2.5, 7.5):
                for z in (2.5, 7.5):
                    block_rows.append({"x": x, "y": y, "z": z, "grade": grade, "lith": lith})

        pd.DataFrame(block_rows).to_csv(blocks_path, index=False)
        pd.DataFrame(
            [
                {"sx": 5.0, "sy": 5.0, "sz": 5.0},
                {"sx": 15.0, "sy": 5.0, "sz": 5.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_samples_with_block_values_from_blocks(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_value_cols=["grade", "lith"],
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)

        assert result["matched_samples"] == 2
        assert result["transferred_columns"] == ["grade", "lith"]
        assert result["column_modes"] == {"grade": "numeric", "lith": "categorical"}
        assert output_df.loc[0, "grade"] == pytest.approx(2.0)
        assert output_df.loc[0, "lith"] == "ore"
        assert output_df.loc[1, "grade"] == pytest.approx(6.0)
        assert output_df.loc[1, "lith"] == "waste"


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


def test_export_block_domain_sample_metrics_can_emit_distance_band_counts_for_explicit_domains():
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
                {"sx": 40.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A"},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B"},
                {"sx": 31.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B"},
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
            distance_count_step=5.0,
            distance_count_max_factor=2,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=None,
        )

        output_df = pd.read_csv(output_path)
        expected_columns = [
            "dom_Sample_Count_0_5",
            "dom_Sample_Count_5_10",
            "dom_Sample_Count_GE_10",
        ]

        assert result["distance_count_columns"] == expected_columns
        assert result["distance_count_step"] == 5.0
        assert result["distance_count_max_factor"] == 2

        assert output_df.loc[0, expected_columns].tolist() == [1, 0, 2]
        assert output_df.loc[1, expected_columns].tolist() == [1, 0, 2]
        assert output_df.loc[2, expected_columns].tolist() == [1, 1, 0]


def test_export_block_domain_sample_metrics_can_emit_distance_band_counts_for_inferred_domains():
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
                {"sx": 2.0, "sy": 5.0, "sz": 5.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0},
                {"sx": 28.0, "sy": 5.0, "sz": 5.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            distance_count_step=5.0,
            distance_count_max_factor=2,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)
        expected_columns = [
            "dom_Sample_Count_0_5",
            "dom_Sample_Count_5_10",
            "dom_Sample_Count_GE_10",
        ]

        assert result["matched_samples"] == 4
        assert output_df.loc[0, expected_columns].tolist() == [1, 0, 1]
        assert output_df.loc[1, expected_columns].tolist() == [1, 0, 1]
        assert output_df.loc[2, expected_columns].tolist() == [2, 0, 0]


def test_export_domain_interpolation_confidence_metrics_summarizes_per_domain_distances():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "domain_confidence.csv")

        pd.DataFrame(
            [
                {"x": 0.0, "y": 0.0, "z": 0.0, "dom": "A"},
                {"x": 10.0, "y": 0.0, "z": 0.0, "dom": "A"},
                {"x": 100.0, "y": 0.0, "z": 0.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 0.0, "sy": 0.0, "sz": 0.0, "sample_dom": "A"},
                {"sx": 10.0, "sy": 0.0, "sz": 0.0, "sample_dom": "A"},
                {"sx": 100.0, "sy": 0.0, "sz": 0.0, "sample_dom": "B"},
            ]
        ).to_csv(samples_path, index=False)

        result = export_domain_interpolation_confidence_metrics(
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
        )

        output_df = pd.read_csv(output_path)

        assert result["domain_count"] == 2
        assert result["matched_samples"] == 3
        assert result["processed_blocks"] == 3

        domain_a = output_df.loc[output_df["Domain"] == "A"].iloc[0]
        assert domain_a["Source_Sample_Count"] == 2
        assert domain_a["Domain_Block_Count"] == 2
        assert domain_a["Avg_Source_Sample_Distance"] == pytest.approx(10.0)
        assert domain_a["Avg_Block_To_Source_Sample_Distance"] == pytest.approx(5.0)
        assert domain_a["Avg_Domain_Block_Distance"] == pytest.approx(10.0)
        assert domain_a["Sample_To_Block_Distance_Ratio"] == pytest.approx(1.0)

        domain_b = output_df.loc[output_df["Domain"] == "B"].iloc[0]
        assert domain_b["Source_Sample_Count"] == 1
        assert domain_b["Domain_Block_Count"] == 1
        assert pd.isna(domain_b["Avg_Source_Sample_Distance"])
        assert domain_b["Avg_Block_To_Source_Sample_Distance"] == pytest.approx(0.0)
        assert pd.isna(domain_b["Avg_Domain_Block_Distance"])
        assert pd.isna(domain_b["Sample_To_Block_Distance_Ratio"])


def test_export_domain_interpolation_confidence_metrics_includes_axiswise_distances():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "domain_confidence_axes.csv")

        pd.DataFrame(
            [
                {"x": 0.0, "y": 0.0, "z": 0.0, "dom": "A"},
                {"x": 10.0, "y": 20.0, "z": 30.0, "dom": "A"},
                {"x": 50.0, "y": 50.0, "z": 50.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 0.0, "sy": 0.0, "sz": 0.0, "sample_dom": "A"},
                {"sx": 10.0, "sy": 20.0, "sz": 30.0, "sample_dom": "A"},
                {"sx": 50.0, "sy": 50.0, "sz": 50.0, "sample_dom": "B"},
            ]
        ).to_csv(samples_path, index=False)

        export_domain_interpolation_confidence_metrics(
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
        )

        output_df = pd.read_csv(output_path)

        domain_a = output_df.loc[output_df["Domain"] == "A"].iloc[0]
        assert domain_a["Avg_Source_Sample_Distance_X"] == pytest.approx(10.0)
        assert domain_a["Avg_Source_Sample_Distance_Y"] == pytest.approx(20.0)
        assert domain_a["Avg_Source_Sample_Distance_Z"] == pytest.approx(30.0)
        assert domain_a["Avg_Block_To_Source_Sample_Distance_X"] == pytest.approx(5.0)
        assert domain_a["Avg_Block_To_Source_Sample_Distance_Y"] == pytest.approx(10.0)
        assert domain_a["Avg_Block_To_Source_Sample_Distance_Z"] == pytest.approx(15.0)
        assert domain_a["Avg_Domain_Block_Distance_X"] == pytest.approx(10.0)
        assert domain_a["Avg_Domain_Block_Distance_Y"] == pytest.approx(20.0)
        assert domain_a["Avg_Domain_Block_Distance_Z"] == pytest.approx(30.0)

        domain_b = output_df.loc[output_df["Domain"] == "B"].iloc[0]
        assert pd.isna(domain_b["Avg_Source_Sample_Distance_X"])
        assert domain_b["Avg_Block_To_Source_Sample_Distance_X"] == pytest.approx(0.0)
        assert pd.isna(domain_b["Avg_Domain_Block_Distance_X"])


def test_export_block_volume_weighted_average_infers_subblock_volumes_from_centroids():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        output_path = os.path.join(td, "blocks_weighted.csv")

        pd.DataFrame(
            [
                {"x": 2.5, "y": 5.0, "z": 5.0, "grade": 2.0},
                {"x": 7.5, "y": 5.0, "z": 5.0, "grade": 4.0},
                {"x": 15.0, "y": 5.0, "z": 5.0, "grade": 10.0},
            ]
        ).to_csv(blocks_path, index=False)

        result = export_block_volume_weighted_average(
            blocks_path,
            value_col="grade",
            output_file=output_path,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)

        assert result["processed_rows"] == 3
        assert result["exported_rows"] == 1
        assert result["rows_with_numeric_value"] == 3
        assert result["total_volume"] == 2000.0
        assert result["weighted_sum"] == 13000.0
        assert result["weighted_average"] == 6.5
        assert output_df.loc[0, "Total_Volume"] == 2000.0
        assert output_df.loc[0, "Weighted_Sum"] == 13000.0
        assert output_df.loc[0, "Weighted_Average"] == 6.5


def test_build_export_blocks_dataframe_can_expand_base_blocks_back_to_source_subblocks():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        pd.DataFrame(
            [
                {"x": 2.5, "y": 5.0, "z": 5.0},
                {"x": 7.5, "y": 5.0, "z": 5.0},
                {"x": 12.5, "y": 5.0, "z": 5.0},
                {"x": 17.5, "y": 5.0, "z": 5.0},
            ]
        ).to_csv(blocks_path, index=False)

        block_rows = [
            {
                "_Grid_Index": (0, 0, 0),
                "x": 0.0,
                "y": 5.0,
                "z": 5.0,
                "Value": 10.0,
                "Domain": "A",
                "Source": "First Pass",
            },
            {
                "_Grid_Index": (1, 0, 0),
                "x": 10.0,
                "y": 5.0,
                "z": 5.0,
                "Value": 20.0,
                "Domain": "B",
                "Source": "First Pass",
            },
        ]

        class FakeBlocks:
            pass

        fake_blocks = FakeBlocks()
        fake_blocks._block_info = {
            "min_bounds": [0.0, 0.0, 0.0],
            "block_size": [10.0, 10.0, 10.0],
            "rotation_matrix": None,
            "rotation_center": None,
            "source_blocks_file": blocks_path,
            "source_blocks_delimiter": ",",
            "source_blocks_header_line": 1,
            "source_block_x_col": "x",
            "source_block_y_col": "y",
            "source_block_z_col": "z",
            "source_block_filters": [],
            "expand_interpolation_exports_to_subblocks": True,
        }

        output_df = _build_export_blocks_dataframe(fake_blocks, block_rows)

        assert len(output_df) == 4
        assert output_df["x"].tolist() == [2.5, 7.5, 12.5, 17.5]
        assert output_df["Value"].tolist() == [10.0, 10.0, 20.0, 20.0]
        assert output_df["Domain"].tolist() == ["A", "A", "B", "B"]


def test_build_export_blocks_dataframe_keeps_empty_interpolation_export_empty():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        pd.DataFrame(
            [
                {"x": 2.5, "y": 5.0, "z": 5.0},
                {"x": 7.5, "y": 5.0, "z": 5.0},
            ]
        ).to_csv(blocks_path, index=False)

        class FakeBlocks:
            pass

        fake_blocks = FakeBlocks()
        fake_blocks._block_info = {
            "min_bounds": [0.0, 0.0, 0.0],
            "block_size": [10.0, 10.0, 10.0],
            "rotation_matrix": None,
            "rotation_center": None,
            "source_blocks_file": blocks_path,
            "source_blocks_delimiter": ",",
            "source_blocks_header_line": 1,
            "source_block_x_col": "x",
            "source_block_y_col": "y",
            "source_block_z_col": "z",
            "source_block_filters": [],
            "expand_interpolation_exports_to_subblocks": True,
        }

        output_df = _build_export_blocks_dataframe(fake_blocks, [])

        assert output_df.empty
        assert list(output_df.columns) == []


def test_export_block_volume_weighted_average_can_compute_domain_sensitive_summaries():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        output_path = os.path.join(td, "blocks_weighted.csv")

        pd.DataFrame(
            [
                {"x": 2.5, "y": 5.0, "z": 5.0, "grade": 2.0, "dom": "A"},
                {"x": 7.5, "y": 5.0, "z": 5.0, "grade": 4.0, "dom": "A"},
                {"x": 12.5, "y": 5.0, "z": 5.0, "grade": 10.0, "dom": "B"},
                {"x": 17.5, "y": 5.0, "z": 5.0, "grade": 14.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        result = export_block_volume_weighted_average(
            blocks_path,
            value_col="grade",
            output_file=output_path,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)

        assert result["domain_column"] == "dom"
        assert result["exported_rows"] == 2
        assert result["domain_summaries"]["A"]["total_volume"] == 1000.0
        assert result["domain_summaries"]["A"]["weighted_average"] == 3.0
        assert result["domain_summaries"]["B"]["total_volume"] == 1000.0
        assert result["domain_summaries"]["B"]["weighted_average"] == 12.0
        assert output_df["dom"].tolist() == ["A", "B"]
        assert output_df["Weighted_Average"].tolist() == [3.0, 12.0]
        assert output_df["Total_Volume"].tolist() == [1000.0, 1000.0]


def test_export_block_volume_weighted_average_can_use_custom_weight_column():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        output_path = os.path.join(td, "blocks_weighted.csv")

        pd.DataFrame(
            [
                {"grade": 2.0, "tonnage": 1.0},
                {"grade": 4.0, "tonnage": 3.0},
                {"grade": 10.0, "tonnage": 6.0},
            ]
        ).to_csv(blocks_path, index=False)

        result = export_block_volume_weighted_average(
            blocks_path,
            value_col="grade",
            weight_col="tonnage",
            output_file=output_path,
        )

        output_df = pd.read_csv(output_path)

        assert result["weight_column"] == "tonnage"
        assert result["processed_rows"] == 3
        assert result["rows_with_numeric_value"] == 3
        assert result["total_weight"] == 10.0
        assert result["weighted_sum"] == 74.0
        assert result["weighted_average"] == 7.4
        assert pd.isna(result["total_volume"])
        assert output_df.loc[0, "Weight_Column"] == "tonnage"
        assert output_df.loc[0, "Total_Weight"] == 10.0
        assert output_df.loc[0, "Weighted_Average"] == 7.4
        assert pd.isna(output_df.loc[0, "Total_Volume"])


def test_export_block_volume_weighted_average_applies_block_filters():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        output_path = os.path.join(td, "blocks_weighted.csv")

        pd.DataFrame(
            [
                {"x": 2.5, "y": 5.0, "z": 5.0, "grade": 2.0, "dom": "A"},
                {"x": 7.5, "y": 5.0, "z": 5.0, "grade": 4.0, "dom": "A"},
                {"x": 12.5, "y": 5.0, "z": 5.0, "grade": 10.0, "dom": "B"},
                {"x": 17.5, "y": 5.0, "z": 5.0, "grade": 14.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        result = export_block_volume_weighted_average(
            blocks_path,
            value_col="grade",
            output_file=output_path,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_size=(10, 10, 10),
            block_filters=[
                {"field": "dom", "type": "categorical", "values": ["A"]},
            ],
        )

        output_df = pd.read_csv(output_path)

        assert result["input_blocks"] == 4
        assert result["filtered_blocks"] == 2
        assert result["filters_applied"][0]["summary"] == "dom in [A]"
        assert result["total_volume"] == 1000.0
        assert result["weighted_average"] == 3.0
        assert output_df.loc[0, "Weighted_Average"] == 3.0
        assert output_df.loc[0, "Total_Volume"] == 1000.0


def test_filter_data_source_can_read_csv_columns_on_demand():
    with tempfile.TemporaryDirectory() as td:
        csv_path = os.path.join(td, "filters.csv")
        pd.DataFrame(
            [
                {"dom": "A", "grade": 2.0, "weight": 1.0},
                {"dom": "B", "grade": 5.0, "weight": 3.0},
                {"dom": "A", "grade": 8.0, "weight": 6.0},
            ]
        ).to_csv(csv_path, index=False)

        source = FilterDataSource(csv_path)

        assert source.has_field("dom") is True
        values, total_values, truncated = source.get_categorical_values("dom")
        assert values == ["A", "B"]
        assert total_values == 2
        assert truncated is False
        assert source.get_numeric_range("grade") == (2.0, 8.0)


def test_filter_data_source_supports_in_memory_dataframes():
    source = FilterDataSource(
        pd.DataFrame(
            [
                {"category": "X", "value": 1.0},
                {"category": "Y", "value": 4.0},
                {"category": "X", "value": 9.0},
            ]
        )
    )

    assert source.columns == ["category", "value"]
    assert source.get_numeric_range("value") == (1.0, 9.0)


def test_detect_grid_rotation_accepts_integer_points():
    points = [
        [0, 0, 0],
        [10, 0, 0],
        [0, 10, 0],
        [10, 10, 0],
        [0, 0, 10],
        [10, 0, 10],
        [0, 10, 10],
        [10, 10, 10],
        [20, 0, 0],
        [20, 10, 0],
    ]

    rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(points, block_size_hint=(10, 10, 10))

    assert rotation_matrix.shape == (3, 3)
    assert rotation_center.shape == (3,)
    assert bool(is_rotated) is False


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