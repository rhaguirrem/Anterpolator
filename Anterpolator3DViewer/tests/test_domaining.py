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
    _restore_list_widget_selection,
    FilterDataSource,
    apply_blank_sample_domain_behavior,
    compute_domain_sensitive_assignment_mask,
    create_blocks,
    detect_grid_rotation,
    ensure_sample_domains_for_domain_operations,
    export_blocks_to_csv,
    export_block_volume_weighted_average,
    export_block_domain_sample_metrics,
    export_blocks_with_source_block_values,
    export_domain_interpolation_confidence_metrics,
    export_domained_samples_from_blocks,
    export_samples_with_block_values_from_blocks,
    load_block_domain_catalog,
    load_large_blocks_metadata,
    normalize_selected_sample_domain_column,
    parse_leapfrog_block_metadata,
    resolve_effective_csv_header_line,
    sync_csv_header_line_widget,
)


_QT_APP = None


def _get_qt_app():
    global _QT_APP
    app = viewer_module.QtWidgets.QApplication.instance()
    if app is None:
        _QT_APP = viewer_module.QtWidgets.QApplication([])
        app = _QT_APP
    return app


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


def test_streaming_metadata_warns_when_block_rows_do_not_match_parent_grid():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "leapfrog_blocks_outside_grid.csv")
        with open(blocks_path, "w", encoding="utf-8") as handle:
            handle.write("#   rotation type: Leapfrog\n")
            handle.write("#   azimuth: 0 degrees\n")
            handle.write("#   dip: 0 degrees\n")
            handle.write("#   pitch: 0 degrees\n")
            handle.write("#   parent block size: 10 10 10\n")
            handle.write("#   size in parent blocks: 1 1 1 = 1\n")
            handle.write("#   minimum corner: 0 0 0\n")
            handle.write("#   maximum corner: 10 10 10\n")
            handle.write("x,y,z,dom\n")
            handle.write("5,5,5,A\n")
            handle.write("15,5,5,B\n")

        with pytest.warns(RuntimeWarning) as warning_records:
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

        warning_messages = [str(record.message) for record in warning_records]
        assert any("map outside metadata parent grid" in message for message in warning_messages)
        assert any("extend above metadata maximum corner" in message for message in warning_messages)
        assert metadata["domain_mapping"][(0, 0, 0)] == "A"


def test_leapfrog_metadata_parser_tolerates_missing_or_unrecognized_comments():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        with open(blocks_path, "w", encoding="utf-8") as handle:
            handle.write("# arbitrary export note\n")
            handle.write("# another comment without a known key\n")
            handle.write("x,y,z,dom\n")
            handle.write("1,2,3,A\n")

        metadata = parse_leapfrog_block_metadata(blocks_path)

        assert metadata["raw_lines"] == ["arbitrary export note", "another comment without a known key"]
        assert "parent_block_size" not in metadata
        assert resolve_effective_csv_header_line(blocks_path, 1) == 3

    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "plain.csv")
        with open(blocks_path, "w", encoding="utf-8") as handle:
            handle.write("x,y,z,dom\n")
            handle.write("1,2,3,A\n")

        assert parse_leapfrog_block_metadata(blocks_path) == {}
        assert resolve_effective_csv_header_line(blocks_path, 1) == 1


def test_sync_csv_header_line_widget_updates_for_leapfrog_metadata():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "leapfrog_blocks.csv")
        with open(blocks_path, "w", encoding="utf-8") as handle:
            handle.write("# arbitrary export note\n")
            handle.write("# another comment without a known key\n")
            handle.write("x,y,z,dom\n")
            handle.write("1,2,3,A\n")

        app = _get_qt_app()
        spin_box = viewer_module.QtWidgets.QSpinBox()
        spin_box.setRange(1, 1_000_000)
        spin_box.setValue(1)

        effective_line = sync_csv_header_line_widget(spin_box, blocks_path)
        app.processEvents()

        assert effective_line == 3
        assert spin_box.value() == 3


def test_streaming_metadata_uses_leapfrog_minimum_corner_for_parent_indexing():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "leapfrog_blocks.csv")
        with open(blocks_path, "w", encoding="utf-8") as handle:
            handle.write("# BM DG May2026 v06 20260410.csv\n")
            handle.write("#   exported from Leapfrog Geo\n")
            handle.write("#   rotation type: Leapfrog\n")
            handle.write("#   azimuth: 0 degrees (rotate clockwise around the Z axis when looking down)\n")
            handle.write("#   dip: 0 degrees\n")
            handle.write("#   pitch: 0 degrees\n")
            handle.write("#   parent block size: 25 25 15\n")
            handle.write("#   size in parent blocks: 2 1 1 = 2\n")
            handle.write("#   minimum parent centroid: 11547.5 99032.5 707.5\n")
            handle.write("#   maximum parent centroid: 11572.5 99032.5 707.5\n")
            handle.write("#   minimum corner: 11535 99020 700\n")
            handle.write("#   maximum corner: 11585 99045 715\n")
            handle.write("#   sub-blocks: octree 4 4 4\n")
            handle.write("x,y,z,dom\n")
            handle.write("11538.125,99023.125,701.875,A\n")
            handle.write("11556.875,99023.125,701.875,A\n")
            handle.write("11563.125,99023.125,701.875,B\n")
            handle.write("11581.875,99023.125,701.875,B\n")

        with pytest.warns(RuntimeWarning, match="configured block size"):
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

        assert metadata["unified_dims"].tolist() == [25.0, 25.0, 15.0]
        assert metadata["grid_index_origin"].tolist() == [11535.0, 99020.0, 700.0]
        assert metadata["source_blocks_header_line"] == 14
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


def test_apply_sample_filters_matches_decimal_string_equivalents():
    df = pd.DataFrame(
        {
            "Valido 2": [1.0, 0.0, 1, "1.0", "01", "A"],
            "value": [10, 20, 30, 40, 50, 60],
        }
    )

    filtered_df, applied_filters = viewer_module.apply_sample_filters(
        df,
        sample_filters=[{"field": "Valido 2", "type": "categorical", "values": [1]}],
    )

    assert filtered_df["Valido 2"].tolist() == [1.0, 1, "1.0"]
    assert applied_filters[0]["summary"] == "Valido 2 in [1]"


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


def test_export_blocks_with_source_values_preserves_target_rows_and_exact_matches():
    with tempfile.TemporaryDirectory() as td:
        source_path = os.path.join(td, "source.csv")
        target_path = os.path.join(td, "target.csv")
        output_path = os.path.join(td, "target_enriched.csv")
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "grade": 1.5, "lith": "ore"},
            {"x": 15.0, "y": 5.0, "z": 5.0, "grade": 3.5, "lith": "waste"},
        ]).to_csv(source_path, index=False)
        pd.DataFrame([
            {"tx": 5.0, "ty": 5.0, "tz": 5.0, "target_id": "A"},
            {"tx": 15.0, "ty": 5.0, "tz": 5.0, "target_id": "B"},
        ]).to_csv(target_path, index=False)

        result = export_blocks_with_source_block_values(
            source_path,
            target_path,
            output_file=output_path,
            source_x_col="x", source_y_col="y", source_z_col="z",
            target_x_col="tx", target_y_col="ty", target_z_col="tz",
            source_value_cols=["grade", "lith"],
            source_block_size=(10, 10, 10),
            target_block_size=(10, 10, 10),
        )
        output_df = pd.read_csv(output_path)

        assert len(output_df) == 2
        assert output_df["target_id"].tolist() == ["A", "B"]
        assert output_df["grade"].tolist() == pytest.approx([1.5, 3.5])
        assert output_df["lith"].tolist() == ["ore", "waste"]
        assert result["overlap_matched_blocks"] == 2
        assert result["nearest_matched_blocks"] == 0


def test_export_blocks_with_source_values_uses_exact_grid_fast_path(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        source_path = os.path.join(td, "source.csv")
        target_path = os.path.join(td, "target.csv")
        output_path = os.path.join(td, "target_enriched.csv")
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "grade": 1.5},
            {"x": 15.0, "y": 5.0, "z": 5.0, "grade": 3.5},
        ]).to_csv(source_path, index=False)
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "target_id": "A"},
            {"x": 15.0, "y": 5.0, "z": 5.0, "target_id": "B"},
        ]).to_csv(target_path, index=False)

        def fail_geometry(*args, **kwargs):
            raise AssertionError("geometry path should not run for exact grid matches")

        monkeypatch.setattr(viewer_module, "_resolve_block_row_geometry", fail_geometry)

        result = viewer_module.export_blocks_with_source_block_values(
            source_path,
            target_path,
            output_file=output_path,
            source_x_col="x", source_y_col="y", source_z_col="z",
            target_x_col="x", target_y_col="y", target_z_col="z",
            source_value_cols=["grade"],
            source_block_size=(10, 10, 10),
            target_block_size=(10, 10, 10),
        )
        output_df = pd.read_csv(output_path)

        assert output_df["grade"].tolist() == pytest.approx([1.5, 3.5])
        assert result["source_geometry_mode"] == "exact-grid"
        assert result["target_geometry_mode"] == "exact-grid"


def test_export_blocks_with_source_values_exact_prematch_skips_matched_targets(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        source_path = os.path.join(td, "source.csv")
        target_path = os.path.join(td, "target.csv")
        output_path = os.path.join(td, "target_enriched.csv")
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "grade": 1.5},
            {"x": 15.0, "y": 5.0, "z": 5.0, "grade": 3.5},
        ]).to_csv(source_path, index=False)
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "target_id": "A"},
            {"x": 15.0, "y": 5.0, "z": 5.0, "target_id": "B"},
            {"x": 105.0, "y": 5.0, "z": 5.0, "target_id": "C"},
        ]).to_csv(target_path, index=False)

        original_resolve = viewer_module._resolve_block_row_geometry
        target_geometry_sizes = []

        def wrapped_resolve(df, coordinate_columns, base_block_size, size_columns=None,
                            progress_callback=None, progress_label='Resolving block geometry'):
            if 'target_id' in df.columns:
                target_geometry_sizes.append(len(df))
                assert len(df) == 1
            return original_resolve(df, coordinate_columns, base_block_size, size_columns, progress_callback, progress_label)

        monkeypatch.setattr(viewer_module, "_resolve_block_row_geometry", wrapped_resolve)

        result = viewer_module.export_blocks_with_source_block_values(
            source_path,
            target_path,
            output_file=output_path,
            source_x_col="x", source_y_col="y", source_z_col="z",
            target_x_col="x", target_y_col="y", target_z_col="z",
            source_value_cols=["grade"],
            source_block_size=(10, 10, 10),
            target_block_size=(10, 10, 10),
        )
        output_df = pd.read_csv(output_path)

        assert target_geometry_sizes == [1]
        assert output_df["grade"].tolist() == pytest.approx([1.5, 3.5, 3.5])
        assert result["overlap_matched_blocks"] == 2
        assert result["nearest_matched_blocks"] == 1
        assert result["target_geometry_mode"].startswith("exact-prematch + ")


def test_export_blocks_with_source_values_streams_csv_target(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        source_path = os.path.join(td, "source.csv")
        target_path = os.path.join(td, "target.csv")
        output_path = os.path.join(td, "target_enriched.csv")
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "grade": 1.5},
            {"x": 15.0, "y": 5.0, "z": 5.0, "grade": 3.5},
        ]).to_csv(source_path, index=False)
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "target_id": "A"},
            {"x": 15.0, "y": 5.0, "z": 5.0, "target_id": "B"},
        ]).to_csv(target_path, index=False)

        original_loader = viewer_module.load_full_blocks_dataframe

        def wrapped_loader(path, *args, **kwargs):
            if os.path.abspath(path) == os.path.abspath(target_path):
                raise AssertionError("CSV target should be streamed, not loaded wholesale")
            return original_loader(path, *args, **kwargs)

        monkeypatch.setattr(viewer_module, "load_full_blocks_dataframe", wrapped_loader)

        result = viewer_module.export_blocks_with_source_block_values(
            source_path,
            target_path,
            output_file=output_path,
            source_x_col="x", source_y_col="y", source_z_col="z",
            target_x_col="x", target_y_col="y", target_z_col="z",
            source_value_cols=["grade"],
            source_block_size=(10, 10, 10),
            target_block_size=(10, 10, 10),
        )
        output_df = pd.read_csv(output_path)

        assert output_df["grade"].tolist() == pytest.approx([1.5, 3.5])
        assert result["target_geometry_mode"] == "exact-grid"


def test_restore_list_widget_selection_restores_multi_select_items():
    _get_qt_app()
    widget = viewer_module.QtWidgets.QListWidget()
    widget.setSelectionMode(viewer_module.QtWidgets.QAbstractItemView.MultiSelection)
    for column in ["grade", "lith", "density"]:
        widget.addItem(viewer_module.QtWidgets.QListWidgetItem(column))

    _restore_list_widget_selection(widget, ["grade", "lith"])

    assert [item.text() for item in widget.selectedItems()] == ["grade", "lith"]


def test_export_blocks_with_source_values_uses_overlap_volume_for_subblocks():
    with tempfile.TemporaryDirectory() as td:
        source_path = os.path.join(td, "source_subblocks.csv")
        target_path = os.path.join(td, "target.csv")
        output_path = os.path.join(td, "target_enriched.csv")
        pd.DataFrame([
            {"x": 2.0, "y": 5.0, "z": 5.0, "dx": 4.0, "dy": 10.0, "dz": 10.0, "grade": 2.0, "lith": "A"},
            {"x": 7.0, "y": 5.0, "z": 5.0, "dx": 6.0, "dy": 10.0, "dz": 10.0, "grade": 6.0, "lith": "B"},
        ]).to_csv(source_path, index=False)
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "target_id": 1},
        ]).to_csv(target_path, index=False)

        result = export_blocks_with_source_block_values(
            source_path,
            target_path,
            output_file=output_path,
            source_x_col="x", source_y_col="y", source_z_col="z",
            target_x_col="x", target_y_col="y", target_z_col="z",
            source_value_cols=["grade", "lith"],
            source_block_size=(10, 10, 10),
            target_block_size=(10, 10, 10),
            source_size_cols=("dx", "dy", "dz"),
        )
        output_df = pd.read_csv(output_path)

        assert output_df.loc[0, "grade"] == pytest.approx(4.4)
        assert output_df.loc[0, "lith"] == "B"
        assert result["source_geometry_mode"] == "explicit-size-columns"
        assert result["overlap_matched_blocks"] == 1


def test_export_blocks_with_source_values_uses_nearest_fallback_without_creating_rows():
    with tempfile.TemporaryDirectory() as td:
        source_path = os.path.join(td, "source.csv")
        target_path = os.path.join(td, "target.csv")
        output_path = os.path.join(td, "target_enriched.csv")
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "grade": 8.0},
        ]).to_csv(source_path, index=False)
        pd.DataFrame([
            {"x": 105.0, "y": 5.0, "z": 5.0, "target_id": "far"},
        ]).to_csv(target_path, index=False)

        result = export_blocks_with_source_block_values(
            source_path,
            target_path,
            output_file=output_path,
            source_x_col="x", source_y_col="y", source_z_col="z",
            target_x_col="x", target_y_col="y", target_z_col="z",
            source_value_cols=["grade"],
            source_block_size=(10, 10, 10),
            target_block_size=(10, 10, 10),
        )
        output_df = pd.read_csv(output_path)

        assert len(output_df) == 1
        assert output_df.loc[0, "target_id"] == "far"
        assert output_df.loc[0, "grade"] == pytest.approx(8.0)
        assert result["overlap_matched_blocks"] == 0
        assert result["nearest_matched_blocks"] == 1


def test_export_blocks_with_source_values_respects_max_nearest_distance():
    with tempfile.TemporaryDirectory() as td:
        source_path = os.path.join(td, "source.csv")
        target_path = os.path.join(td, "target.csv")
        output_path = os.path.join(td, "target_enriched.csv")
        pd.DataFrame([
            {"x": 5.0, "y": 5.0, "z": 5.0, "grade": 8.0},
        ]).to_csv(source_path, index=False)
        pd.DataFrame([
            {"x": 105.0, "y": 5.0, "z": 5.0, "target_id": "far"},
        ]).to_csv(target_path, index=False)

        result = export_blocks_with_source_block_values(
            source_path,
            target_path,
            output_file=output_path,
            source_x_col="x", source_y_col="y", source_z_col="z",
            target_x_col="x", target_y_col="y", target_z_col="z",
            source_value_cols=["grade"],
            source_block_size=(10, 10, 10),
            target_block_size=(10, 10, 10),
            max_nearest_distance=25.0,
        )
        output_df = pd.read_csv(output_path)

        assert len(output_df) == 1
        assert pd.isna(output_df.loc[0, "grade"])
        assert result["overlap_matched_blocks"] == 0
        assert result["nearest_matched_blocks"] == 0
        assert result["unmatched_blocks"] == 1


def test_export_domained_samples_preserves_empty_columns():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "domained.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A"},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "B"},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 5.0, "sy": 5.0, "sz": 5.0, "empty_col": ""},
                {"sx": 15.0, "sy": 5.0, "sz": 5.0, "empty_col": ""},
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
        )

        output_df = pd.read_csv(output_path)

        assert result["matched_samples"] == 2
        assert "empty_col" in output_df.columns
        assert output_df["empty_col"].isna().all()


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
            sample_value_col="grade",
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
        nearest_distance_column = result["nearest_distance_column"]

        assert result["filtered_samples"] == 2
        assert result["matched_samples"] == 2
        assert len(result["filters_applied"]) == 2
        assert "sample_count_column" not in result
        assert "dom_Sample_Count" not in output_df.columns

        assert output_df.loc[0, nearest_distance_column] == 13.0
        assert output_df.loc[0, "dom_Avg_Distance"] == 13.0

        assert output_df.loc[1, nearest_distance_column] == 3.0
        assert output_df.loc[1, "dom_Avg_Distance"] == 3.0

        assert output_df.loc[2, nearest_distance_column] == 0.0
        assert output_df.loc[2, "dom_Avg_Distance"] == 0.0


def test_export_block_domain_sample_metrics_can_emit_closest_sample_id_from_multiple_columns():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 1.0},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 2.0},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "B", "block_grade": 3.0},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-01", "from_m": 10, "lith": "A", "assay": 1.5},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-02", "from_m": 20, "lith": "A", "assay": 2.5},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-03", "from_m": 30, "lith": "B", "assay": 3.5},
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
            sample_value_col="assay",
            closest_sample_id_cols=["hole", "from_m"],
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_value_col="block_grade",
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)

        assert result["closest_sample_id_column"] == "block_grade_dom_Closest_Sample_ID"
        assert result["closest_sample_id_source_columns"] == ["hole", "from_m"]
        assert output_df.loc[0, "block_grade_dom_Closest_Sample_ID"] == "DDH-01 | 10"
        assert output_df.loc[1, "block_grade_dom_Closest_Sample_ID"] == "DDH-02 | 20"
        assert output_df.loc[2, "block_grade_dom_Closest_Sample_ID"] == "DDH-03 | 30"


def test_export_block_domain_sample_metrics_can_disable_prefix_for_closest_sample_id_column():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 1.0},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-01", "lith": "A", "assay": 1.5},
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
            sample_value_col="assay",
            selected_metrics=["closest_sample_id"],
            closest_sample_id_cols=["hole"],
            use_block_value_prefix=False,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_value_col="block_grade",
        )

        output_df = pd.read_csv(output_path)

        assert result["closest_sample_id_column"] == "dom_Closest_Sample_ID"
        assert "dom_Closest_Sample_ID" in output_df.columns
        assert "block_grade_dom_Closest_Sample_ID" not in output_df.columns
        assert output_df.loc[0, "dom_Closest_Sample_ID"] == "DDH-01"


def test_export_block_domain_sample_metrics_closest_sample_id_ignores_nearer_samples_without_values():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 1.0},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 4.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-NULL", "lith": "A", "assay": np.nan},
                {"sx": 9.0, "sy": 5.0, "sz": 5.0, "hole": "DDH-VALID", "lith": "A", "assay": 9.0},
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
            sample_value_col="assay",
            selected_metrics=["closest_sample_id"],
            closest_sample_id_cols=["hole"],
            use_block_value_prefix=False,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_value_col="block_grade",
        )

        output_df = pd.read_csv(output_path)

        assert result["closest_sample_id_column"] == "dom_Closest_Sample_ID"
        assert output_df.loc[0, "dom_Closest_Sample_ID"] == "DDH-VALID"


def test_export_block_domain_sample_metrics_can_emit_nearest_sample_value_residual_metrics():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 5.0},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 9.0},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 11.0},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": 4.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": 10.0},
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
            sample_value_col="assay",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_value_col="block_grade",
        )

        output_df = pd.read_csv(output_path)
        nearest_sample_value_column = result["nearest_sample_value_column"]
        nearest_sample_residual_column = result["nearest_sample_residual_column"]
        nearest_sample_abs_residual_column = result["nearest_sample_abs_residual_column"]
        nearest_sample_group_block_count_column = result["nearest_sample_group_block_count_column"]
        nearest_sample_group_mean_residual_column = result["nearest_sample_group_mean_residual_column"]
        nearest_sample_group_rms_residual_column = result["nearest_sample_group_rms_residual_column"]
        nearest_sample_group_std_residual_column = result["nearest_sample_group_std_residual_column"]

        assert nearest_sample_value_column == "block_grade_Nearest_Sample_Value"
        assert nearest_sample_group_std_residual_column == "block_grade_Nearest_Sample_Group_StdResidual"

        assert output_df.loc[0, nearest_sample_value_column] == pytest.approx(4.0)
        assert output_df.loc[0, nearest_sample_residual_column] == pytest.approx(-1.0)
        assert output_df.loc[0, nearest_sample_abs_residual_column] == pytest.approx(1.0)
        assert output_df.loc[0, nearest_sample_group_block_count_column] == pytest.approx(1.0)
        assert output_df.loc[0, nearest_sample_group_mean_residual_column] == pytest.approx(-1.0)
        assert output_df.loc[0, nearest_sample_group_rms_residual_column] == pytest.approx(1.0)
        assert output_df.loc[0, nearest_sample_group_std_residual_column] == pytest.approx(-1.0)

        assert output_df.loc[1, nearest_sample_value_column] == pytest.approx(10.0)
        assert output_df.loc[1, nearest_sample_residual_column] == pytest.approx(1.0)
        assert output_df.loc[1, nearest_sample_abs_residual_column] == pytest.approx(1.0)
        assert output_df.loc[1, nearest_sample_group_block_count_column] == pytest.approx(2.0)
        assert output_df.loc[1, nearest_sample_group_mean_residual_column] == pytest.approx(0.0)
        assert output_df.loc[1, nearest_sample_group_rms_residual_column] == pytest.approx(1.0)
        assert output_df.loc[1, nearest_sample_group_std_residual_column] == pytest.approx(1.0)

        assert output_df.loc[2, nearest_sample_value_column] == pytest.approx(10.0)
        assert output_df.loc[2, nearest_sample_residual_column] == pytest.approx(-1.0)
        assert output_df.loc[2, nearest_sample_abs_residual_column] == pytest.approx(1.0)
        assert output_df.loc[2, nearest_sample_group_block_count_column] == pytest.approx(2.0)
        assert output_df.loc[2, nearest_sample_group_mean_residual_column] == pytest.approx(0.0)
        assert output_df.loc[2, nearest_sample_group_rms_residual_column] == pytest.approx(1.0)
        assert output_df.loc[2, nearest_sample_group_std_residual_column] == pytest.approx(-1.0)


def test_export_block_domain_sample_metrics_skips_closer_samples_without_values_for_value_metrics():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 5.0},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 4.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": np.nan, "hole": "DDH-NULL"},
                {"sx": 9.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": 9.0, "hole": "DDH-VALID"},
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
            sample_value_col="assay",
            selected_metrics=["nearest_distance", "closest_sample_id", "nearest_sample_value", "nearest_sample_residual"],
            closest_sample_id_cols=["hole"],
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_value_col="block_grade",
        )

        output_df = pd.read_csv(output_path)
        nearest_sample_value_column = result["nearest_sample_value_column"]
        nearest_sample_residual_column = result["nearest_sample_residual_column"]
        nearest_distance_column = result["nearest_distance_column"]

        assert nearest_distance_column == "block_grade_dom_NN_Distance"
        assert output_df.loc[0, nearest_distance_column] == pytest.approx(4.0)
        assert output_df.loc[0, "block_grade_dom_Closest_Sample_ID"] == "DDH-VALID"
        assert output_df.loc[0, nearest_sample_value_column] == pytest.approx(9.0)
        assert output_df.loc[0, nearest_sample_residual_column] == pytest.approx(4.0)
        assert nearest_sample_value_column == "block_grade_Nearest_Sample_Value"


def test_export_block_domain_sample_metrics_can_select_knn_average_without_exporting_unchecked_metrics():
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        samples_path = os.path.join(td, "samples.csv")
        output_path = os.path.join(td, "block_metrics.csv")

        pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 5.0},
                {"x": 15.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 7.0},
                {"x": 25.0, "y": 5.0, "z": 5.0, "dom": "A", "block_grade": 9.0},
            ]
        ).to_csv(blocks_path, index=False)

        pd.DataFrame(
            [
                {"sx": 0.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "hole": "DDH-01", "assay": 4.0},
                {"sx": 10.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "hole": "DDH-02", "assay": 6.0},
                {"sx": 30.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "hole": "DDH-03", "assay": 10.0},
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
            sample_value_col="assay",
            selected_metrics=["nearest_distance", "average_distance_knn"],
            average_distance_knn_k=2,
            closest_sample_id_cols=["hole"],
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_value_col="block_grade",
        )

        output_df = pd.read_csv(output_path)
        nearest_distance_column = result["nearest_distance_column"]

        assert result["selected_metrics"] == ["nearest_distance", "average_distance_knn"]
        assert result["average_distance_column"] is None
        assert result["average_distance_knn_column"] == "dom_Avg_Distance_KNN"
        assert result["average_distance_knn_k"] == 2
        assert result["closest_sample_id_column"] is None
        assert result["nearest_sample_value_column"] is None

        assert nearest_distance_column == "block_grade_dom_NN_Distance"
        assert nearest_distance_column in output_df.columns
        assert "dom_Avg_Distance_KNN" in output_df.columns
        assert "dom_Avg_Distance" not in output_df.columns
        assert "dom_Closest_Sample_ID" not in output_df.columns
        assert "Nearest_Sample_Value" not in output_df.columns

        assert output_df[nearest_distance_column].tolist() == pytest.approx([5.0, 5.0, 5.0])
        assert output_df["dom_Avg_Distance_KNN"].tolist() == pytest.approx([5.0, 10.0, 10.0])


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
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "sample_grade": 1.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "sample_grade": 2.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B", "sample_grade": 3.0},
                {"sx": 30.0, "sy": 5.0, "sz": 5.0, "sample_dom": "C", "sample_grade": 4.0},
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
            sample_value_col="sample_grade",
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
        )

        output_df = pd.read_csv(output_path)
        nearest_distance_column = result["nearest_distance_column"]

        assert result["matched_samples"] == 3
        assert result["unmatched_samples"] == 1
        assert output_df.loc[0, nearest_distance_column] == 3.0
        assert output_df.loc[1, nearest_distance_column] == 3.0
        assert output_df.loc[2, nearest_distance_column] == 0.0


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
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": 1.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": 2.0},
                {"sx": 40.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": 3.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B", "assay": 4.0},
                {"sx": 31.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B", "assay": 5.0},
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
            sample_value_col="assay",
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

        assert result["distance_count_columns"] == []
        assert result["distance_summary_thresholds"] == [5.0, 10.0]
        assert result["distance_count_step"] == 5.0
        assert result["distance_count_max_factor"] == 2
        assert result["summary_output_file"]

        for column_name in expected_columns:
            assert column_name not in output_df.columns

        summary_df = pd.read_csv(result["summary_output_file"])
        assert len(summary_df) == 2

        domain_a = summary_df.loc[summary_df["Domain"] == "A"].iloc[0]
        assert domain_a["Domain_Sample_Count"] == 3
        assert domain_a["Domain_Block_Volume"] == pytest.approx(2000.0)
        assert domain_a["Domain_Sample_Density"] == pytest.approx(0.0015)
        assert domain_a["Covered_Block_Count_LEQ_5"] == 2
        assert domain_a["Covered_Block_Volume_LEQ_5"] == pytest.approx(2000.0)
        assert domain_a["Coverage_Fraction_LEQ_5"] == pytest.approx(1.0)
        assert domain_a["Coverage_Density_LEQ_5"] == pytest.approx(0.0015)
        assert domain_a["Covered_Block_Count_LEQ_10"] == 2
        assert domain_a["Coverage_Density_LEQ_10"] == pytest.approx(0.0015)


def test_export_block_domain_sample_metrics_summary_works_for_inferred_domains():
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
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "assay": 1.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "assay": 2.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "assay": 3.0},
                {"sx": 28.0, "sy": 5.0, "sz": 5.0, "assay": 4.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            sample_value_col="assay",
            distance_count_step=5.0,
            distance_count_max_factor=2,
            block_x_col="x",
            block_y_col="y",
            block_z_col="z",
            block_domain_col="dom",
            block_size=(10, 10, 10),
        )

        summary_df = pd.read_csv(result["summary_output_file"])
        domain_b = summary_df.loc[summary_df["Domain"] == "B"].iloc[0]

        assert domain_b["Domain_Sample_Count"] == 2
        assert domain_b["Domain_Block_Volume"] == pytest.approx(1000.0)
        assert domain_b["Domain_Sample_Density"] == pytest.approx(0.002)
        assert domain_b["Covered_Block_Count_LEQ_5"] == 1
        assert domain_b["Covered_Block_Volume_LEQ_5"] == pytest.approx(1000.0)
        assert domain_b["Coverage_Density_LEQ_5"] == pytest.approx(0.002)
        assert domain_b["Covered_Block_Count_LEQ_10"] == 1
        assert domain_b["Coverage_Density_LEQ_10"] == pytest.approx(0.002)


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
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "assay": 1.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "assay": 2.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "assay": 3.0},
                {"sx": 28.0, "sy": 5.0, "sz": 5.0, "assay": 4.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            sample_value_col="assay",
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
        for column_name in expected_columns:
            assert column_name not in output_df.columns


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
        assert domain_a["Sample_To_Block_Distance_Ratio"] == pytest.approx(2.0)

        domain_b = output_df.loc[output_df["Domain"] == "B"].iloc[0]
        assert domain_b["Source_Sample_Count"] == 1
        assert domain_b["Domain_Block_Count"] == 1
        assert pd.isna(domain_b["Avg_Source_Sample_Distance"])
        assert domain_b["Avg_Block_To_Source_Sample_Distance"] == pytest.approx(0.0)
        assert pd.isna(domain_b["Sample_To_Block_Distance_Ratio"])
        assert "Avg_Domain_Block_Distance" not in output_df.columns


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
        assert "Avg_Domain_Block_Distance_X" not in output_df.columns
        assert "Avg_Domain_Block_Distance_Y" not in output_df.columns
        assert "Avg_Domain_Block_Distance_Z" not in output_df.columns

        domain_b = output_df.loc[output_df["Domain"] == "B"].iloc[0]
        assert pd.isna(domain_b["Avg_Source_Sample_Distance_X"])
        assert domain_b["Avg_Block_To_Source_Sample_Distance_X"] == pytest.approx(0.0)


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


def test_export_blocks_to_csv_streams_large_subblock_expansion(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        output_path = os.path.join(td, "expanded.csv")
        pd.DataFrame(
            [
                {"x": 2.5, "y": 5.0, "z": 5.0},
                {"x": 7.5, "y": 5.0, "z": 5.0},
                {"x": 12.5, "y": 5.0, "z": 5.0},
                {"x": 17.5, "y": 5.0, "z": 5.0},
            ]
        ).to_csv(blocks_path, index=False)

        progress_labels = []

        def _fake_iterate_csv_path_chunks_with_progress(path, progress_label, progress_callback=None,
                                                        header_line=1, **read_csv_kwargs):
            progress_labels.append(progress_label)
            yield from pd.read_csv(path, **read_csv_kwargs)

        monkeypatch.setattr(
            viewer_module,
            "iterate_csv_path_chunks_with_progress",
            _fake_iterate_csv_path_chunks_with_progress,
        )
        monkeypatch.setattr(viewer_module, "LARGE_BLOCK_FILE_THRESHOLD", 1)

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
        fake_blocks._ant_colony = type("FakeInterpolator", (), {"blocks": {}})()

        original_collect = viewer_module._collect_export_block_data
        monkeypatch.setattr(
            viewer_module,
            "_collect_export_block_data",
            lambda blocks: [
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
            ],
        )

        try:
            export_blocks_to_csv(fake_blocks, output_path)
        finally:
            monkeypatch.setattr(viewer_module, "_collect_export_block_data", original_collect)

        output_df = pd.read_csv(output_path)

        assert progress_labels == ["Reading source blocks for export expansion"]
        assert len(output_df) == 4
        assert output_df["x"].tolist() == [2.5, 7.5, 12.5, 17.5]
        assert output_df["Value"].tolist() == [10.0, 10.0, 20.0, 20.0]
        assert output_df["Domain"].tolist() == ["A", "A", "B", "B"]


def test_export_blocks_to_csv_streaming_leaves_unmatched_rows_blank(monkeypatch):
    with tempfile.TemporaryDirectory() as td:
        blocks_path = os.path.join(td, "blocks.csv")
        output_path = os.path.join(td, "expanded.csv")
        pd.DataFrame(
            [
                {"x": 2.5, "y": 5.0, "z": 5.0},
                {"x": 7.5, "y": 5.0, "z": 5.0},
                {"x": 45.0, "y": 5.0, "z": 5.0},
            ]
        ).to_csv(blocks_path, index=False)

        def _fake_iterate_csv_path_chunks_with_progress(path, progress_label, progress_callback=None,
                                                        header_line=1, **read_csv_kwargs):
            yield from pd.read_csv(path, **read_csv_kwargs)

        monkeypatch.setattr(
            viewer_module,
            "iterate_csv_path_chunks_with_progress",
            _fake_iterate_csv_path_chunks_with_progress,
        )
        monkeypatch.setattr(viewer_module, "LARGE_BLOCK_FILE_THRESHOLD", 1)

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
        fake_blocks._ant_colony = type("FakeInterpolator", (), {"blocks": {}})()

        original_collect = viewer_module._collect_export_block_data
        monkeypatch.setattr(
            viewer_module,
            "_collect_export_block_data",
            lambda blocks: [
                {
                    "_Grid_Index": (0, 0, 0),
                    "x": 0.0,
                    "y": 5.0,
                    "z": 5.0,
                    "Value": 10.0,
                    "Domain": "A",
                    "Source": "First Pass",
                },
            ],
        )

        try:
            export_blocks_to_csv(fake_blocks, output_path)
        finally:
            monkeypatch.setattr(viewer_module, "_collect_export_block_data", original_collect)

        output_df = pd.read_csv(output_path)

        assert len(output_df) == 3
        assert output_df["x"].tolist() == [2.5, 7.5, 45.0]
        assert output_df["Value"].tolist()[:2] == [10.0, 10.0]
        assert pd.isna(output_df.loc[2, "Value"])
        assert pd.isna(output_df.loc[2, "Domain"])


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


def test_detect_grid_rotation_ignores_small_jitter_on_axis_aligned_grid():
    rng = np.random.default_rng(2)
    xs = np.arange(0, 16) * 10.0
    ys = np.arange(0, 16) * 10.0
    zs = np.array([0.0, 10.0])
    points = np.array([(x, y, z) for z in zs for y in ys for x in xs], dtype=float)
    points[:, :2] += rng.normal(0.0, 0.6, size=points[:, :2].shape)

    rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(points, block_size_hint=(10, 10, 10))

    assert rotation_matrix.shape == (3, 3)
    assert rotation_center.shape == (3,)
    assert bool(is_rotated) is False


def test_detect_grid_rotation_detects_clear_rotation_with_noise():
    rng = np.random.default_rng(7)
    xs = np.arange(0, 16) * 10.0
    ys = np.arange(0, 16) * 10.0
    zs = np.array([0.0, 10.0])
    points = np.array([(x, y, z) for z in zs for y in ys for x in xs], dtype=float)

    angle_deg = 8.0
    theta = np.radians(angle_deg)
    rotation = np.array(
        [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ],
        dtype=float,
    )
    center_xy = points[:, :2].mean(axis=0)
    points[:, :2] = (points[:, :2] - center_xy) @ rotation.T + center_xy
    points[:, :2] += rng.normal(0.0, 0.35, size=points[:, :2].shape)

    rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(points, block_size_hint=(10, 10, 10))

    assert rotation_matrix.shape == (3, 3)
    assert rotation_center.shape == (3,)
    assert bool(is_rotated) is True


def test_detect_grid_rotation_large_unrotated_grid_does_not_trust_oversized_fallback_length():
    nx, ny, nz = 100, 100, 25
    xs = np.arange(nx, dtype=float) * 10.0
    ys = np.arange(ny, dtype=float) * 10.0
    zs = np.arange(nz, dtype=float) * 10.0
    points = np.array([(x, y, z) for z in zs for y in ys for x in xs], dtype=float)

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
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "sample_dom": "A", "assay": 1.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "sample_dom": "", "assay": 2.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "sample_dom": "B", "assay": 3.0},
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
            sample_value_col="assay",
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
                {"sx": 2.0, "sy": 5.0, "sz": 5.0, "assay": 1.0},
                {"sx": 18.0, "sy": 5.0, "sz": 5.0, "assay": 2.0},
                {"sx": 25.0, "sy": 5.0, "sz": 5.0, "assay": 3.0},
            ]
        ).to_csv(samples_path, index=False)

        result = export_block_domain_sample_metrics(
            samples_path,
            blocks_path,
            output_file=output_path,
            sample_x_col="sx",
            sample_y_col="sy",
            sample_z_col="sz",
            sample_value_col="assay",
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