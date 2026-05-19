"""Round-trip regression checks for BMF export backends."""
import os
import sys
import tempfile

import pandas as pd

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import bmf_standalone_exporter as bmf_exporter
from bmf_standalone_exporter import export_bmf, load_bmf_table


def run_tbms_config_text_round_trip():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "grid.csv")
        bmf_path = os.path.join(tmpdir, "grid.bmf")

        frame = pd.DataFrame(
            [
                {"x": 0.5, "y": 0.5, "z": 0.5, "grade": 1.25, "domain": "ore"},
                {"x": 1.5, "y": 0.5, "z": 0.5, "grade": 2.5, "domain": "waste"},
                {"x": 0.5, "y": 1.5, "z": 0.5, "grade": None, "domain": "ore"},
            ]
        )
        frame.to_csv(csv_path, index=False)

        export_result = export_bmf(csv_path, bmf_path, backend="tbms-config-text", dry_run=False)
        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]

        assert export_result["backend_summary"]["backend"] == "tbms-config-text"
        assert loaded["reader_mode"] == "tbms-config"
        assert loaded["report"].get("tbms_variant_family") == "tbms-config-text"
        assert loaded["row_count"] == 4
        assert loaded["rows_loaded"] == 4
        assert list(rows.columns) == ["grid_i", "grid_j", "grid_k", "x", "y", "z", "grade", "domain"]

        first = rows.iloc[0].to_dict()
        assert abs(first["grade"] - 1.25) < 1e-9
        assert first["domain"] == "ore"

        third = rows.iloc[2].to_dict()
        assert abs(third["grade"] - (-99.0)) < 1e-9
        assert third["domain"] == "ore"

        fourth = rows.iloc[3].to_dict()
        assert abs(fourth["grade"] - (-99.0)) < 1e-9
        assert fourth["domain"] == ""

    print("BMF export regression suite passed")


def run_dense_guard_backend_selection():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "wide_grid.csv")
        dense_bmf_path = os.path.join(tmpdir, "dense_grid.bmf")
        indexed_bmf_path = os.path.join(tmpdir, "indexed_grid.bmf")

        frame = pd.DataFrame(
            [
                {"x": 0.5, "y": 0.5, "z": 0.5, "grade": 1.0},
                {"x": 1000.5, "y": 0.5, "z": 0.5, "grade": 2.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        original_limit = bmf_exporter.MAX_DENSE_EXPORT_BYTES
        bmf_exporter.MAX_DENSE_EXPORT_BYTES = 1024
        try:
            try:
                export_bmf(
                    csv_path,
                    dense_bmf_path,
                    backend="tbms-config-text",
                    cell_size=[1.0, 1.0, 1.0],
                )
            except MemoryError as exc:
                assert "tbms-config-text" in str(exc)
                assert "tbms-experimental" in str(exc)
            else:
                raise AssertionError("Expected dense tbms-config-text export to fail the dense-size guard")

            result = export_bmf(
                csv_path,
                indexed_bmf_path,
                backend="tbms-experimental",
                cell_size=[1.0, 1.0, 1.0],
            )
            assert result["backend_summary"]["backend"] == "tbms-experimental"
            assert os.path.isfile(indexed_bmf_path)
        finally:
            bmf_exporter.MAX_DENSE_EXPORT_BYTES = original_limit


def run_coarse_cell_size_alignment_diagnostic():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "subblocks.csv")
        bmf_path = os.path.join(tmpdir, "subblocks.bmf")

        frame = pd.DataFrame(
            [
                {"x": 0.0, "y": 0.0, "z": 0.0, "grade": 1.0},
                {"x": 12.5, "y": 0.0, "z": 0.0, "grade": 2.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        try:
            export_bmf(
                csv_path,
                bmf_path,
                backend="tbms-config-text",
                cell_size=[25.0, 25.0, 15.0],
            )
        except ValueError as exc:
            message = str(exc)
            assert "smallest CSV coordinate/centroid increment" in message
            assert "likely minimum physical cell" not in message
            assert "sub-block rows" in message
            assert "aggregate" in message
        else:
            raise AssertionError("Expected coarse explicit cell size to fail with an alignment diagnostic")


def run_regularize_to_base_cell_export():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "subblock_centroids.csv")
        bmf_path = os.path.join(tmpdir, "regularized.bmf")

        frame = pd.DataFrame(
            [
                {"x": 1.0, "y": 1.0, "z": 1.0, "grade": 1.0, "domain": "ore"},
                {"x": 2.0, "y": 1.0, "z": 1.0, "grade": 3.0, "domain": "ore"},
                {"x": 11.0, "y": 1.0, "z": 1.0, "grade": 10.0, "domain": "waste"},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            cell_size=[10.0, 10.0, 10.0],
            regularize_to_base_block=True,
        )
        regularization = result["summary"]["regularization"]
        assert regularization["enabled"] is True
        assert regularization["input_rows"] == 3
        assert regularization["output_rows"] == 2
        assert regularization["aggregation"]["grade"] == "mean"
        assert regularization["aggregation"]["domain"] == "mode"

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 2
        assert abs(rows.iloc[0]["x"] - 5.0) < 1e-9
        assert abs(rows.iloc[0]["grade"] - 2.0) < 1e-9
        assert rows.iloc[0]["domain"] == "ore"
        assert abs(rows.iloc[1]["x"] - 15.0) < 1e-9
        assert abs(rows.iloc[1]["grade"] - 10.0) < 1e-9
        assert rows.iloc[1]["domain"] == "waste"


def run_regularize_ignores_value_exception_replacements():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "exception_subblocks.csv")
        bmf_path = os.path.join(tmpdir, "exception_regularized.bmf")

        frame = pd.DataFrame(
            [
                {"x": 1.0, "y": 1.0, "z": 1.0, "grade": "1.0"},
                {"x": 2.0, "y": 1.0, "z": 1.0, "grade": "3.0"},
                {"x": 3.0, "y": 1.0, "z": 1.0, "grade": "WOV"},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            cell_size=[10.0, 10.0, 10.0],
            value_cols=["grade"],
            column_types={"grade": "double"},
            value_exceptions={"grade": {"WOV": "-98"}},
            regularize_to_base_block=True,
        )

        assert result["summary"]["regularization"]["aggregation"]["grade"] == "mean"
        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 1
        assert abs(rows.iloc[0]["grade"] - 2.0) < 1e-9


def run_regularize_can_include_value_exception_replacements():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "included_exception_subblocks.csv")
        bmf_path = os.path.join(tmpdir, "included_exception_regularized.bmf")

        frame = pd.DataFrame(
            [
                {"x": 1.0, "y": 1.0, "z": 1.0, "grade": "1.0"},
                {"x": 2.0, "y": 1.0, "z": 1.0, "grade": "3.0"},
                {"x": 3.0, "y": 1.0, "z": 1.0, "grade": "WOV"},
            ]
        )
        frame.to_csv(csv_path, index=False)

        export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            cell_size=[10.0, 10.0, 10.0],
            value_cols=["grade"],
            column_types={"grade": "double"},
            value_exceptions={
                "grade": {
                    "WOV": {
                        "replacement": "-98",
                        "include_in_regularization": True,
                    }
                }
            },
            regularize_to_base_block=True,
        )

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 1
        assert abs(rows.iloc[0]["grade"] - ((1.0 + 3.0 - 98.0) / 3.0)) < 1e-9


if __name__ == "__main__":
    run_tbms_config_text_round_trip()
    run_dense_guard_backend_selection()
    run_coarse_cell_size_alignment_diagnostic()
    run_regularize_to_base_cell_export()
    run_regularize_ignores_value_exception_replacements()
    run_regularize_can_include_value_exception_replacements()
