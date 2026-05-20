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
        assert loaded["row_count"] == 3
        assert loaded["rows_loaded"] == 3
        for column_name in ["x", "y", "z", "__lower_x", "__upper_x", "grade", "domain"]:
            assert column_name in rows.columns

        first = rows.iloc[0].to_dict()
        assert abs(first["grade"] - 1.25) < 1e-9
        assert first["domain"] == "ore"
        assert abs(first["__lower_x"] - 0.0) < 1e-9
        assert abs(first["__upper_x"] - 1.0) < 1e-9

        third = rows.iloc[2].to_dict()
        assert abs(third["grade"] - (-99.0)) < 1e-9
        assert third["domain"] == "ore"

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
                    regularize_to_base_block=True,
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


def run_grid_index_columns_drive_csv_export_grid():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "indexed_source.csv")
        bmf_path = os.path.join(tmpdir, "indexed_source.bmf")

        frame = pd.DataFrame(
            [
                {"grid_i": 0, "grid_j": 0, "grid_k": 0, "x": 100.0, "y": 10.0, "z": 5.0, "grade": 1.0},
                {"grid_i": 1, "grid_j": 0, "grid_k": 0, "x": 500.0, "y": 10.0, "z": 5.0, "grade": 2.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            value_cols=["grade"],
        )
        assert result["summary"]["grid"]["index_source"] == "grid_i/grid_j/grid_k"
        assert result["summary"]["grid"]["dimensions"] == [2, 1, 1]

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 2
        assert list(rows["grade"]) == [1.0, 2.0]


def run_dense_metadata_prefers_parent_block_size():
    metadata = {
        "parent_block_size": [25.0, 25.0, 15.0],
        "subblock_factors": [4, 4, 4],
        "minimum_corner": [11535.0, 99020.0, 700.0],
    }

    dense_geometry = bmf_exporter.infer_bmf_export_geometry_from_metadata(
        metadata,
        regularize_to_base_block=False,
        dense_regular_grid=True,
    )
    assert dense_geometry["cell_size"] == [25.0, 25.0, 15.0]
    assert dense_geometry["cell_size_source"] == "parent_block_size"
    assert "dense_metadata_note" in dense_geometry

    indexed_geometry = bmf_exporter.infer_bmf_export_geometry_from_metadata(
        metadata,
        regularize_to_base_block=False,
        dense_regular_grid=False,
    )
    assert indexed_geometry["cell_size"] == [6.25, 6.25, 3.75]
    assert indexed_geometry["cell_size_source"] == "parent_block_size/subblock_factors"


def run_tbms_config_text_irregular_subblocks():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "mixed_subblocks.csv")
        bmf_path = os.path.join(tmpdir, "mixed_subblocks.bmf")

        with open(csv_path, "w", encoding="utf-8", newline="") as handle:
            handle.write("# Parent block size: 10, 10, 10\n")
            handle.write("# Size in parent blocks: 2, 1, 1\n")
            handle.write("# Minimum corner: 0, 0, 0\n")
            handle.write("# Sub-blocks: octree 2, 2, 2\n")
            handle.write("x,y,z,grade,domain\n")
            handle.write("5,5,5,1,parent\n")
            handle.write("12.5,2.5,2.5,2,sub_a\n")
            handle.write("17.5,2.5,2.5,3,sub_b\n")

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            header_line=1,
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["backend_summary"]["row_count"] == 3

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["reader_mode"] == "tbms-config"
        assert loaded["row_count"] == 3
        assert "__lower_x" in rows.columns
        assert "__upper_x" in rows.columns
        assert list(rows["domain"]) == ["parent", "sub_a", "sub_b"]
        assert abs(rows.iloc[0]["__lower_x"] - 0.0) < 1e-6
        assert abs(rows.iloc[0]["__upper_x"] - 10.0) < 1e-6
        assert abs(rows.iloc[1]["__lower_x"] - 10.0) < 1e-6
        assert abs(rows.iloc[1]["__upper_x"] - 15.0) < 1e-6
        assert abs(rows.iloc[2]["__lower_x"] - 15.0) < 1e-6
        assert abs(rows.iloc[2]["__upper_x"] - 20.0) < 1e-6


def run_tbms_config_text_irregular_explicit_extents():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "explicit_extents.csv")
        bmf_path = os.path.join(tmpdir, "explicit_extents.bmf")

        with open(csv_path, "w", encoding="utf-8", newline="") as handle:
            handle.write("x,y,z,__lower_x,__upper_x,__lower_y,__upper_y,__lower_z,__upper_z,grade\n")
            handle.write("5,5,5,0,10,0,10,0,10,1\n")
            handle.write("12.5,2.5,2.5,10,15,0,5,0,5,2\n")

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            header_line=1,
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["index_source"] == "irregular_lower_upper_columns"
        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert list(rows["grade"]) == [1.0, 2.0]
        assert abs(rows.iloc[1]["__upper_x"] - 15.0) < 1e-6


def run_tbms_config_text_irregular_dimension_columns():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "dimension_columns.csv")
        bmf_path = os.path.join(tmpdir, "dimension_columns.bmf")

        frame = pd.DataFrame(
            [
                {"x": 5.0, "y": 5.0, "z": 5.0, "dX": 10.0, "dY": 10.0, "dZ": 10.0, "grade": 1.0},
                {"x": 12.5, "y": 2.5, "z": 2.5, "dX": 5.0, "dY": 5.0, "dZ": 5.0, "grade": 2.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["index_source"] == "irregular_dimension_columns"
        assert result["summary"]["grid"]["irregular_width_source"] == "dimension_columns"
        assert result["summary"]["value_columns"] == ["grade"]
        determination = result["summary"]["block_size_determination"]
        assert determination["method"] == "explicit_dimension_columns"
        assert determination["source"] == "csv_dimension_columns"
        assert determination["geometry_column_mapping"]["size_cols"] == {"dx": "dX", "dy": "dY", "dz": "dZ"}

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert list(rows["grade"]) == [1.0, 2.0]
        assert "dX" not in rows.columns
        assert list(rows["__lower_x"].round(6)) == [0.0, 10.0]
        assert list(rows["__upper_x"].round(6)) == [10.0, 15.0]


def run_tbms_config_text_irregular_custom_extent_columns():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "custom_extents.csv")
        bmf_path = os.path.join(tmpdir, "custom_extents.bmf")

        frame = pd.DataFrame(
            [
                {"cx": 5.0, "cy": 5.0, "cz": 5.0, "lx": 0.0, "ux": 10.0, "ly": 0.0, "uy": 10.0, "lz": 0.0, "uz": 10.0, "grade": 1.0},
                {"cx": 12.5, "cy": 2.5, "cz": 2.5, "lx": 10.0, "ux": 15.0, "ly": 0.0, "uy": 5.0, "lz": 0.0, "uz": 5.0, "grade": 2.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            x_col="cx",
            y_col="cy",
            z_col="cz",
            extent_cols=["lx", "ux", "ly", "uy", "lz", "uz"],
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["index_source"] == "irregular_mapped_lower_upper_columns"
        assert result["summary"]["grid"]["irregular_width_source"] == "mapped_lower_upper_columns"
        determination = result["summary"]["block_size_determination"]
        assert determination["method"] == "explicit_row_extents"
        assert determination["source"] == "mapped_csv_lower_upper_extent_columns"
        assert determination["geometry_column_mapping"]["extent_cols"]["__lower_x"] == "lx"

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert "lx" not in rows.columns
        assert list(rows["__lower_x"].round(6)) == [0.0, 10.0]
        assert list(rows["__upper_x"].round(6)) == [10.0, 15.0]


def run_tbms_config_text_irregular_centroid_hierarchy_fallback():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "centroid_only_subblocks.csv")
        bmf_path = os.path.join(tmpdir, "centroid_only_subblocks.bmf")

        frame = pd.DataFrame(
            [
                {"x": 0.5, "y": 4.0, "z": 4.0, "grade": 1.0},
                {"x": 1.5, "y": 4.0, "z": 4.0, "grade": 2.0},
                {"x": 4.5, "y": 4.0, "z": 4.0, "grade": 3.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            cell_size=[8.0, 8.0, 8.0],
            origin=[0.0, 0.0, 0.0],
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["index_source"] == "irregular_xyz_hierarchy_inference"
        assert result["summary"]["grid"]["irregular_width_source"] == "centroid_hierarchy_inference"
        hierarchy = result["summary"]["grid"]["irregular_hierarchy"]
        assert hierarchy["axes"]["x"]["quantum"] == 1.0
        assert hierarchy["axes"]["x"]["base"] == 8.0
        assert [row["size"] for row in hierarchy["axes"]["x"]["levels"]] == [1.0, 2.0, 4.0, 8.0]
        assert "Axis X: quantum=1, base=8" in result["summary"]["grid"]["irregular_hierarchy_text"]

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        x_widths = rows["__upper_x"] - rows["__lower_x"]
        assert list(x_widths.round(6)) == [1.0, 1.0, 1.0]
        assert 3.0 not in set(float(value) for value in x_widths.round(6))
        assert list(rows["__lower_x"].round(6)) == [0.0, 1.0, 4.0]
        assert list(rows["__upper_x"].round(6)) == [1.0, 2.0, 5.0]


def run_tbms_config_text_irregular_statistical_base_size_fallback():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "centroid_only_no_metadata.csv")
        bmf_path = os.path.join(tmpdir, "centroid_only_no_metadata.bmf")

        frame = pd.DataFrame(
            [
                {"x": 0.5, "y": 0.5, "z": 0.5, "grade": 1.0},
                {"x": 1.5, "y": 1.5, "z": 1.5, "grade": 2.0},
                {"x": 4.5, "y": 4.5, "z": 4.5, "grade": 3.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        statistical_base = result["summary"]["csv_statistical_base_block"]
        assert statistical_base["status"] == "ok"
        assert statistical_base["cell_size"] == [1.0, 1.0, 1.0]
        assert statistical_base["origin"] == [0.0, 0.0, 0.0]
        assert statistical_base["axes"]["x"]["parent_base_size_ambiguous"] is True
        assert statistical_base["axes"]["x"]["enclosing_power_of_two_size"] == 8.0
        assert statistical_base["axes"]["x"]["hierarchy_spacing_histogram"] == [{"spacing": 1.0, "count": 1}]
        assert statistical_base["axes"]["x"]["discarded_spacing_histogram"] == [{"spacing": 3.0, "count": 1}]
        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["cell_size"] == [1.0, 1.0, 1.0]
        assert result["summary"]["grid"]["index_source"] == "irregular_xyz_hierarchy_inference"
        determination = result["summary"]["block_size_determination"]
        assert determination["method"] == "statistical_centroid_spacing"
        assert determination["source"] == "centroid_spacing_quantum_fallback"
        assert "ambiguous" in determination["message"]
        assert determination["cell_size"] == [1.0, 1.0, 1.0]

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 3
        assert list((rows["__upper_x"] - rows["__lower_x"]).round(6)) == [1.0, 1.0, 1.0]


def run_tbms_config_text_statistical_frequency_break_infers_base_size():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "centroid_frequency_break.csv")
        bmf_path = os.path.join(tmpdir, "centroid_frequency_break.bmf")

        diffs = [1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 8]
        centers = [0.5]
        for diff in diffs:
            centers.append(centers[-1] + diff)
        frame = pd.DataFrame(
            {"x": centers, "y": centers, "z": centers, "grade": list(range(len(centers)))}
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        statistical_base = result["summary"]["csv_statistical_base_block"]
        assert statistical_base["status"] == "ok"
        assert statistical_base["cell_size"] == [4.0, 4.0, 4.0]
        assert result["summary"]["grid"]["cell_size"] == [4.0, 4.0, 4.0]
        x_inference = statistical_base["axes"]["x"]["base_size_inference"]
        assert x_inference["status"] == "inferred"
        assert x_inference["candidate_size"] == 4.0
        assert x_inference["drop_to_size"] == 8.0
        assert x_inference["drop_from_count"] == 5
        assert x_inference["drop_to_count"] == 1
        assert statistical_base["axes"]["x"]["parent_base_size_ambiguous"] is False
        determination = result["summary"]["block_size_determination"]
        assert determination["method"] == "statistical_centroid_spacing"
        assert determination["source"] == "centroid_spacing_frequency_break"
        assert "frequency break" in determination["message"]
        assert determination["cell_size"] == [4.0, 4.0, 4.0]
        assert determination["statistical_inference"]["inferred_axes"] == ["X", "Y", "Z"]

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == len(centers)
        assert list((rows["__upper_x"] - rows["__lower_x"]).round(6)) == [1.0] * len(centers)


def run_tbms_config_text_regular_centroids_still_preserve_source_rows():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "regular_centroids_no_metadata.csv")
        bmf_path = os.path.join(tmpdir, "regular_centroids_no_metadata.bmf")

        frame = pd.DataFrame(
            [
                {"x": 0.5, "y": 0.5, "z": 0.5, "grade": 1.0},
                {"x": 1.5, "y": 0.5, "z": 0.5, "grade": 2.0},
                {"x": 2.5, "y": 0.5, "z": 0.5, "grade": 3.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["irregular_width_source"] == "centroid_hierarchy_inference"
        assert result["summary"]["grid"]["cell_size"] == [1.0, 1.0, 1.0]

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 3
        assert "__lower_x" in rows.columns
        assert "__upper_x" in rows.columns
        assert list(rows["__lower_x"].round(6)) == [0.0, 1.0, 2.0]
        assert list(rows["__upper_x"].round(6)) == [1.0, 2.0, 3.0]


def run_tbms_config_text_metadata_parent_still_infers_row_extents():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "metadata_parent_no_extents.csv")
        bmf_path = os.path.join(tmpdir, "metadata_parent_no_extents.bmf")

        with open(csv_path, "w", encoding="utf-8", newline="") as handle:
            handle.write("# Parent block size: 8 8 8\n")
            handle.write("# Sub-blocks: Regular 8 8 8\n")
            handle.write("# Minimum corner: 0 0 0\n")
            handle.write("x,y,z,grade\n")
            handle.write("0.5,4,4,1\n")
            handle.write("1.5,4,4,2\n")
            handle.write("4.5,4,4,3\n")

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            header_line=1,
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["cell_size"] == [8.0, 8.0, 8.0]
        assert result["summary"]["grid"]["index_source"] == "irregular_xyz_metadata_parent"
        assert result["summary"]["grid"]["irregular_width_source"] == "centroid_hierarchy_inference"
        determination = result["summary"]["block_size_determination"]
        assert determination["method"] == "metadata_parent_size"
        assert determination["source"] == "parent_block_size"
        assert determination["cell_size"] == [8.0, 8.0, 8.0]

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 3
        assert list((rows["__upper_x"] - rows["__lower_x"]).round(6)) == [1.0, 1.0, 1.0]
        assert list(rows["__lower_x"].round(6)) == [0.0, 1.0, 4.0]
        assert list(rows["__upper_x"].round(6)) == [1.0, 2.0, 5.0]


def run_tbms_config_text_manual_parent_still_infers_row_extents():
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "manual_parent_no_extents.csv")
        bmf_path = os.path.join(tmpdir, "manual_parent_no_extents.bmf")

        frame = pd.DataFrame(
            [
                {"x": 0.5, "y": 4.0, "z": 4.0, "grade": 1.0},
                {"x": 1.5, "y": 4.0, "z": 4.0, "grade": 2.0},
                {"x": 2.5, "y": 4.0, "z": 4.0, "grade": 3.0},
            ]
        )
        frame.to_csv(csv_path, index=False)

        result = export_bmf(
            csv_path,
            bmf_path,
            backend="tbms-config-text",
            cell_size=[8.0, 8.0, 8.0],
            origin=[0.0, 0.0, 0.0],
            value_cols=["grade"],
            regularize_to_base_block=False,
        )

        assert result["summary"]["grid"]["is_irregular"] is True
        assert result["summary"]["grid"]["cell_size"] == [8.0, 8.0, 8.0]
        assert result["summary"]["grid"]["index_source"] == "irregular_xyz_hierarchy_inference"
        assert result["summary"]["grid"]["irregular_width_source"] == "centroid_hierarchy_inference"
        determination = result["summary"]["block_size_determination"]
        assert determination["method"] == "manual_parent_size"
        assert determination["source"] == "user_cell_size"
        assert determination["cell_size"] == [8.0, 8.0, 8.0]

        loaded = load_bmf_table(bmf_path)
        rows = loaded["dataframe"]
        assert loaded["row_count"] == 3
        assert list((rows["__upper_x"] - rows["__lower_x"]).round(6)) == [1.0, 1.0, 1.0]
        assert list(rows["__lower_x"].round(6)) == [0.0, 1.0, 2.0]
        assert list(rows["__upper_x"].round(6)) == [1.0, 2.0, 3.0]


if __name__ == "__main__":
    run_tbms_config_text_round_trip()
    run_dense_guard_backend_selection()
    run_coarse_cell_size_alignment_diagnostic()
    run_regularize_to_base_cell_export()
    run_regularize_ignores_value_exception_replacements()
    run_regularize_can_include_value_exception_replacements()
    run_grid_index_columns_drive_csv_export_grid()
    run_dense_metadata_prefers_parent_block_size()
    run_tbms_config_text_irregular_subblocks()
    run_tbms_config_text_irregular_explicit_extents()
    run_tbms_config_text_irregular_dimension_columns()
    run_tbms_config_text_irregular_custom_extent_columns()
    run_tbms_config_text_irregular_centroid_hierarchy_fallback()
    run_tbms_config_text_irregular_statistical_base_size_fallback()
    run_tbms_config_text_statistical_frequency_break_infers_base_size()
    run_tbms_config_text_regular_centroids_still_preserve_source_rows()
    run_tbms_config_text_metadata_parent_still_infers_row_extents()
    run_tbms_config_text_manual_parent_still_infers_row_extents()
