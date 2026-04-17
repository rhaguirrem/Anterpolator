import os
import sys
import tempfile

import pandas as pd

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from anterpolator3DViewer import export_domained_samples_from_blocks, load_large_blocks_metadata


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