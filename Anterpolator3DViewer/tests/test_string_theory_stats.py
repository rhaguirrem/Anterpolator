import csv
import os
import sys
import tempfile

os.environ.setdefault("MPLBACKEND", "Agg")

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from string_theory_interpolator import StringTheoryInterpolator


def test_generate_statistics_exports_sample_connection_histogram_and_csv():
    sample_blocks = {
        (0, 0, 0): 1.0,
        (2, 0, 0): 2.0,
        (4, 0, 0): 3.0,
    }
    interpolator = StringTheoryInterpolator()
    interpolator.initialize_blocks(sample_blocks, (5, 1, 1), (0, 0, 0), (1, 1, 1))
    interpolator.paths = [
        ((0, 0, 0), (2, 0, 0), [(0, 0, 0), (1, 0, 0), (2, 0, 0)]),
        ((0, 0, 0), (4, 0, 0), [(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (4, 0, 0)]),
        ((2, 0, 0), (4, 0, 0), [(2, 0, 0), (3, 0, 0), (4, 0, 0)]),
    ]

    with tempfile.TemporaryDirectory() as td:
        interpolator.generate_statistics(td, "Demo")

        stats_dir = os.path.join(td, "StringTheory_Stats")
        histogram_path = os.path.join(stats_dir, "connection_count_histogram_Demo.png")
        csv_path = os.path.join(stats_dir, "sample_connection_stats_Demo.csv")

        assert os.path.exists(histogram_path)
        assert os.path.exists(csv_path)

        with open(csv_path, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))

        counts_by_value = {
            float(row["sample_value"]): int(row["connection_count"])
            for row in rows
        }
        assert counts_by_value == {
            1.0: 2,
            2.0: 1,
            3.0: 0,
        }