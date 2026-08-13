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


def test_string_theory_can_skip_same_hole_connections():
    sample_blocks = {
        (0, 0, 0): 1.0,
        (2, 0, 0): 2.0,
        (4, 0, 0): 3.0,
    }
    interpolator = StringTheoryInterpolator(
        distance_threshold=5.0,
        grade_difference=10.0,
        avoid_same_hole_id=True,
    )
    interpolator.initialize_blocks(
        sample_blocks,
        (5, 1, 1),
        (0, 0, 0),
        (1, 1, 1),
        sample_hole_id_mapping={
            (0, 0, 0): 'DH_01',
            (2, 0, 0): 'DH_01',
            (4, 0, 0): 'DH_02',
        },
    )

    interpolator.run_iteration((5, 1, 1))

    assert interpolator.stats['same_hole_rejected'] >= 1
    assert ((0, 0, 0), (2, 0, 0), [(0, 0, 0), (1, 0, 0), (2, 0, 0)]) not in interpolator.paths
    assert interpolator.stats['connections_made'] == 1


def test_generate_statistics_exports_per_domain_breakdown_for_domain_mode():
    sample_blocks = {
        (0, 0, 0): 1.0,
        (4, 0, 0): 2.0,
        (0, 2, 0): 3.0,
        (4, 2, 0): 4.0,
    }
    interpolator = StringTheoryInterpolator(
        distance_threshold=5.0,
        interpolation_target='domain',
    )
    interpolator.initialize_blocks(
        sample_blocks,
        (5, 3, 1),
        (0, 0, 0),
        (1, 1, 1),
        sample_domain_mapping={
            (0, 0, 0): 'Alpha',
            (4, 0, 0): 'Alpha',
            (0, 2, 0): 'Beta',
            (4, 2, 0): 'Beta',
        },
    )

    interpolator.run_iteration((5, 3, 1))

    with tempfile.TemporaryDirectory() as td:
        interpolator.generate_statistics(td, 'Global')

        stats_dir = os.path.join(td, 'StringTheory_Stats')
        global_csv_path = os.path.join(stats_dir, 'path_stats_Global.csv')
        alpha_csv_path = os.path.join(stats_dir, 'path_stats_Alpha.csv')
        beta_csv_path = os.path.join(stats_dir, 'path_stats_Beta.csv')
        alpha_sample_csv_path = os.path.join(stats_dir, 'sample_connection_stats_Alpha.csv')
        beta_sample_csv_path = os.path.join(stats_dir, 'sample_connection_stats_Beta.csv')

        assert os.path.exists(global_csv_path)
        assert os.path.exists(alpha_csv_path)
        assert os.path.exists(beta_csv_path)
        assert os.path.exists(alpha_sample_csv_path)
        assert os.path.exists(beta_sample_csv_path)

        with open(global_csv_path, newline='', encoding='utf-8') as handle:
            global_rows = list(csv.DictReader(handle))
        assert {row['Category'] for row in global_rows} == {'Alpha', 'Beta'}

        with open(alpha_csv_path, newline='', encoding='utf-8') as handle:
            alpha_rows = list(csv.DictReader(handle))
        with open(beta_csv_path, newline='', encoding='utf-8') as handle:
            beta_rows = list(csv.DictReader(handle))

        assert alpha_rows
        assert beta_rows
        assert {row['Category'] for row in alpha_rows} == {'Alpha'}
        assert {row['Category'] for row in beta_rows} == {'Beta'}

        with open(alpha_sample_csv_path, newline='', encoding='utf-8') as handle:
            alpha_sample_rows = list(csv.DictReader(handle))
        with open(beta_sample_csv_path, newline='', encoding='utf-8') as handle:
            beta_sample_rows = list(csv.DictReader(handle))

        assert {(int(row['sample_x']), int(row['sample_y']), int(row['sample_z'])) for row in alpha_sample_rows} == {
            (0, 0, 0),
            (4, 0, 0),
        }
        assert {(int(row['sample_x']), int(row['sample_y']), int(row['sample_z'])) for row in beta_sample_rows} == {
            (0, 2, 0),
            (4, 2, 0),
        }