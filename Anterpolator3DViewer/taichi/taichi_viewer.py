import argparse
import os
import sys


CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from runtime import create_viewer_from_files


def parse_args():
    parser = argparse.ArgumentParser(description='Anterpolator Taichi viewer')
    parser.add_argument('--samples', required=True, help='Samples CSV file')
    parser.add_argument('--blocks', help='Blocks CSV file')
    parser.add_argument('--interpolation', help='Interpolation CSV file')
    parser.add_argument('--color', help='LFC color file')
    parser.add_argument('--samples-delimiter', help='Samples delimiter override')
    parser.add_argument('--blocks-delimiter', help='Blocks/interpolation delimiter override')
    parser.add_argument('--samples-header-line', type=int, default=1, help='1-based samples header line')
    parser.add_argument('--blocks-header-line', type=int, default=1, help='1-based blocks header line')
    parser.add_argument('--sample-x-col', help='Samples X column name')
    parser.add_argument('--sample-y-col', help='Samples Y column name')
    parser.add_argument('--sample-z-col', help='Samples Z column name')
    parser.add_argument('--sample-value-col', help='Samples value column name')
    parser.add_argument('--block-size', nargs=3, type=float, default=[10.0, 10.0, 10.0], metavar=('SX', 'SY', 'SZ'))
    parser.add_argument('--value-filter', type=float, default=0.0)
    return parser.parse_args()


def main():
    args = parse_args()
    viewer = create_viewer_from_files(
        samples_file=args.samples,
        color_file=args.color,
        interpolation_file=args.interpolation,
        blocks_file=args.blocks,
        samples_delimiter=args.samples_delimiter,
        blocks_delimiter=args.blocks_delimiter,
        samples_header_line=args.samples_header_line,
        blocks_header_line=args.blocks_header_line,
        sample_x_col=args.sample_x_col,
        sample_y_col=args.sample_y_col,
        sample_z_col=args.sample_z_col,
        sample_value_col=args.sample_value_col,
        block_size=args.block_size,
        value_filter=args.value_filter,
    )
    viewer.run()


if __name__ == '__main__':
    main()
