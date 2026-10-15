import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from anterpolator3DViewer import export_matching_value_sample_distances
import anterpolator3DViewer as viewer_module


def _export(tmp_path, blocks, samples, mappings=None, **kwargs):
    blocks_path = tmp_path / 'blocks.csv'
    samples_path = tmp_path / 'samples.csv'
    pd.DataFrame(blocks).to_csv(blocks_path, index=False)
    pd.DataFrame(samples).to_csv(samples_path, index=False)
    result = export_matching_value_sample_distances(
        str(samples_path), str(blocks_path),
        mappings or [{'block_col': 'RQDpc', 'sample_col': 'RQD'}],
        sample_domain_col='sample_UG', block_domain_col='UG', **kwargs,
    )
    return pd.read_csv(result['output_file']), result


def test_multiple_fields_match_values_not_just_nearest_samples(tmp_path):
    blocks = [
        {'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 20, 'lith': ' granite ', 'id': 2},
        {'x': 0, 'y': 0, 'z': 0, 'UG': 'B', 'RQDpc': 20, 'lith': 'basalt', 'id': 1},
    ]
    samples = [
        {'x': 0, 'y': 0, 'z': 0, 'sample_UG': 'A', 'RQD': 10, 'rock': 'basalt'},
        {'x': 3, 'y': 4, 'z': 12, 'sample_UG': 'A', 'RQD': 20, 'rock': 'granite'},
        {'x': 1, 'y': 0, 'z': 0, 'sample_UG': 'B', 'RQD': 20, 'rock': 'granite'},
        {'x': 2, 'y': 0, 'z': 0, 'sample_UG': 'C', 'RQD': None, 'rock': 'basalt'},
    ]
    events = []
    output, result = _export(
        tmp_path, blocks, samples,
        [{'block_col': 'RQDpc', 'sample_col': 'RQD'}, {'block_col': 'lith', 'sample_col': 'rock'}],
        progress_callback=lambda value, maximum, message: events.append((value, maximum)),
    )
    assert result['output_columns'] == ['RQDpc_UG_MinDist', 'RQDpc_AbsMinDist', 'lith_UG_MinDist', 'lith_AbsMinDist']
    assert output['RQDpc_UG_MinDist'].tolist() == [13, 1]
    assert output['RQDpc_AbsMinDist'].tolist() == [1, 1]
    assert output.loc[0, 'lith_UG_MinDist'] == 13
    assert output.loc[0, 'lith_AbsMinDist'] == 1
    assert pd.isna(output.loc[1, 'lith_UG_MinDist'])
    assert output.loc[1, 'lith_AbsMinDist'] == 0
    pd.testing.assert_frame_equal(output[list(blocks[0])], pd.DataFrame(blocks))
    assert events[0][0] == 0 and events[-1] == (100, 100)
    assert all(first[0] <= second[0] for first, second in zip(events, events[1:]))


def test_missing_matches_domains_and_invalid_coordinates_are_reported(tmp_path):
    blocks = [
        {'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 20},
        {'x': 0, 'y': 0, 'z': 0, 'UG': None, 'RQDpc': 20},
        {'x': np.inf, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 20},
        {'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': None},
        {'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 99},
    ]
    samples = [
        {'x': np.inf, 'y': 0, 'z': 0, 'sample_UG': 'A', 'RQD': 20},
        {'x': 'bad', 'y': 0, 'z': 0, 'sample_UG': 'A', 'RQD': 20},
        {'x': 1, 'y': 0, 'z': 0, 'sample_UG': None, 'RQD': 20},
        {'x': 0, 'y': 0, 'z': 0, 'sample_UG': 'A', 'RQD': None},
        {'x': 2, 'y': 0, 'z': 0, 'sample_UG': 'not in blocks', 'RQD': 20},
    ]
    output, result = _export(tmp_path, blocks, samples)
    assert output['RQDpc_UG_MinDist'].isna().all()
    assert output['RQDpc_AbsMinDist'].iloc[:2].tolist() == [1, 1]
    assert output['RQDpc_AbsMinDist'].iloc[2:].isna().all()
    assert result['invalid_coordinate_samples'] == 2
    assert result['invalid_coordinate_blocks'] == 1
    assert result['field_results'][0]['same_domain_unmatched_blocks'] == 5
    assert result['field_results'][0]['absolute_unmatched_blocks'] == 3


def test_numeric_tokens_and_domains_match_across_csv_types(tmp_path):
    output, _ = _export(
        tmp_path,
        [{'x': 0, 'y': 0, 'z': 0, 'UG': 1, 'RQDpc': 20}],
        [
            {'x': 3, 'y': 0, 'z': 0, 'sample_UG': 1.0, 'RQD': '20.0'},
            {'x': 1, 'y': 0, 'z': 0, 'sample_UG': 'B', 'RQD': 'different'},
        ],
    )
    assert output.loc[0, 'RQDpc_UG_MinDist'] == 3
    assert output.loc[0, 'RQDpc_AbsMinDist'] == 3


def test_matches_brute_force_for_multiple_groups(tmp_path):
    rng = np.random.default_rng(42)
    blocks = pd.DataFrame(rng.normal(size=(80, 3)), columns=['x', 'y', 'z'])
    samples = pd.DataFrame(rng.normal(size=(50, 3)), columns=['x', 'y', 'z'])
    blocks['UG'] = rng.choice(['A', 'B', 'C'], len(blocks))
    samples['sample_UG'] = rng.choice(['A', 'B', 'C', 'D'], len(samples))
    blocks['RQDpc'] = rng.integers(0, 5, len(blocks))
    samples['RQD'] = rng.integers(0, 4, len(samples))
    output, _ = _export(tmp_path, blocks, samples)
    for row_index, block in blocks.iterrows():
        for same_domain in (True, False):
            mask = samples['RQD'].eq(block['RQDpc'])
            if same_domain:
                mask &= samples['sample_UG'].eq(block['UG'])
            coords = samples.loc[mask, ['x', 'y', 'z']].to_numpy()
            column = 'RQDpc_UG_MinDist' if same_domain else 'RQDpc_AbsMinDist'
            if len(coords):
                expected = np.linalg.norm(coords - block[['x', 'y', 'z']].to_numpy(dtype=float), axis=1).min()
                assert output.loc[row_index, column] == pytest.approx(expected)
            else:
                assert pd.isna(output.loc[row_index, column])


def test_filters_headers_delimiters_and_empty_sets(tmp_path):
    blocks_path, samples_path = tmp_path / 'blocks.csv', tmp_path / 'samples.csv'
    blocks_path.write_text('metadata\nx;y;z;UG;RQDpc;keep\n0;0;0;A;20;yes\n1;0;0;B;20;no\n')
    samples_path.write_text('metadata\nsx|sy|sz|sample_UG|RQD|keep\n1|0|0|A|20|no\n5|0|0|A|20|yes\n')
    kwargs = dict(
        samples_file=str(samples_path), blocks_file=str(blocks_path),
        field_mappings=[{'block_col': 'RQDpc', 'sample_col': 'RQD'}],
        sample_domain_col='sample_UG', block_domain_col='UG',
        samples_delimiter='|', blocks_delimiter=';', samples_header_line=2, blocks_header_line=2,
        sample_x_col='sx', sample_y_col='sy', sample_z_col='sz',
        sample_filters=[{'field': 'keep', 'type': 'categorical', 'values': ['yes']}],
        block_filters=[{'field': 'keep', 'type': 'categorical', 'values': ['yes']}],
    )
    result = export_matching_value_sample_distances(**kwargs)
    output = pd.read_csv(result['output_file'])
    assert len(output) == 1 and output.loc[0, 'RQDpc_UG_MinDist'] == 5
    assert result['filtered_samples'] == 1 and result['filters_applied']
    kwargs['sample_filters'][0]['values'] = ['absent']
    result = export_matching_value_sample_distances(**kwargs)
    assert pd.read_csv(result['output_file'])['RQDpc_AbsMinDist'].isna().all()
    kwargs['block_filters'][0]['values'] = ['absent']
    result = export_matching_value_sample_distances(**kwargs)
    assert pd.read_csv(result['output_file']).empty
    assert result['processed_blocks'] == 0


def test_reading_progress_and_projected_sample_columns(tmp_path, monkeypatch):
    reads = []
    events = []
    original_read = viewer_module.read_selected_columns_with_header

    def capture_read(path, delimiter, header_line, selected_columns, **kwargs):
        reads.append((os.path.basename(path), selected_columns, kwargs.get('progress_label')))
        return original_read(path, delimiter, header_line, selected_columns, **kwargs)

    monkeypatch.setattr(viewer_module, 'read_selected_columns_with_header', capture_read)
    output, _ = _export(
        tmp_path,
        [{'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 20, 'unused': 'preserve', 'empty': None}],
        [{'x': 3, 'y': 4, 'z': 0, 'sample_UG': 'A', 'RQD': 20, 'unused': 'do not load', 'keep': 'yes'}],
        sample_filters=[{'field': 'keep', 'type': 'categorical', 'values': ['yes']}],
        progress_callback=lambda value, maximum, message: events.append((value, maximum, message)),
    )
    assert any(0 < value <= 20 and 'Reading samples' in message for value, _, message in events)
    assert any(20 < value <= 35 and 'Reading blocks' in message for value, _, message in events)
    assert any('Applying sample filters' in message for _, _, message in events)
    assert any('Preparing' in message for _, _, message in events)
    assert all(first[0] <= second[0] for first, second in zip(events, events[1:]))
    assert reads[0] == ('samples.csv', ['x', 'y', 'z', 'sample_UG', 'RQD', 'keep'], 'Reading samples')
    assert reads[1] == ('blocks.csv', ['x', 'y', 'z', 'UG', 'RQDpc', 'unused', 'empty'], 'Reading blocks')
    assert output.loc[0, 'unused'] == 'preserve'
    assert output['empty'].isna().all()
    assert output.loc[0, 'RQDpc_UG_MinDist'] == 5


@pytest.mark.parametrize('series', [
    pd.Series([20.0, 20.0, 30.0, np.nan], name='value'),
    pd.Series([' 20.0 ', '20', 'A', 'a', '', None, 'A'], name='value'),
    pd.Series([True, 1, 1.0, False, 0, '1', None], dtype=object, name='value'),
    pd.Series([1, 2, None, 1], dtype='Int64', name='value'),
    pd.Series(['A', None, 'B', 'A'], dtype='category', name='value'),
])
def test_distinct_key_normalization_preserves_scalar_matching(series):
    expected = series.map(viewer_module._normalize_categorical_filter_token)
    # Categorical.map can retain categorical dtype; matching keys use object dtype.
    pd.testing.assert_series_equal(
        viewer_module._normalize_matching_value_keys(series), expected.astype(object),
    )


def test_duplicate_block_attributes_are_not_lost(tmp_path):
    blocks_path, samples_path = tmp_path / 'blocks.csv', tmp_path / 'samples.csv'
    blocks_path.write_text('x,y,z,UG,RQDpc,attr,attr\n0,0,0,A,20,first,second\n')
    samples_path.write_text('x,y,z,sample_UG,RQD\n3,4,0,A,20\n')
    result = export_matching_value_sample_distances(
        str(samples_path), str(blocks_path), [{'block_col': 'RQDpc', 'sample_col': 'RQD'}],
        'sample_UG', 'UG',
    )
    output = pd.read_csv(result['output_file'])
    assert output.loc[0, 'attr'] == 'first'
    assert output.loc[0, 'attr_1'] == 'second'
    assert output.loc[0, 'RQDpc_UG_MinDist'] == 5


def test_quoted_block_header_with_39_columns(tmp_path):
    blocks_path, samples_path = tmp_path / 'blocks.csv', tmp_path / 'samples.csv'
    block = {'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 20}
    block.update({f'attribute{i}': i for i in range(33)})
    block['description, units'] = 'preserve, this'
    assert len(block) == 39
    pd.DataFrame([block]).to_csv(blocks_path, index=False)
    pd.DataFrame([{'x': 3, 'y': 4, 'z': 0, 'sample_UG': 'A', 'RQD': 20}]).to_csv(samples_path, index=False)
    result = export_matching_value_sample_distances(
        str(samples_path), str(blocks_path), [{'block_col': 'RQDpc', 'sample_col': 'RQD'}],
        'sample_UG', 'UG',
    )
    output = pd.read_csv(result['output_file'])
    assert output.columns.tolist() == list(block) + ['RQDpc_UG_MinDist', 'RQDpc_AbsMinDist']
    assert output.loc[0, 'description, units'] == 'preserve, this'
    assert output.loc[0, 'RQDpc_UG_MinDist'] == 5


@pytest.mark.parametrize('delimiter', [',', ';', '\t', '|'])
def test_quoted_sample_headers_and_blank_block_headers(tmp_path, delimiter):
    blocks_path, samples_path = tmp_path / 'blocks.csv', tmp_path / 'samples.csv'
    quoted_name = f'RQD{delimiter}percent "sample"'
    pd.DataFrame([{'x': 3, 'y': 4, 'z': 0, 'sample_UG': 'A', quoted_name: 20}]).to_csv(
        samples_path, index=False, sep=delimiter,
    )
    blocks_path.write_text(
        '# metadata\n' + delimiter.join(['x', 'y', 'z', 'UG', 'RQDpc', '', 'attribute']) + '\n'
        + delimiter.join(['0', '0', '0', 'A', '20', 'unnamed value', 'last value']) + '\n',
    )
    result = export_matching_value_sample_distances(
        str(samples_path), str(blocks_path), [{'block_col': 'RQDpc', 'sample_col': quoted_name}],
        'sample_UG', 'UG', samples_delimiter=delimiter, blocks_delimiter=delimiter,
    )
    output = pd.read_csv(result['output_file'])
    assert output.loc[0, 'Unnamed'] == 'unnamed value'
    assert output.loc[0, 'attribute'] == 'last value'
    assert output.loc[0, 'RQDpc_UG_MinDist'] == 5


@pytest.mark.parametrize('override, message', [
    ({'sample_domain_col': '(None)'}, 'explicit sample domain'),
    ({'block_domain_col': None}, 'explicit block domain'),
    ({'field_mappings': []}, 'at least one'),
    ({'field_mappings': [{'block_col': 'RQDpc', 'sample_col': '(None)'}]}, 'every mapping'),
    ({'field_mappings': [{'block_col': 'RQDpc', 'sample_col': 'missing'}]}, 'not found'),
    ({'field_mappings': [{'block_col': 'missing', 'sample_col': 'RQD'}]}, 'not found'),
    ({'field_mappings': [{'block_col': 'RQDpc', 'sample_col': 'RQD'}] * 2}, 'more than once'),
])
def test_invalid_selections_raise(tmp_path, override, message):
    blocks_path, samples_path = tmp_path / 'blocks.csv', tmp_path / 'samples.csv'
    pd.DataFrame([{'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 20}]).to_csv(blocks_path, index=False)
    pd.DataFrame([{'x': 0, 'y': 0, 'z': 0, 'sample_UG': 'A', 'RQD': 20}]).to_csv(samples_path, index=False)
    kwargs = dict(
        samples_file=str(samples_path), blocks_file=str(blocks_path),
        field_mappings=[{'block_col': 'RQDpc', 'sample_col': 'RQD'}],
        sample_domain_col='sample_UG', block_domain_col='UG',
    )
    kwargs.update(override)
    with pytest.raises(ValueError, match=message):
        export_matching_value_sample_distances(**kwargs)
    with pytest.raises(ValueError, match='different from the input'):
        export_matching_value_sample_distances(
            str(samples_path), str(blocks_path), [{'block_col': 'RQDpc', 'sample_col': 'RQD'}],
            'sample_UG', 'UG', output_file=str(blocks_path),
        )


def test_gui_configuration_and_operation_wiring(tmp_path):
    blocks_path, samples_path = tmp_path / 'blocks.csv', tmp_path / 'samples.csv'
    pd.DataFrame([{'x': 0, 'y': 0, 'z': 0, 'UG': 'A', 'RQDpc': 20}]).to_csv(blocks_path, index=False)
    pd.DataFrame([{'x': 3, 'y': 4, 'z': 0, 'sample_UG': 'A', 'RQD': 20}]).to_csv(samples_path, index=False)
    script = r'''
import os, runpy, sys
from PyQt5 import QtWidgets
import pandas as pd
sys.path.insert(0, sys.argv[1])
blocks, samples, output = sys.argv[2:5]
def exercise(app):
    dialog = next(widget for widget in app.topLevelWidgets() if type(widget).__name__ == 'ConfigDialog')
    config = dialog.to_dict()
    config.update(
        samples_file=samples, blocks_file=blocks, sample_domain_col='sample_UG', block_domain_col='UG',
        sample_x_col='x', sample_y_col='y', sample_z_col='z', block_x_col='x', block_y_col='y', block_z_col='z',
        matching_distance_file=output,
        matching_distance_field_mappings=[{'block_col': 'RQDpc', 'sample_col': 'RQD'}],
    )
    dialog.from_dict(config)
    restored = dialog.to_dict()
    assert restored['matching_distance_field_mappings'] == config['matching_distance_field_mappings']
    assert restored['matching_distance_file'] == output
    dialog._refresh_matching_distance_fields()
    assert dialog.to_dict()['matching_distance_field_mappings'] == config['matching_distance_field_mappings']
    dialog._set_operation_buttons_enabled(False)
    assert not dialog.start_matching_distance_btn.isEnabled()
    dialog._set_operation_buttons_enabled(True)
    def execute(title, message, operation, kwargs, success_handler, error_context):
        result = operation(**kwargs)
        assert pd.read_csv(result['output_file']).loc[0, 'RQDpc_UG_MinDist'] == 5
        success_handler(result)
    dialog._run_operation_with_progress = execute
    notifications = []
    QtWidgets.QMessageBox.information = lambda *args: notifications.append(args)
    QtWidgets.QMessageBox.critical = lambda *args: (_ for _ in ()).throw(AssertionError(args))
    dialog.run_matching_value_sample_distances_only()
    assert notifications
    dialog.matching_distance_fields.selectRow(0)
    dialog._remove_matching_distance_fields()
    assert dialog.to_dict()['matching_distance_field_mappings'] == []
    dialog._add_matching_distance_field(config['matching_distance_field_mappings'][0])
    assert dialog.to_dict()['matching_distance_field_mappings'] == config['matching_distance_field_mappings']
    dialog.close()
    return 0
QtWidgets.QApplication.exec_ = exercise
sys.argv = [os.path.join(sys.argv[1], 'anterpolator3DViewer.py')]
runpy.run_path(sys.argv[0], run_name='__main__')
'''
    result = subprocess.run(
        [sys.executable, '-c', script, PROJECT_ROOT, str(blocks_path), str(samples_path), str(tmp_path / 'output.csv')],
        capture_output=True, text=True, timeout=60, env={**os.environ, 'QT_QPA_PLATFORM': 'offscreen'},
    )
    assert result.returncode == 0, result.stdout + result.stderr
