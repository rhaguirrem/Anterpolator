import ast
import os
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import anterpolator3DViewer as app
from ant_colony import AntColonyInterpolator, Block
from ant_colony_ii import AntColonyIIInterpolator
from block_postprocess import fill_isolated_blank_blocks
from provenance_utils import get_export_provenance


def make_model(interpolator_type=AntColonyInterpolator, target='value'):
    positions = {(0, 0, 0), (0, 2, 0), (4, 0, 0), (5, 0, 0), (6, 0, 0)}
    domains = {pos: 'A' for pos in positions}
    domains[(6, 0, 0)] = 'B'
    interp = interpolator_type(interpolation_target=target)
    interp.allowed_grid_override = positions
    interp.domain_mapping = domains
    interp.initialize_blocks(
        {(0, 0, 0): 10.0, (0, 2, 0): 30.0}, (7, 3, 1),
        np.zeros(3), [2.0, 1.0, 1.0], use_domain_mapping=True,
        sample_domain_mapping={(0, 0, 0): 'A', (0, 2, 0): 'A'},
    )
    return SimpleNamespace(
        _ant_colony=interp,
        _block_info={
            'min_bounds': np.zeros(3), 'dims': (7, 3, 1),
            'block_size': [2.0, 1.0, 1.0],
            'allowed_grid': list(positions), 'domain_mapping': domains,
        },
    )


@pytest.mark.parametrize('interpolator_type', [AntColonyInterpolator, AntColonyIIInterpolator])
@pytest.mark.parametrize('n', [1, 2, 8])
def test_isolated_fill_uses_fixed_same_domain_donors_and_world_distances(interpolator_type, n):
    model = make_model(interpolator_type)
    assert fill_isolated_blank_blocks(
        model, {'fill_isolated_blanks': True, 'isolated_fill_neighbors': n},
    ) == 2
    interp = model._ant_colony
    for pos in [(4, 0, 0), (5, 0, 0)]:
        if n == 1:
            expected = 10.0
        else:
            distance_1 = pos[0] * 2.0
            distance_2 = np.hypot(distance_1, 2.0)
            expected = (10.0 / distance_1**2 + 30.0 / distance_2**2) / (1 / distance_1**2 + 1 / distance_2**2)
        assert interp.blocks[pos].value == pytest.approx(expected)
        assert interp.blocks[pos].domain == 'A'
        assert get_export_provenance(interp, pos) == (
            'Post-process', 'Same-domain IDW Fill', 'same_domain_idw_fill',
        )
    assert interp.blocks[(0, 0, 0)].value == 10.0
    assert interp.blocks[(0, 2, 0)].value == 30.0
    assert (6, 0, 0) not in interp.blocks


def test_isolated_fill_is_optional_and_does_not_change_domain_interpolation():
    model = make_model()
    assert fill_isolated_blank_blocks(model, {}) == 0
    domain_model = make_model(target='domain')
    assert fill_isolated_blank_blocks(domain_model, {'fill_isolated_blanks': True}) == 0
    assert len(domain_model._ant_colony.blocks) == 2
    with pytest.raises(ValueError, match='positive integer'):
        fill_isolated_blank_blocks(model, {'fill_isolated_blanks': True, 'isolated_fill_neighbors': 0})


def test_isolated_fill_replaces_nonfinite_blanks_without_cross_domain_donors():
    model = make_model()
    interp = model._ant_colony
    interp.blocks[(4, 0, 0)] = Block(value=np.nan, is_sample=False, mark_class=0.0, domain='A')
    interp.blocks[(6, 0, 0)] = Block(value=9999.0, is_sample=True, mark_class=9999.0, domain='B')
    fill_isolated_blank_blocks(model, {'fill_isolated_blanks': True, 'isolated_fill_neighbors': 1})
    assert interp.blocks[(4, 0, 0)].value == 10.0
    assert interp.blocks[(5, 0, 0)].value == 10.0
    assert interp.blocks[(6, 0, 0)].value == 9999.0
    assert get_export_provenance(interp, (6, 0, 0)) is None


@pytest.mark.parametrize('streamed', [False, True])
@pytest.mark.parametrize('target', ['value', 'domain'])
@pytest.mark.parametrize('empty', [False, True])
def test_export_retains_regular_source_rows_and_value_mode_source_domains(tmp_path, monkeypatch, streamed, target, empty):
    model = make_model(target=target)
    source_path = tmp_path / 'source.csv'
    output_path = tmp_path / 'output.csv'
    source = pd.DataFrame([
        {'X': 8.0, 'Y': 0.0, 'Z': 0.0, 'UG': 'A'},
        {'X': 0.0, 'Y': 0.0, 'Z': 0.0, 'UG': 'original'},
        {'X': 12.0, 'Y': 0.0, 'Z': 0.0, 'UG': 'B'},
    ])
    source.to_csv(source_path, index=False)
    model._block_info.update({
        'source_blocks_file': str(source_path), 'source_blocks_delimiter': ',',
        'source_blocks_header_line': 1, 'source_block_x_col': 'X',
        'source_block_y_col': 'Y', 'source_block_z_col': 'Z',
        'source_block_domain_col': 'UG',
        'expand_interpolation_exports_to_subblocks': False,
    })
    if empty:
        model._ant_colony.blocks.clear()
    monkeypatch.setattr(app, 'LARGE_BLOCK_FILE_THRESHOLD', 1 if streamed else 10**9)
    app.export_blocks_to_csv(model, str(output_path), value_column_name='RQD_pctg')
    exported = pd.read_csv(output_path)
    assert len(exported) == len(source)
    np.testing.assert_array_equal(exported[['x', 'y', 'z']], source[['X', 'Y', 'Z']])
    assert exported.loc[[0, 2], 'RQD_pctg'].isna().all()
    if target == 'value':
        assert exported.Domain.tolist() == source.UG.tolist()
    elif not empty:
        assert exported.loc[1, 'Domain'] == 'A'
        assert exported.loc[[0, 2], 'Domain'].isna().all()
    assert 'interpolation_value_column_name' not in model._block_info


@pytest.mark.parametrize('streamed', [False, True])
def test_export_preserves_each_subblock_domain_not_parent_majority(tmp_path, monkeypatch, streamed):
    model = make_model()
    source_path = tmp_path / 'subblocks.csv'
    output_path = tmp_path / 'out.csv'
    source = pd.DataFrame([
        {'x': 0.25, 'y': 0.0, 'z': 0.0, 'UG': 'A'},
        {'x': 0.75, 'y': 0.0, 'z': 0.0, 'UG': 'B'},
        {'x': 8.25, 'y': 0.0, 'z': 0.0, 'UG': 'C'},
    ])
    source.to_csv(source_path, index=False)
    model._block_info.update({
        'source_blocks_file': str(source_path), 'source_blocks_delimiter': ',',
        'source_block_x_col': 'x', 'source_block_y_col': 'y', 'source_block_z_col': 'z',
        'source_block_domain_col': 'UG',
    })
    monkeypatch.setattr(app, 'LARGE_BLOCK_FILE_THRESHOLD', 1 if streamed else 10**9)
    app.export_blocks_to_csv(model, str(output_path))
    output = pd.read_csv(output_path)
    assert output.x.tolist() == source.x.tolist()
    assert output.Domain.tolist() == source.UG.tolist()
    assert output.Value.iloc[:2].tolist() == [10.0, 10.0]
    assert pd.isna(output.loc[2, 'Value'])


def test_isolated_fill_follows_average_fill_and_includes_its_values():
    model = make_model()
    model._block_info['allowed_grid'].append((1, 0, 0))
    model._block_info['domain_mapping'][(1, 0, 0)] = 'A'
    interp = model._ant_colony
    interp.allowed_positions.add((1, 0, 0))
    created, assigned = app._run_domain_post_process(interp, (7, 3, 1), 'fill_with_average', 'A')
    assert created + assigned > 0
    average_value = interp.blocks[(1, 0, 0)].value
    assert (4, 0, 0) not in interp.blocks
    fill_isolated_blank_blocks(model, {'fill_isolated_blanks': True, 'isolated_fill_neighbors': 1})
    assert interp.blocks[(4, 0, 0)].value == pytest.approx(average_value)


def test_silent_interpolation_runs_isolated_fill_after_average_stage(tmp_path, monkeypatch):
    model = make_model()
    interp = model._ant_colony
    model._interpolators = {'A': [interp]}
    model._block_info['dims'] = (7, 3, 1)
    events = []
    monkeypatch.setattr(app, '_run_interpolator_with_progress', lambda *args, **kwargs: None)
    monkeypatch.setattr(app, '_run_interpolator_statistics_with_retry', lambda *args, **kwargs: None)
    monkeypatch.setattr(app, '_run_domain_post_process', lambda *args, **kwargs: events.append('average') or (0, 0))
    monkeypatch.setattr(app, 'fill_isolated_blank_blocks', lambda *args: events.append('idw') or 0)
    monkeypatch.setattr(app, 'export_blocks_to_file', lambda *args, **kwargs: events.append('export') or str(tmp_path / 'out.csv'))
    plotter = SimpleNamespace(_blocks_data=model, _fill_isolated_blanks=True)
    app.silent_interpolation(plotter, 1, str(tmp_path / 'out.csv'))
    assert events == ['average', 'idw', 'export']


def test_interactive_viewer_runs_isolated_fill_after_average_stage(monkeypatch):
    runtime = app._load_taichi_runtime_module()
    model = make_model()
    model._interpolators = {'A': [model._ant_colony]}
    viewer = runtime.TaichiInterpolationViewer.__new__(runtime.TaichiInterpolationViewer)
    viewer.blocks_data = model
    viewer.block_info = model._block_info
    viewer.config = {'iterations': 1, 'fill_isolated_blanks': True}
    viewer._reset_interpolation_progress()
    events = []
    monkeypatch.setattr(viewer, '_run_single_interpolator_iteration', lambda *args: (False, False, True))
    monkeypatch.setattr(viewer, '_run_domain_post_process', lambda *args: events.append('average') or (0, 0))
    monkeypatch.setattr(runtime, 'fill_isolated_blank_blocks', lambda *args: events.append('idw') or 1)
    result = viewer._run_interpolators_once()
    assert result['all_finished']
    assert result['changed']
    assert events == ['average', 'idw']


@pytest.mark.parametrize('algorithm', ['ant_colony', 'ant_colony_ii'])
def test_advanced_settings_roundtrip_and_gui_run_export(tmp_path, monkeypatch, algorithm):
    source_path = tmp_path / 'blocks.csv'
    sample_path = tmp_path / 'samples.csv'
    output_path = tmp_path / 'out.csv'
    pd.DataFrame([
        {'X': x, 'Y': 0.5, 'Z': 0.5, 'UG': 'A'} for x in [0.5, 1.5, 4.5]
    ]).to_csv(source_path, index=False)
    pd.DataFrame([
        {'sx': 0.5, 'sy': 0.5, 'sz': 0.5, 'grade': 10.0}
    ]).to_csv(sample_path, index=False)

    # Define the actual GUI classes without running the module's application loop.
    with open(app.__file__, encoding='utf-8-sig') as source_file:
        module_ast = ast.parse(source_file.read(), filename=app.__file__)
    main = next(node for node in module_ast.body if isinstance(node, ast.If) and ast.unparse(node.test) == "__name__ == '__main__'")
    classes = ast.Module(body=[node for node in main.body if isinstance(node, ast.ClassDef)], type_ignores=[])
    namespace = dict(vars(app))
    exec(compile(classes, app.__file__, 'exec'), namespace)
    monkeypatch.setenv('QT_QPA_PLATFORM', 'offscreen')
    qt_app = app.QtWidgets.QApplication.instance() or app.QtWidgets.QApplication([])
    errors = []
    monkeypatch.setattr(app.QtWidgets.QMessageBox, 'critical', lambda *args: errors.append(args[-1]))
    monkeypatch.setattr(app.QtWidgets.QMessageBox, 'information', lambda *args: None)
    dialog = namespace['ConfigDialog']()
    try:
        assert not dialog.fill_isolated_blanks.isChecked()
        assert dialog.isolated_fill_neighbors.value() == 8
        dialog.from_dict({
            'samples_file': str(sample_path), 'blocks_file': str(source_path),
            'interpolation_file': str(output_path), 'color_file': '',
            'samples_delimiter': ',', 'blocks_delimiter': ',',
            'samples_header_line': 1, 'blocks_header_line': 1,
            'sample_x_col': 'sx', 'sample_y_col': 'sy', 'sample_z_col': 'sz',
            'sample_value_col': 'grade', 'block_x_col': 'X', 'block_y_col': 'Y',
            'block_z_col': 'Z', 'block_domain_col': 'UG', 'block_size': [1, 1, 1],
            'algorithm': algorithm, 'iterations': 1,
            'process_domains_sequentially': True,
            'domain_algorithm_overrides': {'A': {'post_process': 'fill_with_average'}},
            'fill_isolated_blanks': True, 'isolated_fill_neighbors': 3,
        })
        saved = dialog.to_dict()
        assert saved['fill_isolated_blanks']
        assert saved['isolated_fill_neighbors'] == 3
        assert dialog.isolated_fill_neighbors.isEnabled()
        dialog.run_interpolation_only()
        assert errors == []
        output = pd.read_csv(output_path)
        assert output.x.tolist() == [0.5, 1.5, 4.5]
        assert output.Domain.tolist() == ['A', 'A', 'A']
        assert output.loc[0, 'grade'] == 10.0
        expected = (output.loc[0, 'grade'] / 4.0**2 + output.loc[1, 'grade'] / 3.0**2) / (1 / 4.0**2 + 1 / 3.0**2)
        assert output.loc[2, 'grade'] == pytest.approx(expected)
        assert output.loc[2, 'Algorithm'] == 'Same-domain IDW Fill'
    finally:
        dialog._prepare_for_shutdown()
        dialog.close()
        qt_app.processEvents()
