# Porting Guide: Viewer Loading And Startup Improvements

This file is a manual rebuild guide for porting the loading and startup optimizations from this workspace into another workspace that already has its own interpolation changes.

Reference commit in this workspace:

- `dd411d0` - `Improve large file loading and viewer startup`

Primary source file in this workspace:

- `Anterpolator3DViewer/anterpolator3DViewer.py`

Do not blindly overwrite the target workspace file. Port these changes by responsibility area.

## What This Port Includes

This guide covers four change groups:

1. Loading-path changes for large CSV files.
2. True single-pass large-grid aggregation after rotation sampling.
3. Fast explicit-column loading for the samples file.
4. Deferred initial block rendering so the viewer window opens sooner.

## Recommended Port Order

Apply changes in this order:

1. Add the new helper functions and constants.
2. Replace the old rotation detection dependency path.
3. Replace the large `blocks_file` loading path.
4. Replace the sample-file loading path.
5. Add deferred block rendering.
6. Run the validation checklist at the end.

## 1. Add Constants And Lazy Imports

Add these module-level constants near the imports:

- `LARGE_BLOCK_FILE_THRESHOLD = 512 * 1024 * 1024`
- `INITIAL_BLOCK_RENDER_THRESHOLD = 5000`

Add lazy `pyvista` loading:

- Module global: `pv = None`
- Helper: `_require_pyvista()`

Intent:

- Avoid importing `pyvista` and `vtk` at module import time.
- Keep startup fast until the viewer actually needs 3D rendering.

Also change interpolator creation so `ant_colony` and `molecular_clock_interpolator` are imported inside `create_interpolator()` instead of at module scope.

## 2. Add CSV Progress And Selected-Column Helpers

Port these helpers as a group:

- `prepare_csv_read_kwargs(...)`
- `ProgressTextReader`
- `read_csv_with_progress(...)`
- `iterate_csv_with_progress(...)`
- `build_unique_column_names(...)`
- `read_selected_columns_with_header(...)`

Important implementation details:

- Progress must use the underlying binary file handle position via `self._raw.tell()`.
- Do not use the text wrapper `tell()` for progress; it returns an opaque cookie and produces nonsense values.
- For progress-wrapped reads, remove `memory_map` from `read_csv_kwargs`.
- For non-progress direct path reads, keep `memory_map=True` when using the pandas C engine on local files.

Intent:

- Keep the fast pandas parser.
- Show real progress on large files.
- Read only the columns you actually need when mappings are known.

## 3. Replace `sklearn` Rotation Detection With NumPy-Only Logic

Remove the dependency on `NearestNeighbors` inside `detect_grid_rotation(...)`.

Add helper:

- `_collect_same_level_vectors(...)`

Then update `detect_grid_rotation(...)` so it:

- Collects a bounded sample of points.
- Groups by approximate Z level.
- Computes 2D nearest neighbors within each level using NumPy distance matrices.
- Builds candidate XY vectors from those neighbors.
- Continues with the existing angle histogram logic.

Important behavior:

- Keep the return contract the same: `(rotation_matrix, center, is_rotated)`.
- Preserve the existing angle-selection logic downstream so the rest of the code does not need to change.

Intent:

- Avoid hanging on `sklearn.neighbors` import during large-file processing.
- Keep rotation detection cheap and self-contained.

## 4. Add Sample-File Fast Path

Add helper:

- `load_samples_dataframe(...)`

This helper should:

- Detect whether explicit sample column mapping exists.
- If explicit mapping exists, read only those four columns:
  - `sample_x_col`
  - `sample_y_col`
  - `sample_z_col`
  - `sample_value_col`
- Rename them immediately to `x`, `y`, `z`, `Value`.
- Fall back to the old general CSV loading path if mappings are not configured.

Replace the direct sample loading logic in both places:

1. `load_and_visualize_samples(...)`
2. `run_interpolation_only(...)`

Intent:

- Avoid parsing unrelated sample-file columns.
- Make sample loading behave more like the optimized block-file loading path.

## 5. Add Large-Grid Streaming Helpers

Port these helpers:

- `plan_block_file_columns(...)`
- `normalize_block_chunk(...)`
- `quantize_grid_indices(...)`
- `resolve_base_block_domains_from_counts(...)`
- `load_large_blocks_metadata(...)`

These are the core of the large-grid optimization.

### 5a. Rotation Sample Pass

Inside `load_large_blocks_metadata(...)`:

- Read only selected block columns.
- Use `chunksize=250_000`.
- Collect only enough valid rows for rotation detection.
- Stop early after `rotation_sample_target = 10000` rows.

Intent:

- The rotation angle does not need a full-file scan.

### 5b. True Single-Pass Bounds + Domain Mapping

This is the key algorithmic change.

Old behavior:

1. One full pass to compute final bounds.
2. Another full pass to compute domain mapping.

New behavior:

1. One early-stop rotation sample pass.
2. One combined full pass that computes both bounds and domain-domain counts.

Implementation idea:

- After rotation is known, pick a stable `grid_reference` from the first rotated row seen.
- Quantize every rotated coordinate row to an absolute integer grid key with:

  - `quantize_grid_indices(coords, grid_reference, unified_dims)`

- While streaming chunks:

  - Update `min_quantized_idx`
  - Update `max_quantized_idx`
  - Aggregate domain counts into a dictionary keyed by absolute grid index tuple

- At EOF:

  - Recover final `all_min_bounds` from `grid_reference + min_quantized_idx * unified_dims`
  - Recover `all_max_bounds` similarly from `max_quantized_idx`
  - Shift all absolute grid keys by subtracting `min_quantized_idx`
  - Convert those shifted counts into the final `domain_mapping`

Important note:

- This avoids rereading the source file for domain mapping.
- It works because the final grid indices differ from the absolute keys only by a constant shift.

## 6. Wire The Large-File Branch Into `create_blocks(...)`

Update `create_blocks(...)` so that when:

- `os.path.getsize(blocks_file) >= LARGE_BLOCK_FILE_THRESHOLD`

it uses `load_large_blocks_metadata(...)` instead of materializing the whole CSV DataFrame.

Keep the existing smaller-file code path intact.

The large-file branch should return at least these derived values back into `create_blocks(...)`:

- `all_min_bounds`
- `all_max_bounds`
- `unified_dims`
- `domain_mapping`
- `subblock_counts`
- `mixed_domain_blocks`
- `rotation_matrix`
- `rotation_center`
- `is_rotated`

Intent:

- Large files remain memory-safe.
- Smaller files keep the simpler in-memory behavior.

## 7. Add Deferred Initial Block Rendering

This change is separate from the CSV loading optimization. It solves the long pause after loading finishes but before the window appears.

Port these helpers:

- `_get_block_raw_value(...)`
- `_classify_block_value(...)`
- `_set_block_display_value(...)`
- `_get_block_grid_position(...)`
- `_build_block_lookup(...)`
- `_should_display_block(...)`
- `_create_visible_blocks(...)`
- `_get_blocks_mesh_kwargs(...)`
- `_prepare_blocks_for_display(...)`
- `_ensure_blocks_actor(...)`
- `_rebuild_blocks_actor(...)`
- `_mark_block_datasets_modified(...)`

Then update `load_and_visualize_samples(...)` so that:

- It sets these plotter fields initially to `None` or falsey values:
  - `_block_lookup`
  - `_visible_blocks`
  - `_visible_positions`
  - `_blocks_actor`
  - `_blocks_display_prepared`

- If `len(blocks) <= INITIAL_BLOCK_RENDER_THRESHOLD`, keep eager block rendering.
- Otherwise, do not call `add_mesh(...)` for blocks during startup.
- Instead print a message telling the user that pressing `b` will build/show block meshes.

Update `toggle_blocks(plotter)` so that:

- If `_blocks_actor` is missing, it calls `_ensure_blocks_actor(plotter, visible=True)`.
- Otherwise it performs the old visibility toggle behavior.

Intent:

- The viewer window opens sooner on very large models.
- Heavy block mesh construction becomes user-triggered instead of startup-blocking.

## 8. Keep `update_interpolation(...)` Compatible With Lazy Rendering

Because blocks may not be rendered yet, `update_interpolation(...)` must tolerate:

- missing `_block_lookup`
- missing `_visible_blocks`
- missing `_visible_positions`
- missing `_blocks_actor`

Port the guard logic that initializes these lazily before using them.

Also keep the block-update logic that:

- updates existing blocks in place when possible
- appends new blocks only when needed
- rebuilds the actor only when visibility membership changed

## 9. Minimal Search Anchors In The Target Workspace

In the target workspace, search for these symbols or messages and patch around them:

- `create_interpolator`
- `detect_grid_rotation`
- `create_blocks`
- `load_and_visualize_samples`
- `run_interpolation_only`
- `toggle_blocks`
- `update_interpolation`
- `Loading predefined cells from blocks_file...`
- `Loading sample file from`
- `plotter.show()`

## 10. Expected Runtime Behavior After Port

For very large grid files, startup should look like this:

1. `Reading grid file (rotation sample)`
2. `Detecting rotation on ... points...`
3. `Building bounds and domain mapping...`
4. `Reading grid file (bounds + domain mapping)`
5. Domain initialization logs
6. Viewer window appears without waiting for eager block mesh creation when block count is large

For mapped sample files, startup should show a single progress bar that reads only four columns instead of the full sample table.

## 11. Validation Checklist After Porting

Run these checks in the target workspace:

1. Import time:
   - `python -c "import anterpolator3DViewer"`
   - Confirm startup import is still fast.

2. Large blocks file:
   - Confirm the `rotation sample` pass stops early.
   - Confirm there is only one full pass after rotation sampling.

3. Sample file:
   - Confirm explicit mapped columns produce a DataFrame with only `x`, `y`, `z`, `Value`.

4. Viewer startup:
   - Confirm the window appears before block meshes are fully built on large datasets.
   - Press `b` to build/show blocks.

5. Interpolation update:
   - Press `i` and confirm updates still work with lazy-rendered block actors.

## 12. Practical Merge Advice

Because the target workspace already has interpolation changes, do not treat this as a line-by-line transplant.

Merge by intent:

1. Port helper functions first.
2. Port the large-file control flow in `create_blocks(...)`.
3. Port the sample-file helper usage.
4. Port the deferred-rendering helpers and wiring last.

If the target workspace diverges heavily, prefer copying the helper functions verbatim and then manually reconnecting the call sites rather than trying to replay a raw diff.