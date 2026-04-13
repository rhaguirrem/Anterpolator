# Large Dataset Handling Reference

This document explains how this workspace handles large datasets so the same patterns can be taught or ported into other projects.

The main implementation lives in [../Anterpolator3DViewer/anterpolator3DViewer.py](../Anterpolator3DViewer/anterpolator3DViewer.py). The existing porting note in [PORT_VIEWER_LOADING_IMPROVEMENTS.md](PORT_VIEWER_LOADING_IMPROVEMENTS.md) is useful if you want a change-by-change transplant checklist; this file is meant to explain the design and point to the most relevant code.

## Executive Summary

This project stays usable with large CSV-backed models by combining five tactics:

1. It avoids expensive imports until they are actually needed.
2. It reads only the columns required for the current operation.
3. It streams very large block files in chunks instead of materializing them entirely.
4. It computes grid metadata in a sampling pass plus one full streaming pass, rather than multiple full rereads.
5. It defers heavy 3D block actor creation until the user needs it.

In practice, most of the large-dataset behavior is concentrated in a single module, which makes this workspace a good reference implementation.

## Where The Behavior Lives

Primary file:

- [../Anterpolator3DViewer/anterpolator3DViewer.py](../Anterpolator3DViewer/anterpolator3DViewer.py)

Supporting documentation already in this repo:

- [PORT_VIEWER_LOADING_IMPROVEMENTS.md](PORT_VIEWER_LOADING_IMPROVEMENTS.md)
- [../README.md](../README.md)

Repository notes that also describe the performance work:

- `/memories/repo/csv-performance.md`
- `/memories/repo/startup-performance.md`

## Architecture Pattern

If you want to reproduce this approach in another codebase, copy the pattern, not just the syntax:

1. Keep a cheap startup path.
2. Detect file structure once.
3. Read only the columns you need.
4. Stream large inputs in predictable chunk sizes.
5. Sample early for any heuristic that does not require the full dataset.
6. Aggregate metadata while streaming so the file does not need to be scanned again.
7. Postpone visualization work that can wait until after the window is open.

## Core Techniques

### 1. Lazy imports and startup control

The module avoids paying the cost of heavy visualization and interpolator imports until the code path really needs them.

Relevant entry points:

- [create_interpolator](../Anterpolator3DViewer/anterpolator3DViewer.py#L73)
- [prepare_csv_read_kwargs](../Anterpolator3DViewer/anterpolator3DViewer.py#L202)
- [_require_pyvista](../Anterpolator3DViewer/anterpolator3DViewer.py#L28)

What to notice:

- `pv` starts as `None`, and `_require_pyvista()` imports `pyvista` only on demand.
- `create_interpolator()` imports algorithm implementations inside the factory instead of at module import time.
- This keeps the config dialog and basic startup responsive even when optional dependencies are heavy.

Why this matters for large datasets:

- Large file handling feels worse if startup is already blocked by imports. This project fixes that first.

### 2. Fast CSV reads with the pandas C engine when possible

The loader keeps the default fast parser whenever the delimiter permits it and enables memory mapping for direct file-path reads.

Relevant functions:

- [prepare_csv_read_kwargs](../Anterpolator3DViewer/anterpolator3DViewer.py#L202)
- [read_csv_with_progress](../Anterpolator3DViewer/anterpolator3DViewer.py#L300)
- [iterate_csv_with_progress](../Anterpolator3DViewer/anterpolator3DViewer.py#L310)

What the code does:

- Uses pandas' default engine unless a complex delimiter forces `engine='python'`.
- Enables `memory_map=True` for local path-based reads when the fast engine is still available.
- Removes `memory_map` for wrapped progress reads because those operate on a file handle, not a plain path.

Transferable rule:

- Prefer the fastest parser by default, and fall back only when the input format demands it.

### 3. Real progress reporting without breaking parsing

The progress wrapper is careful about how bytes are counted.

Relevant functions:

- [ProgressTextReader](../Anterpolator3DViewer/anterpolator3DViewer.py#L212)
- [read_csv_with_progress](../Anterpolator3DViewer/anterpolator3DViewer.py#L300)
- [iterate_csv_with_progress](../Anterpolator3DViewer/anterpolator3DViewer.py#L310)

Why it works:

- `ProgressTextReader` opens the file in binary mode, wraps it with `io.TextIOWrapper`, and tracks bytes through the raw handle.
- Progress is read from `self._raw.tell()`, not the text wrapper position.
- Updates are throttled so UI feedback does not dominate the parse time.

This is a good pattern when another project needs both throughput and trustworthy progress for multi-hundred-megabyte or multi-gigabyte text files.

### 4. Selected-column sample loading

The samples file has a fast path when the user already mapped the `x`, `y`, `z`, and `Value` columns.

Relevant functions:

- [read_selected_columns_with_header](../Anterpolator3DViewer/anterpolator3DViewer.py#L337)
- [load_samples_dataframe](../Anterpolator3DViewer/anterpolator3DViewer.py#L358)
- [load_and_visualize_samples](../Anterpolator3DViewer/anterpolator3DViewer.py#L2756)

What the code does:

- Parses the configured header line.
- Builds stable unique column names.
- Uses `usecols` so pandas reads only the four required sample columns.
- Renames the chosen columns immediately to the normalized internal schema: `x`, `y`, `z`, `Value`.

Why this matters:

- Real-world assay files often carry many auxiliary columns. Parsing all of them wastes memory and CPU when interpolation only needs four fields.

### 5. Domain-only streaming for UI preparation

The configuration dialog can preload just the domain catalog from the block file and cache it, instead of scanning the whole file repeatedly.

Relevant functions:

- [build_domain_catalog_cache_signature](../Anterpolator3DViewer/anterpolator3DViewer.py#L394)
- [load_block_domain_catalog](../Anterpolator3DViewer/anterpolator3DViewer.py#L405)
- [domain catalog loading in the dialog](../Anterpolator3DViewer/anterpolator3DViewer.py#L4103)

What the code does:

- Reads only the selected domain column in chunks.
- Reports progress into a Qt progress dialog.
- Caches the discovered domain list keyed by file path, size, mtime, delimiter, header line, and chosen domain column.

Why this matters:

- In large-data applications, metadata discovery can become an accidental bottleneck. This keeps the UI from redoing expensive scans unnecessarily.

### 6. NumPy-only rotation detection on a bounded sample

The project avoids a heavy dependency path for grid rotation detection and only inspects a bounded set of points.

Relevant functions:

- [_collect_same_level_vectors](../Anterpolator3DViewer/anterpolator3DViewer.py#L546)
- [detect_grid_rotation](../Anterpolator3DViewer/anterpolator3DViewer.py#L601)

What the code does:

- Samples points instead of analyzing the whole grid.
- Groups points by approximate Z level.
- Builds candidate XY neighbor vectors with NumPy operations.
- Preserves the existing rotation result contract: rotation matrix, center, and a rotated/not-rotated flag.

Why this matters:

- Rotation detection is a heuristic. It does not need the full file, so the project treats it as a sampling problem rather than a full-data problem.

### 7. Chunked large-grid streaming instead of full DataFrame materialization

This is the main large-dataset technique in the codebase.

Relevant functions:

- [plan_block_file_columns](../Anterpolator3DViewer/anterpolator3DViewer.py#L769)
- [normalize_block_chunk](../Anterpolator3DViewer/anterpolator3DViewer.py#L797)
- [quantize_grid_indices](../Anterpolator3DViewer/anterpolator3DViewer.py#L437)
- [resolve_base_block_domains_from_counts](../Anterpolator3DViewer/anterpolator3DViewer.py#L736)
- [load_large_blocks_metadata](../Anterpolator3DViewer/anterpolator3DViewer.py#L812)
- [create_blocks](../Anterpolator3DViewer/anterpolator3DViewer.py#L987)

How the algorithm works:

1. If the block file is below the threshold, the code keeps the simpler in-memory path.
2. If the file size reaches `LARGE_BLOCK_FILE_THRESHOLD`, the code switches to streaming mode.
3. The first pass is intentionally small: it reads chunked data only until there are enough valid coordinate rows to estimate rotation.
4. The second pass is the real work: it streams the file chunk by chunk and computes bounds plus domain counts together.

The important idea is in [load_large_blocks_metadata](../Anterpolator3DViewer/anterpolator3DViewer.py#L812):

- Coordinates are normalized and optionally rotated.
- Each row is quantized to an absolute integer grid index relative to a stable `grid_reference`.
- The code tracks running minimum and maximum quantized indices.
- Domain counts are aggregated per absolute base-cell index while the file is still streaming.
- At the end, indices are shifted once to recover final grid-local coordinates and build the final `domain_mapping`.

Why this is effective:

- The code avoids keeping the whole block CSV in memory.
- It also avoids a second full reread just to compute domain mapping after bounds are known.

### 8. Threshold-based branching between small and large file strategies

The code does not force the complex path on every input.

Relevant constants and entry point:

- [LARGE_BLOCK_FILE_THRESHOLD](../Anterpolator3DViewer/anterpolator3DViewer.py#L25)
- [create_blocks](../Anterpolator3DViewer/anterpolator3DViewer.py#L987)

Behavior:

- Small files use the straightforward in-memory normalization and block-building path.
- Large files switch to the streamed metadata path.

This is worth copying into other projects because it keeps maintenance simple: optimize only where scale forces it.

### 9. Deferred initial 3D block rendering

Reading the data is only half the problem. Rendering thousands of blocks can also stall the application.

Relevant functions:

- [_build_block_lookup](../Anterpolator3DViewer/anterpolator3DViewer.py#L1838)
- [_create_visible_blocks](../Anterpolator3DViewer/anterpolator3DViewer.py#L1850)
- [_prepare_blocks_for_display](../Anterpolator3DViewer/anterpolator3DViewer.py#L1888)
- [_ensure_blocks_actor](../Anterpolator3DViewer/anterpolator3DViewer.py#L1909)
- [_rebuild_blocks_actor](../Anterpolator3DViewer/anterpolator3DViewer.py#L1923)
- [load_and_visualize_samples](../Anterpolator3DViewer/anterpolator3DViewer.py#L2756)
- [deferred-render trigger](../Anterpolator3DViewer/anterpolator3DViewer.py#L3057)

What the code does:

- Initializes block-display state on the plotter without immediately building the actor.
- If the block count is below `INITIAL_BLOCK_RENDER_THRESHOLD`, it renders eagerly.
- If the block count is larger, it defers block actor creation and lets the window open first.
- The user can then press `b` to build or show blocks on demand.

Why this matters:

- Large-data applications often solve I/O performance and then still feel slow because first render is too expensive. This project addresses both.

## The Main Data Flow

For another project, the reusable pipeline looks like this:

1. Start cheap: lazy-import the heavy visualization stack.
2. Detect or receive delimiter and header configuration.
3. Read only required columns.
4. For large block files, stream chunks of normalized rows.
5. Run a small early sampling pass for any geometric heuristics.
6. Aggregate bounds and domain statistics during the same full streaming pass.
7. Build the lightweight viewer state first.
8. Create the heavy visual actor only when necessary.

## Relevant Functions By Responsibility

### Startup and dependency control

- [create_interpolator](../Anterpolator3DViewer/anterpolator3DViewer.py#L73)
- [_require_pyvista](../Anterpolator3DViewer/anterpolator3DViewer.py#L28)

### CSV parsing and progress

- [prepare_csv_read_kwargs](../Anterpolator3DViewer/anterpolator3DViewer.py#L202)
- [ProgressTextReader](../Anterpolator3DViewer/anterpolator3DViewer.py#L212)
- [read_csv_with_progress](../Anterpolator3DViewer/anterpolator3DViewer.py#L300)
- [iterate_csv_with_progress](../Anterpolator3DViewer/anterpolator3DViewer.py#L310)
- [read_selected_columns_with_header](../Anterpolator3DViewer/anterpolator3DViewer.py#L337)

### Sample-file optimization

- [load_samples_dataframe](../Anterpolator3DViewer/anterpolator3DViewer.py#L358)
- [load_and_visualize_samples](../Anterpolator3DViewer/anterpolator3DViewer.py#L2756)

### Block-file streaming and aggregation

- [quantize_grid_indices](../Anterpolator3DViewer/anterpolator3DViewer.py#L437)
- [_collect_same_level_vectors](../Anterpolator3DViewer/anterpolator3DViewer.py#L546)
- [detect_grid_rotation](../Anterpolator3DViewer/anterpolator3DViewer.py#L601)
- [resolve_base_block_domains_from_counts](../Anterpolator3DViewer/anterpolator3DViewer.py#L736)
- [plan_block_file_columns](../Anterpolator3DViewer/anterpolator3DViewer.py#L769)
- [normalize_block_chunk](../Anterpolator3DViewer/anterpolator3DViewer.py#L797)
- [load_large_blocks_metadata](../Anterpolator3DViewer/anterpolator3DViewer.py#L812)
- [create_blocks](../Anterpolator3DViewer/anterpolator3DViewer.py#L987)

### Domain metadata preparation and caching

- [build_domain_catalog_cache_signature](../Anterpolator3DViewer/anterpolator3DViewer.py#L394)
- [load_block_domain_catalog](../Anterpolator3DViewer/anterpolator3DViewer.py#L405)
- [dialog-side cache usage](../Anterpolator3DViewer/anterpolator3DViewer.py#L4103)

### Deferred rendering

- [_get_block_raw_value](../Anterpolator3DViewer/anterpolator3DViewer.py#L1804)
- [_classify_block_value](../Anterpolator3DViewer/anterpolator3DViewer.py#L1809)
- [_set_block_display_value](../Anterpolator3DViewer/anterpolator3DViewer.py#L1821)
- [_get_block_grid_position](../Anterpolator3DViewer/anterpolator3DViewer.py#L1834)
- [_build_block_lookup](../Anterpolator3DViewer/anterpolator3DViewer.py#L1838)
- [_should_display_block](../Anterpolator3DViewer/anterpolator3DViewer.py#L1845)
- [_create_visible_blocks](../Anterpolator3DViewer/anterpolator3DViewer.py#L1850)
- [_get_blocks_mesh_kwargs](../Anterpolator3DViewer/anterpolator3DViewer.py#L1864)
- [_prepare_blocks_for_display](../Anterpolator3DViewer/anterpolator3DViewer.py#L1888)
- [_ensure_blocks_actor](../Anterpolator3DViewer/anterpolator3DViewer.py#L1909)
- [_rebuild_blocks_actor](../Anterpolator3DViewer/anterpolator3DViewer.py#L1923)
- [_mark_block_datasets_modified](../Anterpolator3DViewer/anterpolator3DViewer.py#L1933)

## Porting Rules For Other Projects

If another codebase needs to learn from this one, these are the rules worth preserving:

1. Keep the normal path simple for small files.
2. Gate the streaming path behind a file-size threshold.
3. Normalize input columns early into a small internal schema.
4. Sample early for heuristics such as rotation detection.
5. Use one full streaming pass to accumulate as much metadata as possible.
6. Separate data loading from actor creation or other expensive rendering work.
7. Cache metadata scans that users are likely to repeat from the UI.

## When This Pattern Is A Good Fit

This design is a good reference for projects that have:

- CSV or delimited files large enough that full DataFrame materialization is risky.
- Geometric or spatial grids that can be represented with quantized integer indices.
- UI workflows where users need progress feedback during long reads.
- Visualization layers where opening the window quickly is more important than eagerly drawing everything.

It is less appropriate when:

- The data source is already a database or columnar format with native predicate pushdown.
- The full dataset must be in memory for downstream algorithms anyway.
- Rendering cost is negligible compared with computation.

## Short Porting Checklist

If you are teaching another project to do this, the minimum useful transplant is:

1. Copy the lazy-import pattern.
2. Copy the selected-column CSV read path.
3. Copy the chunked large-block metadata pass.
4. Copy the threshold switch in [create_blocks](../Anterpolator3DViewer/anterpolator3DViewer.py#L987).
5. Copy the deferred-rendering pattern in [load_and_visualize_samples](../Anterpolator3DViewer/anterpolator3DViewer.py#L2756).
6. Add the domain-catalog cache if the target project has a similar UI.

That combination is what makes this workspace feel capable of handling large datasets rather than merely being able to parse them.