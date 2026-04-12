Transparent Taichi Blocks

Status

Taichi volume rendering is deprecated. The Taichi viewer now runs in mesh mode only, and any saved config that still requests `volume` or `taichi_transparent_blocks: true` is normalized back to mesh.

Current behavior

The Display tab keeps only mesh-relevant controls. Sample display size is configured with `taichi_sample_diameter`, expressed in model units.

```json
{
  "taichi_block_render_mode": "mesh",
  "taichi_transparent_blocks": false,
  "taichi_sample_diameter": 1.0
}
```

Viewer controls

The viewer dialog still uses explicit buttons:

- `Open Viewer`: starts the Taichi viewer and performs the full initial load.
- `Refresh Viewer`: applies display-oriented changes to the running viewer without rebuilding sample/block data.
- `Reload Data`: re-reads input data, rebuilds blocks/domaining, and refreshes the running viewer.

Compatibility

Older config files may still contain `taichi_volume_*` keys from the deprecated renderer. Those values are ignored by the mesh-only viewer and can be removed by saving the config again.