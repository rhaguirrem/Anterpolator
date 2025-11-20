Ant Colony Fix (2025-09-17)

Summary
- Restored stable weighting behavior from the 2025-04-09 backup while preserving new features (background fill, distance cutoff, average_with_blocks, domain support, and ant backtracking avoidance).

Changes
- calculate_block_value: added nearest-sample proximity factor and ensured deterministic float outputs (avoid numpy scalar leakage). This mirrors logic in backup and improves stability.
- Left movement logic and extended API as-is to remain compatible with current 3D viewer.

Notes
- File compared against backup/ant_colony_20250409.py and selectively merged improvements rather than a full rollback, because the current viewer relies on newer parameters and behavior.

Quick test
A small smoke test was run in the workspace to import and step the interpolator; no syntax/runtime errors observed and values are plain floats.
