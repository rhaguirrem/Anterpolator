from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, Iterable, Optional, Tuple


_PROVENANCE_ATTR = '_export_provenance'


def normalize_provenance_algorithm_name(algo_name: Any) -> Tuple[str, str]:
    name = str(algo_name or '').strip().lower()
    if name in ('sample',):
        return 'Sample', 'sample'
    if 'fill with average' in name or 'fill_with_average' in name or 'fill average' in name:
        return 'Fill with Average', 'fill_with_average'
    if 'ant colony' in name or 'anterpolator' in name:
        return 'Ant Colony', 'ant_colony'
    if 'adaptive octree' in name:
        return 'Adaptive Octree', 'adaptive_octree'
    if 'string theory' in name:
        return 'String Theory', 'string_theory'
    if 'molecular clock' in name or 'phylogeographic' in name or 'biochemical clock' in name:
        return 'Molecular Clock', 'molecular_clock'
    if 'gaussian kernel' in name:
        return 'Gaussian Kernel', 'gaussian_kernel'
    return 'Unknown', 'unknown'


def _append_unique(target: list[str], value: Optional[str]):
    cleaned = str(value or '').strip()
    if cleaned and cleaned not in target:
        target.append(cleaned)


def _get_block_attr(block: Any, key: str, default: Any = None) -> Any:
    if hasattr(block, key):
        return getattr(block, key)
    if isinstance(block, dict):
        return block.get(key, default)
    return default


def _get_block_is_sample(block: Any) -> bool:
    return bool(_get_block_attr(block, 'is_sample', False))


def _normalize_scalar(value: Any) -> Any:
    if value is None:
        return None
    if hasattr(value, 'item'):
        try:
            value = value.item()
        except Exception:
            pass
    try:
        return float(value)
    except Exception:
        return value


def _current_target(interpolator: Any) -> str:
    target = str(getattr(interpolator, 'interpolation_target', 'value') or 'value').strip().lower()
    if target == 'domain':
        return 'domain'
    return 'value'


def get_interpolator_provenance(interpolator: Any) -> Dict[tuple, Dict[str, list[str]]]:
    existing = getattr(interpolator, _PROVENANCE_ATTR, None)
    if isinstance(existing, dict):
        return existing
    provenance: Dict[tuple, Dict[str, list[str]]] = {}
    setattr(interpolator, _PROVENANCE_ATTR, provenance)
    return provenance


def set_interpolator_provenance(interpolator: Any, provenance: Dict[tuple, Dict[str, list[str]]]):
    setattr(interpolator, _PROVENANCE_ATTR, provenance)


def copy_interpolator_provenance(source: Any, target: Any):
    set_interpolator_provenance(target, deepcopy(get_interpolator_provenance(source)))


def seed_original_sample_provenance(interpolator: Any, original_sample_positions: Optional[Iterable[tuple]] = None):
    allowed = None if original_sample_positions is None else set(original_sample_positions)
    provenance = get_interpolator_provenance(interpolator)
    for pos, block in getattr(interpolator, 'blocks', {}).items():
        if not _get_block_is_sample(block):
            continue
        if allowed is not None and pos not in allowed:
            continue
        entry = provenance.setdefault(pos, {'source': [], 'algorithm': []})
        _append_unique(entry['source'], 'Original Sample')
        _append_unique(entry['algorithm'], 'Sample')


def snapshot_interpolator_state(interpolator: Any) -> Dict[tuple, tuple]:
    target = _current_target(interpolator)
    snapshot: Dict[tuple, tuple] = {}
    for pos, block in getattr(interpolator, 'blocks', {}).items():
        if target == 'domain':
            snapshot[pos] = ('domain', str(_get_block_attr(block, 'domain', '') or '').strip())
        else:
            snapshot[pos] = ('value', _normalize_scalar(_get_block_attr(block, 'value')))
    return snapshot


def finalize_phase_provenance(interpolator: Any, source_label: str, algorithm_name: Any, before_snapshot: Optional[Dict[tuple, tuple]] = None) -> int:
    before = before_snapshot or {}
    after = snapshot_interpolator_state(interpolator)
    algorithm_label, _ = normalize_provenance_algorithm_name(algorithm_name)
    provenance = get_interpolator_provenance(interpolator)
    changed = 0
    missing = object()
    for pos, state in after.items():
        previous = before.get(pos, missing)
        if previous is missing or previous != state:
            entry = provenance.setdefault(pos, {'source': [], 'algorithm': []})
            _append_unique(entry['source'], source_label)
            _append_unique(entry['algorithm'], algorithm_label)
            changed += 1
    return changed


def get_export_provenance(interpolator: Any, pos: tuple) -> Optional[Tuple[str, str, str]]:
    provenance = get_interpolator_provenance(interpolator)
    entry = provenance.get(pos)
    if not entry:
        return None
    source_parts = [str(part).strip() for part in entry.get('source', []) if str(part).strip()]
    algorithm_parts = [str(part).strip() for part in entry.get('algorithm', []) if str(part).strip()]
    if not source_parts or not algorithm_parts:
        return None
    _, algo_type = normalize_provenance_algorithm_name(algorithm_parts[-1])
    return ' + '.join(source_parts), ' + '.join(algorithm_parts), algo_type