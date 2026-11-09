import numpy as np
from scipy.spatial import cKDTree

from ant_colony import AntColonyInterpolator, Block
from provenance_utils import get_interpolator_provenance


def fill_isolated_blank_blocks(blocks, config):
    """Fill remaining blank base blocks using fixed same-domain inverse-square IDW donors."""
    if not config.get('fill_isolated_blanks', False):
        return 0
    neighbors = config.get('isolated_fill_neighbors', 8)
    if isinstance(neighbors, bool) or not isinstance(neighbors, int) or neighbors < 1:
        raise ValueError('Isolated blank fill requires a positive integer neighbor count.')

    info = blocks._block_info
    source_domains = info.get('domain_mapping') or {}
    allowed = info.get('allowed_grid')
    if allowed is None:
        dims = info['dims']
        allowed = [(i, j, k) for i in range(dims[0]) for j in range(dims[1]) for k in range(dims[2])]
    groups = (
        [(domain, passes[-1]) for domain, passes in blocks._interpolators.items() if passes]
        if hasattr(blocks, '_interpolators')
        else [(None, blocks._ant_colony)]
    )
    scale = np.asarray(info['block_size'], dtype=float)
    filled = 0
    for group_domain, interp in groups:
        if str(getattr(interp, 'interpolation_target', 'value')).lower() == 'domain':
            continue
        if getattr(interp, 'output_mode', '') == 'adaptive_leaf_cover':
            raise ValueError('Isolated blank fill requires a base-grid output, not Adaptive Leaf Cover.')
        provenance = get_interpolator_provenance(interp)
        donors = {}
        targets = {}
        for pos in allowed:
            pos = tuple(pos)
            if group_domain is not None and group_domain != 'Undomained':
                if source_domains.get(pos, 'Undomained') != group_domain:
                    continue
            block = interp.blocks.get(pos)
            domain = str(source_domains.get(pos, group_domain or 'default'))
            if block is None:
                targets.setdefault(domain, []).append(pos)
                continue
            value = block.get('value') if isinstance(block, dict) else block.value
            if value is not None and np.isfinite(value):
                donors.setdefault(domain, []).append((pos, float(value)))
            elif not (block.get('is_sample', False) if isinstance(block, dict) else block.is_sample):
                targets.setdefault(domain, []).append(pos)

        for domain, positions in targets.items():
            domain_donors = donors.get(domain, [])
            if not domain_donors:
                print(f"Isolated blank fill: leaving {len(positions):,} blocks blank in {domain!r}; no valued donors.")
                continue
            donor_points = np.asarray([pos for pos, _ in domain_donors], dtype=float) * scale
            donor_values = np.asarray([value for _, value in domain_donors])
            tree = cKDTree(donor_points)
            k = min(neighbors, len(domain_donors))
            batch_size = max(1, min(50_000, 1_000_000 // k))
            for start in range(0, len(positions), batch_size):
                batch = positions[start:start + batch_size]
                distances, indices = tree.query(np.asarray(batch, dtype=float) * scale, k=k)
                distances = np.asarray(distances).reshape(len(batch), k)
                indices = np.asarray(indices).reshape(len(batch), k)
                # Normalize before squaring to avoid large inverse-square weights.
                weights = (distances[:, :1] / distances) ** 2
                values = np.sum(weights * donor_values[indices], axis=1) / weights.sum(axis=1)
                for pos, value, index in zip(batch, values, indices[:, 0]):
                    block = interp.blocks.get(pos)
                    if block is None:
                        if isinstance(interp, AntColonyInterpolator):
                            block = Block(
                                value=float(value), is_sample=False,
                                mark_class=interp.get_mark_class(float(value)),
                                block_id=interp.next_block_id, domain=domain,
                                nearest_sample_value=float(donor_values[index]),
                            )
                            interp.next_block_id += 1
                            interp.domain_created_counts[domain] = interp.domain_created_counts.get(domain, 0) + 1
                        else:
                            block = {'value': float(value), 'is_sample': False, 'domain': domain}
                        interp.blocks[pos] = block
                    elif isinstance(block, dict):
                        block['value'] = float(value)
                        block['domain'] = domain
                    else:
                        block.value = float(value)
                        block.domain = domain
                        block.mark_class = interp.get_mark_class(float(value))
                    entry = provenance.setdefault(pos, {'source': [], 'algorithm': []})
                    if 'Post-process' not in entry['source']:
                        entry['source'].append('Post-process')
                    if 'Same-domain IDW Fill' not in entry['algorithm']:
                        entry['algorithm'].append('Same-domain IDW Fill')
                    filled += 1
    print(f"Isolated blank fill: assigned {filled:,} blocks using up to {neighbors} same-domain neighbors.")
    return filled
