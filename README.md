# Anterpolator

**3D geological interpolation and block-model operations toolkit**

Anterpolator is a 3D geological modeling tool focused on interpolation, domain-aware workflows, and block-model operations. It combines experimental bio-inspired interpolators with practical data-preparation utilities for samples and blocks, making it useful both for generating interpolated models and for auditing, exporting, and transforming geological datasets.

## Key Features

### 1. Ant Colony Interpolator
Inspired by the foraging behavior of ants, this algorithm is excellent for tracing continuous, curvilinear structures like:
- High-grade ore shoots
- Vein networks
- Structurally controlled mineralization

Ants release "pheromones" based on sample grades, establishing connected paths through the model space that respect the natural continuity of the deposit.

### 2. String Theory Interpolator
A connectivity-driven interpolator that builds "strings" between compatible samples.
- Supports configurable distance threshold, grade-difference control, collision policy, and processing order.
- Can interpolate either numeric values or categorical domains.
- Includes optional azimuth and dip frequency filtering to keep dominant structural trends.

### 3. Molecular Clock Interpolator (Phylogeographic Approach)
A groundbreaking approach that treats spatial dispersion as evolutionary divergence.
- **Common Ancestor Inference:** Finds "feeder" zones or sources in depth.
- **Multi-Event Detection:** Automatically separates distinct geological events (pulses) using DBSCAN clustering.
- **Tree-Based Connectivity:** Builds minimum spanning trees to connect samples, ideal for magmatic intrusions, dikes, and hydrothermal systems.

### 4. Adaptive Octree Interpolator
A bottom-up domain-aware hierarchical interpolator that aggregates support over the configured finest block grid.
- Supports **Dense Blocks Cover** for regular final grids and **Adaptive Leaf Cover** for non-overlapping mixed-size output blocks.
- Preserves domain boundaries while propagating values upward only within each domain.
- Includes configurable **Support Density Alpha** weighting, using $s / (V^\alpha)$ to penalize large sparsely supported regions when desired.

### 5. Gaussian Kernel Interpolator
A non-iterative interpolation mode that fills blocks from nearby sample influence using configurable kernel bandwidth, cutoff distance, and nearest-sample fallback behavior.

### 6. Two-Pass and Domain-Specific Workflows
Interpolation is not limited to a single algorithm run.
- Configure a first-pass and second-pass algorithm in one workflow.
- Use the output of pass 1 as the input of pass 2.
- Override algorithms per domain so different geological domains can use different interpolation strategies.

### 7. Operations Panel for Samples, Blocks, and Exports
Beyond interpolation, Anterpolator includes a dedicated operations workflow for common geological data tasks:
- **Sample Blocks:** aggregate samples into their containing blocks.
- **Assign Block Columns To Samples:** transfer selected block attributes onto samples. Selecting the domain column here replaces the former dedicated Domain Samples action, and sample-to-block matching now honors sub-block extents when DX/DY/DZ columns are available.
- **Assign Block Columns To Block Model:** enrich a target block model from overlapping source blocks while preserving target rows; numeric fields use overlap-volume-weighted averages, categorical fields use the dominant overlap category, and a nearest-block fallback can be enabled for targets with no overlap.
- **Assign Attributes From Table:** join one or more CSV table columns onto a block model by matching one or more shared key columns, preserving all block-model rows while updating or appending the selected attributes.
- **Block Domain Sample Metrics:** export nearest-distance, k-nearest, residual, and summary metrics by domain.
- **Matching Value Sample Distances:** map multiple block value/category fields to sample fields and calculate the closest matching sample both within each block's domain and across all domains.
- **Domain Interpolation Confidence:** summarize spacing and support metrics for each domain.
- **Block Volume Weighted Average:** compute weighted summaries using metadata-defined or inferred block volumes, or a chosen weight field.
- **CSV Grid To BMF:** export regular CSV grids to the experimental TBMS2.0 BMF backend.
- **Equation Finder By Domain:** run symbolic regression per domain to search for candidate equations.

Block-to-sample and block-to-block transfers use explicit DX/DY/DZ columns first, then
Leapfrog **regular** subdivision metadata (`parent_block_size / subblock_factors`),
and centroid inference only when neither is available. Volume-based summaries also
use regular subdivision metadata before inference. For example, a 15x15x15 parent
with regular 3x3x3 subdivision has 5x5x5 source rows even when rows are missing or
filtered out. Keep the GUI base block size at 15x15x15 in this case.
Octree metadata defines a minimum size, not necessarily each row's actual size,
so it does not override centroid inference for mixed-size models.

#### Matching Value Sample Distances

Choose the samples and blocks files and their XYZ columns under **Files & Data**.
Select an explicit domain field for each file under **Samples Columns** and
**Blocks Columns**; the sample and block field names may differ.
In **Operations > Matching Value Sample Distances**, use **Add Field** to map each
block value/category field to its corresponding sample field, choose an output CSV,
and click **Export Matching Value Distances**.

For a block field `RQDpc` and block domain field `UG`, the operation appends
`RQDpc_UG_MinDist` (same value and same domain) and `RQDpc_AbsMinDist` (same value,
any domain). The domain portion of the output name uses the **block** domain
field. Distances are exact 3D Euclidean distances from block centroids to sample
coordinates, in the input coordinate units; both files must use the same coordinate
system. Matching trims whitespace, remains case-sensitive for categories, and
normalizes equivalent numeric tokens such as `20` and `20.0`; no value tolerance
or nearest nonmatching fallback is used.

App-wide sample and block filters apply. All filtered block rows and original
columns are retained, with one pair of distance columns per mapped block field.
Invalid/nonfinite coordinates, missing values, and absent matching samples produce
blank distances and are reported in the completion summary. Missing domains exclude
rows only from same-domain matching: samples with missing domains or domains not
present in the block model still participate in absolute matching.
The default output is `<blocks>_matching_value_distances.csv`; input files cannot
be overwritten. Field mappings and the output path are saved with the configuration.
Loading reports byte-based progress for each file, followed by filtering and
preparation updates. Only the required sample columns (including filter fields)
are read; all block columns are retained. Repeated values and domains are
normalized once per distinct value rather than once per row.
CSV headers may contain quoted separators or quotes. Unnamed block columns are
retained as `Unnamed` (with numbered suffixes for additional unnamed columns).

### 8. Viewer and Export Utilities
- Interactive 3D viewer for blocks, samples, and interpolation output.
- CSV export of interpolation results.
- Optional sub-block-aware export expansion.
- Standalone BMF export tooling for grid conversion workflows.

## Project Structure

- **Anterpolator3DViewer/**: Core Python application and viewer, including the main GUI, interpolation engines, and BMF export backend.
- **utils/**: Utility scripts (e.g., backup management).
- **DOCS/**: Documentation.

## Getting Started

### Prerequisites
- Python 3.8+
- Recommended: Create a virtual environment.

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/rhaguirrem/Anterpolator.git
   cd Anterpolator
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

### Usage

Run the main viewer application:

```bash
python Anterpolator3DViewer/anterpolator3DViewer.py
```

## Configuration

The application uses `config.json` files to manage interpolation, workflow, and export parameters. You can customize:
- **Ant Colony:** Pheromone levels, ant count, sampling percentage, range size, visited-threshold behavior, and interpolation target.
- **Molecular Clock:** Spatial vs. attribute weights, event detection sensitivity, ancestor depth bias, and interpolation method.
- **Gaussian Kernel:** Bandwidth, cutoff sigma, nearest fallback, and background handling.
- **String Theory:** Distance threshold, grade difference, connection limits, collision policy, processing order, interpolation target, and structural frequency filtering.
- **Adaptive Octree:** Output mode, maximum aggregation levels, support-density alpha weighting, and dense-output provenance columns.
- **Workflow:** Default algorithm, second-pass algorithm, and domain-specific algorithm overrides.
- **Operations:** Filters, selected columns, metrics, export paths, and BMF export settings.

## Author

**Roberto Aguirre Maturana, Geologist**

---
*This project explores geology-focused interpolation, connectivity analysis, and block-model operations through a mix of experimental and production-oriented tooling.*
