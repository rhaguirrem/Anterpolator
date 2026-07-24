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
- **Domain Samples:** assign block-model domains back onto sample rows.
- **Assign Block Columns To Samples:** transfer selected block attributes onto samples.
- **Assign Block Columns To Block Model:** enrich a target block model from overlapping source blocks while preserving target rows; numeric fields use overlap-volume-weighted averages, categorical fields use the dominant overlap category, and a nearest-block fallback can be enabled for targets with no overlap.
- **Assign Attributes From Table:** join one or more CSV table columns onto a block model by matching one or more shared key columns, preserving all block-model rows while updating or appending the selected attributes.
- **Block Domain Sample Metrics:** export nearest-distance, k-nearest, residual, and summary metrics by domain.
- **Domain Interpolation Confidence:** summarize spacing and support metrics for each domain.
- **Block Volume Weighted Average:** compute weighted summaries using inferred block volumes or a chosen weight field.
- **CSV Grid To BMF:** export regular CSV grids to the experimental TBMS2.0 BMF backend.
- **Equation Finder By Domain:** run symbolic regression per domain to search for candidate equations.

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
