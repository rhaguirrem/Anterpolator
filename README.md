# Anterpolator

**Numeric Interpolator based on Ant Colony Optimization and Molecular Clock Phylogeography**

Anterpolator is a novel 3D spatial interpolation tool designed for geological modeling. It moves beyond traditional geostatistical methods (kriging, IDW) and implicit modeling (RBF) by applying bio-inspired algorithms to solve complex connectivity problems in geological structures.

## Key Features

### 1. Ant Colony Interpolator
Inspired by the foraging behavior of ants, this algorithm is excellent for tracing continuous, curvilinear structures like:
- High-grade ore shoots
- Vein networks
- Structurally controlled mineralization

Ants release "pheromones" based on sample grades, establishing connected paths through the model space that respect the natural continuity of the deposit.

### 2. Molecular Clock Interpolator (Phylogeographic Approach)
A groundbreaking approach that treats spatial dispersion as evolutionary divergence.
- **Common Ancestor Inference:** Finds "feeder" zones or sources in depth.
- **Multi-Event Detection:** Automatically separates distinct geological events (pulses) using DBSCAN clustering.
- **Tree-Based Connectivity:** Builds minimum spanning trees to connect samples, ideal for magmatic intrusions, dikes, and hydrothermal systems.

## Project Structure

- **Anterpolator3DViewer/**: Core Python application and viewer.
  - `anterpolator3DViewer.py`: Main entry point and GUI.
  - `ant_colony.py`: Implementation of the Ant Colony Optimization algorithm.
  - `molecular_clock_interpolator.py`: Implementation of the Phylogeographic/Molecular Clock algorithm.
  - `interpolator_base.py`: Abstract base class for interpolators.
- **Ant Hill/**: Legacy/Alternative implementations and data.
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

The application uses `config.json` files to manage parameters for both algorithms. You can customize:
- **Ant Colony:** Pheromone levels, ant count, range size.
- **Molecular Clock:** Spatial vs. attribute weights, event detection sensitivity, ancestor depth bias.

## Author

**Ramon Aguirre M., PhD**

---
*This project explores the intersection of biology and geostatistics to solve complex geological modeling challenges.*
