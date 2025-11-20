"""
Base abstract class for interpolation algorithms.
All interpolators (ant colony, biochemical clock, etc.) inherit from this.
"""
from abc import ABC, abstractmethod
from typing import Dict, Tuple, Any


class InterpolatorBase(ABC):
    """Abstract base class for all interpolation algorithms"""
    
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.blocks: Dict[Tuple[int, int, int], Any] = {}
        self.dims = None
        
    @abstractmethod
    def initialize_blocks(self, sample_blocks: Dict[Tuple[int, int, int], float], 
                         dims: Tuple[int, int, int], min_bounds, block_size, **kwargs):
        """
        Initialize the interpolation grid with sample data.
        
        Parameters:
        -----------
        sample_blocks : dict
            Dictionary mapping (x,y,z) coordinates to sample values
        dims : tuple
            Grid dimensions (nx, ny, nz)
        min_bounds : array-like
            Minimum bounds [xmin, ymin, zmin]
        block_size : array-like
            Block size [dx, dy, dz]
        **kwargs : additional algorithm-specific parameters
        """
        pass
    
    @abstractmethod
    def run_iteration(self, dims: Tuple[int, int, int]) -> bool:
        """
        Execute one iteration of the algorithm.
        
        Parameters:
        -----------
        dims : tuple
            Grid dimensions (nx, ny, nz)
            
        Returns:
        --------
        bool : True if algorithm should continue, False if converged/finished
        """
        pass
    
    @abstractmethod
    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        """
        Return interpolated values for all blocks.
        
        Returns:
        --------
        dict : Mapping of (x,y,z) coordinates to interpolated values
        """
        pass
    
    @abstractmethod
    def get_algorithm_name(self) -> str:
        """Return human-readable algorithm name"""
        pass
    
    def is_converged(self) -> bool:
        """
        Check if algorithm has converged (optional).
        Default implementation returns False (never converges).
        """
        return False
    
    def get_metadata(self) -> Dict[str, Any]:
        """
        Return algorithm-specific metadata (optional).
        Can include convergence info, tree structures, etc.
        """
        return {}
