import numpy as np 
from smt.sampling_methods import LHS

def generate_LHS(dim, b_bounds, Lmax_bounds):
    """
    Generates a Latin Hypercube Sampling (LHS) grid of size dim*dim for the given bounds on b and Lmax.
    
    :int dim: The dimension of the grid.
    :List[float] b_bounds: The lower and upper bounds of b values.
    :List[float] Lmax_bounds: The lower and upper bounds of Lmax values.
    :return np.ndarray: The generated LHS grid.
    """
    grid_size = dim*dim
    xlimits = np.array(b_bounds, Lmax_bounds)
    lhs = LHS(xlimits=xlimits, criterion="maximin")
    grid = lhs(grid_size)
    
    return grid

def generate_uniform_grid(dim, b_bounds, Lmax_bounds):
    """
    Generates a uniform grid of size dim*dim for the given bounds on b and Lmax.
    
    :int dim: The dimension of the grid.
    :List[float] b_bounds: The lower and upper bounds of b values.
    :List[float] Lmax_bounds: The lower and upper bounds of Lmax values.
    :return np.ndarray: The generated grid.
    """
    grid = np.zeros((dim*dim, 2))
    b_grid = np.linspace(b_bounds[0], b_bounds[1], dim)
    Lmax_grid = np.linspace(Lmax_bounds[0], Lmax_bounds[1], dim)
    B, L = np.meshgrid(b_grid, Lmax_grid)

    grid[:, 0] = B.flatten()
    grid[:, 1] = L.flatten()
    
    return grid

if __name__ == "__main__":

    # Define your grid parameters
    dim = 50
    b_bounds = [0.0001, 0.075]
    Lmax_bounds = [0, 0.1]
    
    # Generate a Latin Hypercube Sampling (LHS) grid, optimal if using an interpolation method like kriging
    grid = generate_LHS(dim, b_bounds, Lmax_bounds)
    
    # Generate a uniform grid, optimal if not using any interpolation
    grid = generate_uniform_grid(dim, b_bounds, Lmax_bounds)

    # Save the grid in a .txt file
    np.savetxt("./LHS.txt", grid)