import numpy as np

def central_difference(array, dx, alpha = 0.5, right_boundary = None, left_boundary = None):

    return_array = np.zeros_like(array)

    if left_boundary is None:
        return_array[0] = (array[0] - array[1]) / dx
    else:
        return_array[0] = alpha * (left_boundary - array[0]) / dx + \
                (1 - alpha) * (array[0] - array[1]) / dx

    if right_boundary is None:
        return_array[-1] = (array[-2] - array[-1]) / dx
    else:
        return_array[-1] = alpha * (array[-2] - array[-1]) / dx + \
                (1 - alpha) * (array[-1] - right_boundary) / dx

    for i in range(1,len(array)-1):
        return_array[i] = alpha * (array[i-1] - array[i]) / dx + \
                (1 - alpha) * (array[i] - array[i+1]) / dx

    return return_array



def create_surface_arrays(nodes, dx, S, depth = 0, direction = "downstream"):
    """
    Creates arrays for distance, bed elevation, water surface elevation, and depth.
    
    If direction = "downstream", highest elevations are on the left (low indices)
    If direction = "upstream", integration happens upstream and thus the lowest elevations
    are on the left
    """
    
    L = nodes*dx
    x = np.array(range(nodes)) * dx
    eta = S*L - S*(x+dx)
    eta[-1] = 0
    H = np.zeros_like(x) + depth
    ksi = eta + H
    
    if direction == "upstream":
        eta = np.flipud(eta)
        
    ksi = eta + H
    
    return x, eta, H, ksi
    
