import numpy as np

po = np.array((0.6684, 0.7639, 0.8601, 0.9096, 0.9615, 1.0000, 1.055, 1.108, 1.197, 1.302, 1.407, 1.529, 1.641, 1.702, 1.832, 1.937, 2.044, 2.261, 2.499, 2.732, 2.993, 3.477, 4.075, 4.469, 5.016, 6.158, 7.821, 10.06, 14.38, 19.97, 25.79, 38.57, 68.74, 91.95, 231.2, 2320.0))

oo = np.array((1.011, 1.011, 1.01, 1.008, 1.004, 0.9997, 0.9903, 0.9789, 0.9567, 0.9273, 0.8964, 0.8604, 0.8287, 0.8123, 0.7796, 0.7554, 0.7326, 0.6928, 0.6585, 0.6345, 0.615, 0.5877, 0.564, 0.5523, 0.5395, 0.5209, 0.5045, 0.4917, 0.479, 0.4712, 0.4668, 0.462, 0.4578, 0.4564, 0.4541, 0.4527))

so = np.array((0.8157, 0.8157, 0.8182, 0.8233, 0.8333, 0.8439, 0.8621, 0.8825, 0.9214, 0.9723, 1.025, 1.083, 1.13, 1.153, 1.196, 1.225, 1.25, 1.287, 1.313, 1.333, 1.352, 1.38, 1.403, 1.414, 1.426, 1.444, 1.458, 1.469, 1.48, 1.486, 1.49, 1.493, 1.497, 1.498, 1.499, 1.5))

def Chezy_Cf(Cz):
    
    Cf = 1. / Cz**2
    
    return Cf

def CubicInterpolant(ym1, y, yp1, yp2, xm1, x, xp1, xp2, phisgo):

    a1 = ym1;
    a2 = (y - a1)/(x - xm1);
    a3 = (yp1 - a1 - a2*(xp1 - xm1))/((xp1 - xm1)*(xp1 - x));
    a4 = yp2 - a1 - a2*(xp2 - xm1) - a3*(xp2 - xm1)*(xp2 - x);
    a4 = a4/((xp2 - xm1)*(xp2-x)*(xp2-xp1));

    output = a1 + a2*(phisgo - xm1) + a3*(phisgo-xm1)*(phisgo - x) + a4*(phisgo - xm1)*(phisgo - x)*(phisgo - xp1);

    return output


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
    
