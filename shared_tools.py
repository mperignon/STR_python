import numpy as np

def central_difference(array, dx, alpha = 0.5, right_boundary = None, left_boundary = None):

    return_array = np.zeros_like(array)

    if right_boundary is None:
        return_array[0] = (array[0] - array[1]) / dx
    else:
        return_array[0] = alpha * (right_boundary - array[0]) / dx + \
                (1 - alpha) * (array[0] - array[1]) / dx

    if left_boundary is None:
        return_array[-1] = (array[-2] - array[-1]) / dx
    else:
        return_array[-1] = alpha * (array[-2] - array[-1]) / dx + \
                (1 - alpha) * (array[-1] - left_boundary) / dx

    for i in range(1,len(array)-1):
        return_array[i] = alpha * (array[i-1] - array[i]) / dx + \
                (1 - alpha) * (array[i] - array[i+1]) / dx

    return return_array
        