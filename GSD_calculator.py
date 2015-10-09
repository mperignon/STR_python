from math import sqrt, log10, pow
import numpy as np

class GSD_calculator(object):
    """
    GSD_calculator

    Python version of Gary Parker's 1D Sediment Transport Morphodynamics e-book,
    originally in Visual Basic and converted to C by Andrew Leman (2009)

    MC Perignon
    Oct 2015

    ----------------------------------------------

    Computes, from a grain size distribution, the geometric mean diameter, geometric
    standard deviation and characteristic diameters based on percent finer.

    ----------------------------------------------
    Input:
    datapoints: list of tuples of bound diameters (in mm) and percent finer
    [(diameter, fraction), (diameter, fraction), ...]
    
    If there are no lower and upper bounds in the inputted distribution such that
    the percent finer are equal to 0 and 1 (or 0 and 100), compute these bounds
    with a linear interpolation in the phi scale.

    
    """
    
    
    def __init__(self, verbose = True):
        
        self._verbose = verbose
        if verbose:
            print 'Initializing RouseVanoni_Equilibrium...'
            
        self._datapoints = None



        
        
    @property
    def datapoints(self):
        """Grain bound diameters (mm) and percent finer"""
        return self._datapoints


    @datapoints.setter
    def datapoints(self, new_dp):
        """Set the datapoints.

        Parameters
        ----------
        new_dp: list of tuples of bound diameters (in mm) and percent finer
        """
        self._datapoints = new_dp

            
            
    def run(self):
    
        assert self._datapoints, "Grain size data must be set before running."
        
        # convert the list of tuples of datapoints to separate arrays
        
        
        # re-scale the percent finers so they are 0 to 1
            
        
        
            
    def finalize(self):
    
        pass   
                
            
            
    def tuples_to_arrays(self):
    
        self._diameters = []
        self._fractions = []
            
        for i in self._datapoints:
        
            self._diameters.append(i[0])
            self._fractions.append(i[1])
            
        self._diameters = np.array(self._diameters)
        self._fractions = np.array(self._fractions)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            