from math import sqrt, log
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
            
        self._datapts = [(4, 100), (2, 99), (1, 97), (0.5, 83.4), ...
                         (0.25, 42), (0.125, 10), (0.062, 3.2), (0.031, 2)]
        
        self._D_g = None
        self._sigma_g = None



        
        
    @property
    def datapoints(self):
        """Grain bound diameters (mm) and percent finer"""
        return self._datapts


    @datapoints.setter
    def datapoints(self, new_dp):
        """Set the datapoints.

        Parameters
        ----------
        new_dp: list of tuples of bound diameters (in mm) and percent finer
        """
        
        self._datapoints = new_dp
        self.pre_process_data()
        
        

    @property
    def geometric_mean(self):
        """Geometric mean (mm) of grain size distribution"""
        return self._D_g
        
    @property
    def std(self):
        """Standard deviation (mm) of grain size distribution"""
        return self._sigma_g           
        
        
        
         
    def run(self):
    
        assert self._datapts is not None, "Grain size data must be set before running."
        
        self.grain_statistics()
        self.finalize()
        
        
            
    def finalize(self):
    
        sorting = 'Well sorted'
        if self._sigma_g >= 1.6:
            sorting = 'Poorly sorted'
    
        print 'The statistics of the grain size distribution are:'
        print 'Geometric mean:', self._D_g, 'mm'
        print 'Standard deviation:', self._sigma_g, 'mm    [', sorting, ']'
        
        
    
    
    def pre_process_data(self):
    
        dtype = [('diameters', float), ('fractions', float)]
        self._datapts = np.array(self._datapoints, dtype=dtype)

        # scale
        if max(self._datapts['fractions']) > 1.:
            self._datapts['fractions'] = self._datapts['fractions'] / 100.

        # sort
        self._datapts = np.sort(self._datapts, order='diameters')       

        # interpolate
        if self._datapts['fractions'][0] > 0.0:

            # if smallest fraction is not zero, interpolate linearly
            # to find the missing value
    
            D3 = log(self._datapts['diameters'][1],2)
            D2 = log(self._datapts['diameters'][0],2)
            f3 = self._datapts['fractions'][1]
            f2 = self._datapts['fractions'][0]
            newpsi = D3 + ((D2 - D3)/(f2 - f3))*(0-f3)
            first_entry = np.array((2**newpsi,0), dtype=dtype)
            self._datapts = np.hstack((first_entry, self._datapts))
    
        if self._datapts['fractions'][-1] < 1.0:

            # if largest fraction is not zero, assume that the next larger
            # psi scale grain size had a fraction of 1
    
            psi_scale = np.array(range(-10,20), dtype = np.float)
            psi_scale_mm = 2**psi_scale
    
            next_larger_size = \
                psi_scale_mm[ np.where(psi_scale_mm > \
                    self._datapts['diameters'][-1])[0][0]]
                
            last_entry = np.array([(next_larger_size, 1.0)], dtype = dtype)
            self._datapts = np.hstack((self._datapts, last_entry))    
    

    def characteristic_size(self, pct):
        '''
        Calculate characteristic size based on percent finer
    
        Input
        ------
        datapts_: nd array of grain size diameters and percent finers
        pct: percentage such that D_pct is the size such that pct percent
            of the sample is finer than D_pct (between 0 and 100)
        
        Output
        ------
        D_x (in millimiters), where x percent of the sample is finer than Dx
        
        
        pct can be array-like, in which case D_x is the same length
        '''
        
        if type(pct) is int or type(pct) is float:
            pct = [pct]
            
        assert len(pct), 'pct is the wrong format!'
        
        self.D_x = []
            
        for i in pct:

            pct_ = float(i)
    
            loc_ = np.where(self._datapts['fractions'] >= pct_/100)[0][0]
            assert loc_>0, 'Percent finer appears to be less than 0!'

            D2_ = log(self._datapts['diameters'][loc_],2)
            D1_ = log(self._datapts['diameters'][loc_-1],2)
            f2_ = self._datapts['fractions'][loc_]
            f1_ = self._datapts['fractions'][loc_-1]

            psix_ = D1_ + ((D2_ - D1_)/(f2_ - f1_)) * (pct_/100 - f1_)
            self.D_x.append(2**psix_)
            
        
        


    def grain_statistics(self):
        '''
        Statistical characteristics of grain size distributions
    
        Input
        -----
        datapts: nd array of grain size diameters and percent finers
    
        Output
        ------
        D_g: Geometric mean grain diameter in millimiter
        sigma_g: Standard deviation of grain size in millimiters
            Sediment is well sorted if sigma_g < 1.6
        '''

        psi_b = [log(i,2) for i in self._datapts['diameters']]

        dtype = [('psi', float), ('diameters', float), ('fractions', float)]
        self._centerpts = np.zeros((len(self._datapts)-1), dtype=dtype)

        for i in range(len(self._centerpts)):    
            self._centerpts['psi'][i] = (psi_b[i+1] + psi_b[i]) / 2
            self._centerpts['diameters'][i] = sqrt(self._datapts['diameters'][i] * \
                    self._datapts['diameters'][i+1])
            self._centerpts['fractions'][i] = self._datapts['fractions'][i+1] - \
                    self._datapts['fractions'][i]

        self._psi_bar = sum(self._centerpts['psi'] * self._centerpts['fractions'])
        self._D_g = 2 ** self._psi_bar

        self._sigma = sqrt(sum((self._centerpts['psi'] - self._psi_bar)**2 * \
            self._centerpts['fractions']))
        self._sigma_g = 2 ** self._sigma      
            
            
            
            
            
            
            
            
            
        





    
            
            
            
            
            
            
            