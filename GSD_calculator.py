from math import sqrt, log
import numpy as np
import json


def pre_process_data(datapoints):
    '''
    Processes a grain size distribution (diameters and mm and fractions) for 
    calculating bedload grain size distributions
    '''

    dtype = [('diameters', float), ('fractions', float)]
    datapts = np.array(datapoints, dtype=dtype)

    # scale
    if max(datapts['fractions']) > 1.:
        datapts['fractions'] = datapts['fractions'] / 100.

    # sort
    datapts = np.sort(datapts, order='diameters')       

    # interpolate
    if datapts['fractions'][0] > 0.0:

        # if smallest fraction is not zero, interpolate linearly
        # to find the missing value

        D3 = log(datapts['diameters'][1],2)
        D2 = log(datapts['diameters'][0],2)
        f3 = datapts['fractions'][1]
        f2 = datapts['fractions'][0]
        newpsi = D3 + ((D2 - D3)/(f2 - f3))*(0-f3)
        first_entry = np.array((2**newpsi,0), dtype=dtype)
        datapts = np.hstack((first_entry, datapts))

    if datapts['fractions'][-1] < 1.0:

        # if largest fraction is not zero, assume that the next larger
        # psi scale grain size had a fraction of 1

        psi_scale = np.array(range(-10,20), dtype = np.float)
        psi_scale_mm = 2**psi_scale

        next_larger_size = \
            psi_scale_mm[ np.where(psi_scale_mm > \
                datapts['diameters'][-1])[0][0]]

        last_entry = np.array([(next_larger_size, 1.0)], dtype = dtype)
        datapts = np.hstack((datapts, last_entry))

    return datapts 


def grain_statistics(datapts):
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

    psi_b = [log(i,2) for i in datapts['diameters']]

    dtype = [('psi', float), ('diameters', float), ('fractions', float)]
    centerpts = np.zeros((len(datapts)-1), dtype=dtype)

    for i in range(len(centerpts)):    
        centerpts['psi'][i] = (psi_b[i+1] + psi_b[i]) / 2
        centerpts['diameters'][i] = sqrt(datapts['diameters'][i] * \
                datapts['diameters'][i+1])
        centerpts['fractions'][i] = datapts['fractions'][i+1] - \
                datapts['fractions'][i]

    psi_bar = sum(centerpts['psi'] * centerpts['fractions'])
    D_g = 2 ** psi_bar

    sigma = sqrt(sum((centerpts['psi'] - psi_bar)**2 * \
        centerpts['fractions']))
    sigma_g = 2 ** sigma

    return D_g, sigma_g




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
    ------
    datapoints: list of tuples of bound diameters (in mm) and percent finer
    [(diameter, fraction), (diameter, fraction), ...]
    
    If there are no lower and upper bounds in the inputted distribution such that
    the percent finer are equal to 0 and 1 (or 0 and 100), compute these bounds
    with a linear interpolation in the phi scale.
    
    Output:
    -------
    json file of grain size distribution statistics

    """
    
    
    def __init__(self, datapts = None, verbose = False, save_output = True):
        
        self._verbose = verbose
        if verbose:
            print 'Initializing RouseVanoni_Equilibrium...'
            
        self._datapoints = datapts
        
        if not datapts:
            self._datapoints = [(4, 100), (2, 99), (1, 97), (0.5, 83.4), 
                             (0.25, 42), (0.125, 10), (0.062, 3.2), (0.031, 2)]
            
        self.save_output = save_output

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
        
        

    @property
    def geometric_mean(self):
        """Geometric mean (mm) of grain size distribution"""
        return self._D_g
        
    @property
    def std(self):
        """Standard deviation (mm) of grain size distribution"""
        return self._sigma_g           
        
        
        
         
    def run(self):
        '''
        Performs calculations
        '''
        
        self._datapts = pre_process_data(self._datapoints)
        self._D_g, self._sigma_g = grain_statistics(self._datapts)
        self.finalize()
        
        
            
    def finalize(self):
        '''
        Creates output files
        '''
    
        sorting = 'Well sorted'
        if self._sigma_g >= 1.6:
            sorting = 'Poorly sorted'
            
        pct_passing = self.characteristic_size([70,50,30])
    
        if self._verbose:
            print 'The statistics of the grain size distribution are:'
            print 'Geometric mean:', self._D_g, 'mm'
            print 'Standard deviation:', self._sigma_g, 'mm    [', sorting, ']'
            
        
        if self.save_output:
            
            output_dict = {
                'Geometric_mean' : self._D_g,
                'Standard_deviation' : self._sigma_g,
                'Sorting' : sorting,
                'D_70' : pct_passing[0],
                'D_50' : pct_passing[1],
                'D_30' : pct_passing[2],
            }
            
            with open('output/GSD_calculator_stats.json', 'w') as f:
                json.dump(output_dict, f, indent=4, sort_keys=True)
    
    

    def characteristic_size(self, pct = None):
        '''
        Calculate characteristic size based on percent finer
    
        Input
        ------
        pct: percentage such that D_pct is the size such that pct percent
            of the sample is finer than D_pct (between 0 and 100)
        
        Output
        ------
        D_x (in millimiters), where x percent of the sample is finer than Dx
        
        
        pct can be array-like, in which case D_x is the same length
        '''
        
        if not pct:
            pct = [70,50,30]
        
        if type(pct) is int or type(pct) is float:
            pct = [pct]
            
        assert len(pct), 'pct is the wrong format!'
        
        D_x = []
        
        datapts = self._datapts
            
        for i in pct:

            pct_ = float(i)
    
            loc_ = np.where(datapts['fractions'] >= pct_/100)[0][0]
            assert loc_>0, 'Percent finer appears to be less than 0!'

            D2_ = log(datapts['diameters'][loc_],2)
            D1_ = log(datapts['diameters'][loc_-1],2)
            f2_ = datapts['fractions'][loc_]
            f1_ = datapts['fractions'][loc_-1]

            psix_ = D1_ + ((D2_ - D1_)/(f2_ - f1_)) * (pct_/100 - f1_)
            D_x.append(2**psix_)
            
        return D_x
            
        
        



            
            
            
            
            
            
            
            
            
        





    
            
            
            
            
            
            
            