import numpy as np

def get_1_sigma_range(x):
    '''computes the 1-sigma range from the data set provided
    
    Parameters
    ----------
    x : array
        Data set for which you want the 1-sigma range
        
    Returns
    -------
    1-sigma range : float
        1-sigma range for the data set
        NOTE: you must divide by 2 to get the 1-sigma value
    ''' 
    x_sorted = np.sort(x)
    x_median = x_sorted[int(0.5*len(x))]
    x_low = x_sorted[int(0.16*len(x))]
    x_high = x_sorted[int(0.84*len(x))]

    return x_high - x_low

def get_m1(m2, C):
    '''computes the bright component mass based on the 
    dark component mass and mass function parameter:
    C = m_f/(sin i)^3
    
    Parameters
    ----------
    m2 : float or array
        Dark component mass(es)
    
    C : float or array
        Mass function parameter(s)
    Returns
    -------
    m1 : float or array
        Bright component mass(es) in the units of m2, C
    ''' 
    m1 = m2*(np.sqrt(m2/C)-1)
    return m1

def dm1_relative(m2, C, dm2, dC):
    '''computes the relative measuremnt error on the 
    bright component mass based on the dark component mass 
    and mass function parameter: C = m_f/(sin i)^3, 
    and their observed measurement errors. This assumes
    gaussian measurement error!
    
    Parameters
    ----------
    m2 : float or array
        Dark component mass(es)
    
    C : float or array
        Mass function parameter(s)
        
    dm2 : float or array
        Measurement error for dark component mass(es)
    
    dC : float or array
        Measurement error for mass function parameter(s)
    Returns
    -------
    dm1_by_m1 : float or array
        Relative measurement error for bright component 
        mass(es) 
    ''' 
    dm1_by_m1 = ( ((1/m2 + 0.5/(M2-C*np.sqrt(M2/C)))*dm2)**2 +
                  ((-0.5*M2/(m2*C - C**2*np.sqrt(M2/C)))*dC)**2 )**0.5
    return dm1_by_m1

def get_m2(m1, C):
    '''computes the dark component mass based on the 
    bright component mass and mass function parameter:
    C = m_f/(sin i)^3
    
    Parameters
    ----------
    m1 : float or array
        Bright component mass(es)
    
    C : float or array
        Mass function parameter(s)
    Returns
    -------
    m2 : float or array
        Dark component mass(es) in the units of m1, C
    ''' 
    xi = 2*C**3 + 18*C**2*m1 + 3*(3**0.5)*(4*C**3*m1**3 + 27*C**2*m1**4)**0.5 + 27*C*m1**2
    
    m2 = C/3.0 + xi**(1./3)/(3*2**(1./3)) + 2**(1./3)*(C**2 + 6*C*m1)/(3*xi**(1./3))
    return m2

def dm2(m1, C, dm1, dC):
    '''computes the measuremnt error on the dark component 
    mass based on the bright component mass 
    and mass function parameter: C = m_f/(sin i)^3, 
    and their observed measurement errors. This assumes
    gaussian measurement error!
    
    Parameters
    ----------
    m1 : float or array
        Dark component mass(es)
    
    C : float or array
        Mass function parameter(s)
        
    dm1 : float or array
        Measurement error for dark component mass(es)
    
    dC : float or array
        Measurement error for mass function parameter(s)
    Returns
    -------
    dm2 : float or array
        Measurement error for dark component mass(es) 
    '''
    xi = 2*C**3 + 18*C**2*m1 + 3*(3**0.5)*(4*C**3*m1**3 + 27*C**2*m1**4)**0.5 + 27*C*m1**2
    alpha = 3*3**0.5*(12*C**3*m1**2 + 108*C**2*m1**3)
    beta = 2*(4*C**3*m1**3 + 27*C**2*m1**4)**0.5
    gamma = 3*(3**0.5)*(12*C**2*m1**3 + 54*C*m1**4)
    zeta = 18*C**2 + alpha/beta + 54*C*m1
    delta = 6*C**2 + gamma/beta + 36*C*m1 + 27*m1**2
    
    term_1 = ( (zeta/(9*2**(1./3)*xi**(2./3))) + 2*2**(1./3)*C/xi**(1./3) - (2**(1./3)*(C**2 + 6*C*m1)*zeta/(9*xi**(4./3))) )**2*(dm1**2)
    term_2 = ( 1./3 + (delta/(9*2**(1./3)*xi**(2./3))) + 2**(1./3)*(2*C + 6*m1)/(3*xi**(1./3)) - (2**(1./3)*(C**2 + 6*C*m1)*delta/(9*xi**(4./3))) )**2*(dC**2)
    
    dm2 = (term_1 + term_2)**(0.5)
    return dm2


C = 5.0
C_err = C*0.001

M1 = 1.0
M1_err = M1*0.001

M1_set = np.random.normal(M1, M1_err, 1000000)
C_set = np.random.normal(C, C_err, 1000000)
M2_obs = get_m2(M1_set, C_set)
M2_obs_err = get_1_sigma_range(M2_obs)/2
M2_out = np.median(M2_obs)

M2 = get_m2(M1, C) 
M2_err = dm2(M1, C, M1_err, C_err)

print('The median value of m2 from sampled m1 and C is {}'.format(M2_out))
print('The 1-sigma error of m2 from sampled m1 and C errors is {}'.format(M2_obs_err))

print('The true value of m2 from sampled m1 and C is {}'.format(M2))
print('The 1-sigma error of m2 is {}'.format(M2_obs_err))
