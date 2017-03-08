import numpy as np
from scipy.optimize import leastsq

def fit_sine(t,data,guesses):
    """
    Guesses: list
      The initial guess values for amplitude, angular freq, phase, and shift
    """
    N = len(data)

    ## define a sine function
    optimize_func = lambda x: x[0]*np.sin(x[1] * t + x[2]) + x[3] - data
    est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, guesses)[0]

    ## recreate the fitted curve using the optimized parameters
    data_fit = est_amp*np.sin(est_freq * t + est_phase) + est_mean

    return data_fit, [est_amp,est_freq,est_phase,est_mean]
