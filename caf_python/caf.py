import numpy as np

def caf(needle, haystack, frequency_offsets, sample_rate=1):
    '''
    Returns the cross ambiguity function surface for a given signal of interest to be localized in a capture.
    Parameters

    Parameters
    ----------
    needle : np.ndarray
        The signal of interest to localize within the haystack.
    haystack : np.ndarray
        The broader capture within which to localize the needle.
    frequency_offsets : np.ndarray
        The frequency offsets to use in computing the CAF.
    sample_rate : float
        The sample rate for both the needle and the haystack.

    Returns
    -------
    caf : np.ndarray
        2D array of correlations of the needle in the haystack over frequency x lag.
    '''
    num_needle_samples = len(needle)
    num_haystack_samples = len(haystack)
    num_frequency_offsets = len(frequency_offsets)
    
    freq_shift_vectors = np.exp(1j*2*np.pi*np.outer(frequency_offsets, np.arange(num_needle_samples) / sample_rate))
    shifted_needle = np.multiply(freq_shift_vectors, needle)

    shift_vectors = np.exp(1j*2*np.pi*np.outer(frequency_offsets, np.arange(num_needle_samples) / sample_rate))
    translated_soi = np.multiply(shift_vectors, soi_samples)
    del shift_vectors
    caf_surface = np.zeros((num_frequency_offsets, num_needle_samples + num_haystack_samples - 1), dtype=np.float64)
    
    for i in range(num_freq_offsets):
        caf_surface[i, :] = np.abs(np.correlate(search_capture_samples, (translated_soi[i, :]), mode='full'))
        
    return caf_surface
