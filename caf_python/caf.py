#!/usr/bin/env python3
import numpy as np
import time

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

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import os
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    start_time = time.time()
    data_dir = '../data'
    search_capture_filename = 'chirp_4_T+70samp_F+82.89Hz.c64'
    soi_filename = 'chirp_4_raw.c64'
    soi_samples = np.fromfile(os.path.join(data_dir, soi_filename), dtype=np.complex64)
    search_capture_samples = np.fromfile(os.path.join(data_dir, search_capture_filename), dtype=np.complex64)

    sample_rate = 48e3

    freq_offset_min = -120
    freq_offset_max = 120
    num_freq_offsets = 400

    freq_offsets = np.linspace(freq_offset_min, freq_offset_max, num_freq_offsets)


    caf_surface = caf(needle=soi_samples,
                      haystack=search_capture_samples,
                      frequency_offsets=freq_offsets,
                      sample_rate=sample_rate)


    peak = np.unravel_index(np.argmax(caf_surface), caf_surface.shape)

    tau_max = peak[1] - len(soi_samples) + 1
    tau = np.arange(-len(soi_samples) + 1, len(search_capture_samples))
    freq_max = freq_offsets[peak[0]]

    stop_time = time.time()

    run_time = stop_time - start_time

    print('CAF compute time: {} sec'.format(run_time))

    plt.figure(dpi=150)
    plt.imshow(np.abs(caf_surface), aspect='auto', interpolation='nearest', extent=[-len(soi_samples)+1, len(search_capture_samples), freq_offsets[-1], freq_offsets[0]])
    plt.ylabel('Frequency offset [Hz]')
    plt.xlabel('Time offset [samples]')
    plt.gca().invert_yaxis()
    plt.plot(tau_max, freq_max, 'x', color='red', alpha=0.75)
    plt.show()

    print('Time lag: {:d} samples'.format(tau_max))
    print('Frequency offset: {:.2f} Hz'.format(freq_max))

    fig = plt.figure(dpi=150)
    ax = fig.add_subplot(111, projection='3d')
    x = tau
    y = freq_offsets
    X, Y = np.meshgrid(x, y)
    Z = caf_surface.reshape(X.shape)

    ax.plot_surface(X, Y, Z, cmap=cm.inferno)

    ax.set_xlabel('Frequency offset [Hz]')
    ax.set_ylabel('Lag [samples]')
    ax.set_zlabel('Correlation')
    plt.show()
