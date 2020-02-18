#!/usr/bin/env python3
import numpy as np
import time
from scipy import signal

def xcor(apple, banana):
    '''1D Cross-Correlation'''
    corr = signal.correlate(apple, banana, mode='same', method='fft')
    return np.abs(corr)

def apply_fdoa(ray, fdoa, samp_rate):
    precache = 2j * np.pi * fdoa / samp_rate
    new_ray = np.empty_like(ray)
    for idx, val in enumerate(ray):
        new_ray[idx] = val * np.exp(precache * idx)
    return new_ray

def amb_surf(needle, haystack, freqs_hz, samp_rate):
    '''
    Returns the cross ambiguity function surface for a pair of signals.

    Parameters
    ----------
    needle : np.ndarray
        The signal of interest to localize within the haystack.
    haystack : np.ndarray
        The broader capture within which to localize the needle.
    freqs_hz : np.ndarray
        The frequency offsets to use in computing the CAF.
    samp_rate : float
        The sample rate for both the needle and the haystack.

    Returns
    -------
    surf : np.ndarray
        2D array of correlations of the needle in the haystack over frequency x lag.
    '''
    len_needle = len(needle)
    len_haystack = len(haystack)
    len_freqs = len(freqs_hz)
    assert len_needle == len_haystack
    surf = np.empty((len_freqs, len_needle))
    for fdx, freq_hz in enumerate(freqs_hz):
        shifted = apply_fdoa(needle, freq_hz, samp_rate)
        surf[fdx] = xcor(shifted, haystack)
    return surf

if __name__ == '__main__':
    # TODO: Ambiguity plot is broken at the moment
    import matplotlib.pyplot as plt
    import os
    # from mpl_toolkits.mplot3d import Axes3D

    start_time = time.time()
    data_dir = '../data'
    needle_filename = 'chirp_4_raw.c64'
    haystack_filename = 'chirp_4_T+70samp_F+82.89Hz.c64'
    print(haystack_filename)
    needle_samples = np.fromfile(os.path.join(data_dir, needle_filename), dtype=np.complex64)
    haystack_samples = np.fromfile(os.path.join(data_dir, haystack_filename), dtype=np.complex64)[0:4096]
    len_needle = len(needle_samples)

    samp_rate = 48e3
    freq_offsets = np.arange(-100, 100, 0.5)

    start = time.time()
    surf = amb_surf(needle_samples, haystack_samples, freq_offsets, samp_rate)
    print('elap',time.time()-start)
    print(surf.shape, surf.dtype)


    fmax, tmax = np.unravel_index(surf.argmax(), surf.shape)

    tau_max = len(needle_samples)//2 - tmax
    freq_max = freq_offsets[fmax]
    print(tau_max, freq_max)

    stop_time = time.time()

    run_time = stop_time - start_time

    print('CAF compute time: {} sec'.format(run_time))

    extents = [
        -len_needle//2, len_needle//2,
        100, -100]

    plt.figure(dpi=150)
    plt.imshow(surf, aspect='auto', interpolation='nearest', extent=extents)
    plt.ylabel('Frequency offset [Hz]')
    plt.xlabel('Time offset [samples]')
    plt.gca().invert_yaxis()
    plt.plot(tau_max, freq_max, 'x', color='red', alpha=0.75)
    plt.show()

    print('Time lag: {:d} samples'.format(tau_max))
    print('Frequency offset: {:.2f} Hz'.format(freq_max))

    # fig = plt.figure(dpi=150)
    # ax = fig.add_subplot(111, projection='3d')
    # x = tau
    # y = freq_offsets
    # X, Y = np.meshgrid(x, y)
    # Z = surf.reshape(X.shape)
    #
    # ax.plot_surface(X, Y, Z, cmap='viridis')
    #
    # ax.set_xlabel('Frequency offset [Hz]')
    # ax.set_ylabel('Lag [samples]')
    # ax.set_zlabel('Correlation')
    # plt.show()
