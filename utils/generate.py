#!/usr/bin/env python3
'''Make some chirpy things for testing'''
import numpy as np
import scipy.signal as sig
import os
import scipy.interpolate as interp
import scipy as sp
import scipy.signal

def apply_offset(signal, dfc, sample_rate):
    '''
    Applies a constant or time-varying frequency offset to the signal.
    Returns a copy of the signal with the frequency offset applied.
    '''
    if type(dfc) in (int, float):
        shift = np.exp(1j*2*np.pi*dfc*np.arange(len(signal))/sample_rate)
    else:
        phi = np.cumsum(2*np.pi*dfc) / sample_rate
        shift = np.exp(1j*(np.arange(len(signal))/sample_rate+phi))
    return shift*signal

def generate_chirp(sample_rate, chirp_length=4096, chirp_order=2, relative_bandwidth=1e-2, sweep_range_Hz=10e3, window=np.hanning):
    # Generate some shaped noise
    kernel = sp.signal.firwin(127, cutoff=0.5*relative_bandwidth, fs=sample_rate)
    srange = np.random.uniform(1e3, 10e3)
    chirp = np.random.normal(0, 1, chirp_length) + 1j*np.random.normal(0, 1, chirp_length)
    chirp = sp.signal.filtfilt(kernel, 1, chirp)

    # Taper the edges
    if window is not None:
        chirp = window(chirp_length)*chirp
    chirp = chirp.astype(np.complex64)

    # Make it move
    shape = np.linspace(-1, 1, chirp_length)**chirp_order
    offset_Hz = shape*sweep_range_Hz
    chirp = apply_offset(chirp, offset_Hz, sample_rate)

    return chirp

if __name__ == '__main__':
    np.random.seed(0)
    data_dir = '../data'
    samp_rate = 48e3

    chirp_length = 4096
    chirp_order = np.random.randint(2,5) # Determines shape of chirp
    relative_bandwidth = np.random.uniform(1e-3, 5e-2) # Width of chirp, relative to sample rate
    sweep_range_Hz = np.random.uniform(1e3, 10e3) # Range of chirp

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    for idx in range(10):
        dfc_range_Hz = 1e2 # Range of frequency offsets for search capture
        lag = np.random.randint(7, 256) # Lag (in samples) of SOI in search capture
        chirp = generate_chirp(chirp_length=chirp_length, chirp_order=chirp_order, relative_bandwidth=relative_bandwidth, sweep_range_Hz=sweep_range_Hz, sample_rate=samp_rate)
        chirp = chirp.astype(np.complex64)
        chirp.tofile(os.path.join(data_dir, 'chirp_{:d}_raw.c64'.format(idx)))

        # Add a random time lag
        foffset = np.random.uniform(-dfc_range_Hz, dfc_range_Hz)
        chirp_search = np.concatenate([np.zeros(lag), chirp, np.zeros(96)])
        chirp_search = apply_offset(chirp_search, foffset, samp_rate)

        # Add some noise
        chirp_search += np.random.normal(0, 1e-5, len(chirp_search)) + 1j*np.random.normal(0, 1e-5, len(chirp_search))
        chirp_search = chirp_search.astype(np.complex64)
        chirp_search.tofile(os.path.join(data_dir, 'chirp_{:d}_T{:+d}samp_F.{:+.2f}Hz.c64'.format(idx, lag, foffset)))
