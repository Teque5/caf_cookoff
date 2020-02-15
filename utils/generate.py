#!/usr/bin/env python3
'''build pulses for testing'''
import numpy as np
from gnuradio.filter import firdes
import scipy.signal as sig
import os
import scipy.interpolate as interp

def gen_burst(size=3000, samp_rate=48e3):
    '''modeled on the gnu radio implementation'''
    times = np.arange(size)/samp_rate
    glfsr_source = np.random.random(size=size) >= .5
    cosine_source = .2 + .8 * np.cos(2*np.pi*samp_rate/2**15*times)
    window = np.bartlett(size//2)
    triangle_source = np.hstack((window, window))
    pulse = glfsr_source * cosine_source * triangle_source
    rrc_taps = firdes.root_raised_cosine(.95, samp_rate, samp_rate*3/8, 0.35, 11*4)
    pulse_filtered = sig.lfilter(rrc_taps, [1], pulse)
    return pulse_filtered

def apply_fdoa(ray, fdoa, samp_rate):
    times = np.arange(len(ray))/samp_rate
    ray *= np.exp(-1j*2*np.pi*samp_rate*fdoa*times).real
    return ray

def apply_tdoa(ray, tdoa, samp_rate):
    times_old = np.arange(len(ray))/samp_rate
    times_new = times_old + tdoa
    kill = interp.interp1d(times_old, ray, fill_value="extrapolate")
    return kill(times_new)
    

if __name__ == '__main__':
    np.random.seed(0)
    data_dir = '../data'
    samp_rate = 48e3
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    for idx in range(10):
        fdoa = np.random.normal(scale=50)/samp_rate
        tdoa = np.random.normal(scale=50)/samp_rate
        amp = np.random.random()
        apple = gen_burst()
        banana = apple + np.random.random(3000)*2-1
        banana = apply_fdoa(banana, fdoa, samp_rate)
        banana = apply_tdoa(banana, tdoa, samp_rate) 
        filename = os.path.join(data_dir,'burst_{:04.0f}_raw.f32'.format(idx))
        apple.astype(np.float32).tofile(filename)
        filename = os.path.join(data_dir,'burst_{:04.0f}_t{:+.6f}_f{:+.6f}.f32'.format(idx,tdoa*samp_rate,fdoa*samp_rate))
        banana.astype(np.float32).tofile(filename)
        print('{}: tdoa, fdoa = ({:+.6f} samples, {:+.6f} Hz)'.format(idx,tdoa*samp_rate, fdoa*samp_rate))     
