#!/usr/bin/env python3
'''build pulses for testing'''
import numpy as np
# from gnuradio.filter import firdes
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
    # rrc_taps = firdes.root_raised_cosine(.95, samp_rate, samp_rate*3/8, 0.35, 11*4)
    rrc_taps = np.array([ 0.00118141,  0.00057044, -0.00139341, -0.00073436,  0.001506  ,
        0.00056174, -0.00208351, -0.00046321,  0.00326889,  0.00113577,
       -0.00413576, -0.00176473,  0.00473822,  0.00040722, -0.00905426,
        0.00101337,  0.02324399,  0.00914119, -0.04807992, -0.05711537,
        0.07358722,  0.28460228,  0.38973358,  0.28460228,  0.07358722,
       -0.05711537, -0.04807992,  0.00914119,  0.02324399,  0.00101337,
       -0.00905426,  0.00040722,  0.00473822, -0.00176473, -0.00413576,
        0.00113577,  0.00326889, -0.00046321, -0.00208351,  0.00056174,
        0.001506  , -0.00073436, -0.00139341,  0.00057044,  0.00118141])
    pulse_filtered = sig.lfilter(rrc_taps, [1], pulse)
    pulse_filtered = pulse
    edge = np.zeros((next_size(size) - size) // 2)
    return np.hstack((edge, pulse_filtered, edge))

def next_size(size):
    pow = np.ceil(np.log2(size)).astype(np.int)
    return 2**pow

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
        banana = apple + np.random.random(len(apple))*2-1
        banana = apply_fdoa(banana, fdoa, samp_rate)
        banana = apply_tdoa(banana, tdoa, samp_rate)
        filename = os.path.join(data_dir,'burst_{:04.0f}_raw.f32'.format(idx))
        apple.astype(np.float32).tofile(filename)
        filename = os.path.join(data_dir,'burst_{:04.0f}_t{:+.6f}_f{:+.6f}.f32'.format(idx,tdoa*samp_rate,fdoa*samp_rate))
        banana.astype(np.float32).tofile(filename)
        print('{}: tdoa, fdoa = ({:+.6f} samples, {:+.6f} Hz)'.format(idx,tdoa*samp_rate, fdoa*samp_rate))
