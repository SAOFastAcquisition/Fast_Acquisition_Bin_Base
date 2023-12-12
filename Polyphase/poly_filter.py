import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

from scipy.fftpack import fft, fftshift
from scipy.signal.windows import kaiser, flattop

# Matplotlib default params
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'white'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
# plt.rcParams['legend.loc'] = 'best'


def freqtime_compare(freq=201, suptitle='', l=1):
    n = 512
    m = n * l

    # Add AWGN
    np.random.seed(42)
    awgn = np.random.normal(0, 9 * 1e-2, m)

    # Input _signal
    t = np.linspace(0, 1, m)

    x = np.cos(2 * np.pi * freq * t) + awgn
    x1 = x.reshape(l, n)
    x2 = x.reshape(l, n).sum(0)

    # FFT
    Xx = np.abs(fft(x[:n], n))[:n // 2]
    X1 = (np.abs(fft(x1, n))).sum(axis=0)[:n // 2]
    X2 = np.abs(fft(x2, n))[:n // 2]

    Xxlog = 20 * np.log10(Xx / Xx.max())
    X1log = 20 * np.log10(X1 / X1.max())
    X2log = 20 * np.log10(X2 / X2.max())

    fn = np.linspace(0, 0.5, n // 2)

    # Plot
    plt.figure(figsize=(14, 2), dpi=120)
    plt.suptitle(suptitle)
    plt.subplot(1, 3, 1)
    plt.plot(fn, Xxlog, color='C0', label=f'N={n}')
    plt.ylim([-60, 0])
    plt.legend(loc='upper right')

    plt.subplot(1, 3, 2)
    plt.plot(fn, X1log, color='C2', label='Avg (freq)')
    plt.ylim([-60, 0])
    plt.legend(loc='upper right')

    plt.subplot(1, 3, 3)
    plt.plot(fn, X2log, color='C3', label='Avg (time)')
    plt.ylim([-60, 0])
    plt.legend(loc='upper right')
    # plt.show()


def comb_filter(x_in, d):
    sample_length = np.size(x_in)
    y = np.zeros(sample_length)
    for n in range(sample_length):
        if n == 0:
            y[n] = x_in[n]
        elif n < d:
            y[n] = x_in[n] + y[n - 1]
        else:
            y[n] = x_in[n] - x_in[n - d] + y[n - 1]

    return y / d


def interpolate_sequence(x_in, r):
    input_length = np.size(x_in)
    y = np.zeros(input_length * r)
    for n in range(input_length):
        y[n * r] = x_in[n]
    return y


def decimate_sequence(x_in, d):
    input_length = np.size(x_in)
    y = np.zeros(input_length // d)
    for n in range(int(input_length // d)):
        y[n] = x_in[n * d]
    return y


def random_signal(n):
    primary_sig = np.random.sample(n) - 0.5
    return primary_sig


m = 2 ** 16
signal_sin = [0.1 * (np.sin(2 * np.pi * 0.3 * i) + 0.75 * np.sin(2 * np.pi * 0.17 * i)) for i in range(m)]
signal_rand = random_signal(m) + signal_sin
sig_mean = np.mean(signal_rand)
sig_var = np.var(signal_rand)
signal_rand = comb_filter(signal_rand, 7)
# r = 4
# signal_rand = interpolate_sequence(signal_rand, r)
# signal_rand = comb_filter(signal_rand, 7)
d = 4
signal_rand = decimate_sequence(signal_rand, d)
signal_rand_sample = np.reshape(signal_rand, (1024, -1))    # Разбиваем на реализации длиной 1024 отсчетов
spectrum_signal_rand = fft(signal_rand_sample, 1024)
spectrum_signal_av = np.average(np.abs(spectrum_signal_rand ** 2), 0)

axes = plt.subplots()
plt.plot(spectrum_signal_av)
plt.show()

freqtime_compare(freq=121 * 1, l=1, suptitle='NumAverage = 1')

l_avg = 15
freqtime_compare(freq=121*l_avg, l=l_avg, suptitle=f'NumAverage = {l_avg}')
pass

N, L = 16, 5
M = N*L

t = np.linspace(0, 1, M)
f1, f2 = L+1, 4*L

x = np.cos(2*np.pi*f1*t) + np.cos(2*np.pi*f2*t)
y1 = x.reshape(L, N)
y = y1.sum(0)

# FFT
XM = np.abs(fft(x, M))[:M//2]
XN = np.abs(fft(y, N))[:N//2]
XM /= XM.max()
XN /= XN.max()

# Interpolated spectrum
XZ = np.zeros(M//2)
#XZ[::L] = XM[::L]
XZ[::L] = XN

plt.figure(figsize=(12, 5), dpi=120)
plt.subplot(2, 2, 1)
plt.plot(x, '-o', markersize=6, color='C0', label='Input _signal')
plt.legend(loc='upper right')
plt.subplot(2, 2, 2)
plt.plot(y, '-*', markersize=6, color='C1', label='Decimated _signal')
plt.legend(loc='upper right')
plt.subplot(2, 2, 3)
plt.stem(XM, use_line_collection=True, linefmt='C2', basefmt='C2', label='Spectrum (_signal)')
plt.legend(loc='upper right')
plt.subplot(2, 2, 4)
plt.stem(XZ, use_line_collection=True, linefmt='C3', basefmt='C3', label='Spectrum (decimated)')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()