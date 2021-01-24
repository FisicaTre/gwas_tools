import numpy as np
from gwpy.timeseries import TimeSeries
from scipy.signal import kaiserord, lfilter, firwin, butter, freqz


def smooth(arr, win):
    if win % 2 == 0:
        win += 1
    out = np.convolve(arr, np.ones(win, dtype=int), "valid") / win
    r = np.arange(1, win-1, 2)
    start = np.cumsum(arr[:win-1])[::2] / r
    stop = (np.cumsum(arr[:-win:-1])[::2] / r)[::-1]

    return np.concatenate((start, out, stop))


def lowpass(x, cutoff):
    x = TimeSeries(x - np.nanmean(x))
    x_low = x.lowpass(1 / cutoff).value

    return x_low
    #passband_ripple = 60
    #passband_ripple_db = (10 ** (passband_ripple / 20) - 1) / (10 ** (passband_ripple / 20) + 1)
    #cutoff_norm = cutoff / (f_samp / 2)
    #width_percentage = -0.98 * steepness + 0.99
    #width = width_percentage * (1 - cutoff_norm)
    #num_taps, beta = kaiserord(passband_ripple, width)
    #taps = firwin(num_taps, cutoff_norm, window=('kaiser', beta))
    #filtered_x = lfilter(taps, 1.0, x)

    #return filtered_x
    
def butter_lowpass(cutoff, f_samp, order=3):
    nyq = 0.5 * f_samp
    normal_cutoff = cutoff / nyq
    response = butter(order, normal_cutoff, btype="lowpass", output="ba", analog=False)

    return response[0], response[1]


def butter_lowpass_filter(x, cutoff, f_samp, order=3):
    b, a = butter_lowpass(cutoff, f_samp, order=order)
    y = lfilter(b, a, x)

    return y


def mean_frequency(x, f):
    spec = np.abs(np.fft.rfft(x))
    freq = np.fft.rfftfreq(len(x), d=1 / f)
    amp = spec / spec.sum()
    mean_f = (freq * amp).sum()

    return mean_f


def mean_amplitude(x, f):
    dft = (2 / len(x)) * np.fft.rfft(x)
    freq = np.fft.rfftfreq(len(x), d=1 / f)
    mean_f_idx = np.argmin(np.abs(freq - mean_frequency(x, f)))

    return np.abs(dft[mean_f_idx])


def predictor(time, ts, N=1, LAMBDA=1.064, smooth_win=None):
    v_mat = np.diff(ts) / np.diff(time)
    if smooth_win is not None:
        v_mat = smooth(v_mat, smooth_win)
    pred = N * (2 / LAMBDA) * np.abs(v_mat)
    
    return pred
