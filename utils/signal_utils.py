#  signal_utils.py - this file is part of the gwadaptive_scattering package.
#  Copyright (C) 2020- Stefano Bianchi, Alessandro Longo
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see <https://www.gnu.org/licenses/>.


import numpy as np
from gwpy.timeseries import TimeSeries
from gwpy.timeseries import TimeSeriesDict
from scipy.signal import kaiserord, lfilter, firwin, butter, freqz
from scipy.signal import hilbert
from scipy.stats import pearsonr
import pytvfemd
from ..common import defines

LAMBDA = 1.064


def smooth(arr, win):
    """Smooth array.
    
    Parameters
    ----------
    arr : numpy array
        input array
    win : int
        smoothing window
        
    Returns
    -------
    numpy array
        smoothed array
    """
    if win % 2 == 0:
        win += 1
    out = np.convolve(arr, np.ones(win, dtype=int), "valid") / win
    r = np.arange(1, win-1, 2)
    start = np.cumsum(arr[:win-1])[::2] / r
    stop = (np.cumsum(arr[:-win:-1])[::2] / r)[::-1]

    return np.concatenate((start, out, stop))


def lowpass(x, cutoff):
    """Lowpass filter.
    
    Parameters
    ----------
    x : numpy array
        input array
    cutoff : float
        cutoff frequency
        
    Returns
    -------
    numpy array
        lowpassed array
    """
    x = TimeSeries(x - np.nanmean(x))
    x_low = x.lowpass(1 / cutoff).value

    return x_low
    # passband_ripple = 60
    # passband_ripple_db = (10 ** (passband_ripple / 20) - 1) / (10 ** (passband_ripple / 20) + 1)
    # cutoff_norm = cutoff / (f_samp / 2)
    # width_percentage = -0.98 * steepness + 0.99
    # width = width_percentage * (1 - cutoff_norm)
    # num_taps, beta = kaiserord(passband_ripple, width)
    # taps = firwin(num_taps, cutoff_norm, window=('kaiser', beta))
    # filtered_x = lfilter(taps, 1.0, x)

    # return filtered_x


def butter_lowpass(cutoff, f_samp, order=3):
    """Butter filter.
    
    Parameters
    ----------
    cutoff : float
        cutoff frequency
    f_samp : float
        sampling frequency
    order : int, optional
        filter order (default : 3)
        
    Returns
    -------
    tuple of float
        filter coefficients
    """
    nyq = 0.5 * f_samp
    normal_cutoff = cutoff / nyq
    response = butter(order, normal_cutoff, btype="lowpass", output="ba", analog=False)

    return response[0], response[1]


def butter_lowpass_filter(x, cutoff, f_samp, order=3):
    """Butter filter.
    
    Parameters
    ----------
    x : numpy array
        input array
    cutoff : float
        cutoff frequency
    f_samp : float
        sampling frequency
    order : int, optional
        filter order (default : 3)
        
    Returns
    -------
    numpy array
        filtered array
    """
    b, a = butter_lowpass(cutoff, f_samp, order=order)
    y = lfilter(b, a, x)

    return y


def mean_frequency_wa(x, f):
    """Array mean frequency.
    
    Parameters
    ----------
    x : numpy array
        input array
    f : float
        sampling frequency
        
    Returns
    -------
    float
        mean frequency
    """
    spec = np.abs(np.fft.rfft(x))
    freq = np.fft.rfftfreq(len(x), d=1 / f)
    amp = spec / spec.sum()
    mean_f = (freq * amp).sum()

    return mean_f


def mean_frequency(channel_name, start, end, bandpass_limits=None, verbose=False):
    """Mean frequency of a signal (highest peak in the spectrum).

    Parameters
    ----------
    channel_name : str
        channel name
    start : int
        gps start
    end : int
        gps end
    bandpass_limits : tuple[float]
        frequencies interval for bandpass filter
    verbose : bool
        verbosity

    Returns
    -------
    float
        mean frequency of the signal
    """
    ts = TimeSeries.get(channel_name, start, end, verbose=verbose)
    if bandpass_limits is not None:
        bp_start = bandpass_limits[0]
        bp_end = bandpass_limits[1]
        ts = ts.bandpass(bp_start, bp_end)
    asd = ts.asd()
    mf_idx = asd == asd.max()
    mf = asd.frequencies[mf_idx]

    return mf.value[0]


def mean_amplitude(x, f):
    """Array mean amplitude.
    
    Parameters
    ----------
    x : numpy array
        input array
    f : float
        sampling frequency
        
    Returns
    -------
    float
        mean amplitude
    """
    dft = (2 / len(x)) * np.fft.rfft(x)
    freq = np.fft.rfftfreq(len(x), d=1 / f)
    mean_f_idx = np.argmin(np.abs(freq - mean_frequency(x, f)))

    return np.abs(dft[mean_f_idx])


def get_predictor(time, ts, N=1, smooth_win=None):
    """Time series predictor.
    
    Parameters
    ----------
    time : numpy array
        time array
    ts : numpy array
        time series
    N : int, optional
        scattering factor (default : 1)
    smooth_win : int, optional
        smoothing window (default : None)
        
    Returns
    -------
    numpy array
        predictor
    """
    v_mat = np.diff(ts) / np.diff(time)
    if smooth_win is not None:
        v_mat = smooth(v_mat, smooth_win)
    pred = N * (2 / LAMBDA) * np.abs(v_mat)
    
    return pred


def get_predictors(channels, fs, smooth_win=None, n_scattering=1):
    """Get predictors from channels.

    Parameters
    ----------
    channels : numpy ndarray
        channels matrix (shape (len(channels_values), num_channels)
    fs : float
        channels common sampling frequency
    smooth_win : int, optional
        smoothing window (default : None)
    n_scattering : int, optional
        scattering factor (default : 1)

    Returns
    -------
    numpy ndarray
        channels predictors
    """
    time = np.arange(0, len(channels) / fs, 1 / fs, dtype=float)
    v_mat = np.diff(channels, axis=0) / np.diff(time)[:, np.newaxis]
    if smooth_win is not None:
        for i in range(len(v_mat[0, :])):
            v_mat[:, i] = smooth(v_mat[:, i], smooth_win)

    predictors = n_scattering * (2 / LAMBDA) * np.abs(v_mat)
    fs_int = int(fs)
    predictors = predictors[(defines.EXTRA_SECONDS * fs_int):-(defines.EXTRA_SECONDS * fs_int), :]
    # for i in range(predictors.shape[1]):
    #     predictors[:, i] = butter_lowpass_filter(predictors[:, i], f_lowpass, fs)

    return predictors


def get_imfs(target_channel, fs, norm=True):
    """Get imfs with pytvfemd.

    Parameters
    ----------
    target_channel : numpy array
        channel from which extract imfs
    fs : float
        `target_channel` sampling frequency
    norm : bool, optional
        normalize imfs (default : True)

    Returns
    -------
    numpy ndarray
        `target_channel` imfs matrix
    """
    imfs = pytvfemd.tvfemd(target_channel)
    fs_int = int(fs)
    imfs = imfs[(defines.EXTRA_SECONDS * fs_int):-(defines.EXTRA_SECONDS * fs_int), :]
    if norm:
        imfs = (imfs - np.nanmean(imfs, axis=0)) / np.nanstd(imfs, axis=0)

    return imfs


def get_data_from_time_series_dict(target_channel_name, channels_list, gps_start, gps_end,
                                   fs=None, verbose=False, frametype=None):
    """Get data from `gwpy` function `TimeSeriesDict.get`.

    Parameters
    ----------
    target_channel_name : str
        name of the target channel
    channels_list : list[str]
        list of auxiliary channels names
    gps_start : int
        starting GPS
    gps_end : int
        ending GPS
    fs : float, optional
        desired sampling frequency for the channels (default : None)
    verbose : bool, optional
        verbosity option of TimeSeriesDict (default : False)
    frametype : str, optional
        frametype for desired data (default : None)

    Returns
    -------
    numpy ndarray
        matrix with channels values, first column
        corresponds to the target channel
    float
        common sampling frequency of the channels
        in the matrix
    """
    if fs is None:
        fs = np.inf

    if frametype is None:
        data_dict = TimeSeriesDict.get([target_channel_name] + channels_list,
                                       gps_start - defines.EXTRA_SECONDS,
                                       gps_end + defines.EXTRA_SECONDS,
                                       verbose=verbose)
    else:
        data_dict = TimeSeriesDict.get([target_channel_name] + channels_list,
                                       gps_start - defines.EXTRA_SECONDS,
                                       gps_end + defines.EXTRA_SECONDS,
                                       verbose=verbose, frametype=frametype)
    dict_fs = np.min([data_dict[ch_name].channel.sample_rate.value for ch_name in channels_list])
    if dict_fs < fs:
        fs = dict_fs
    data_dict.resample(fs)
    data = np.zeros((data_dict[target_channel_name].value.shape[0], len(channels_list) + 1), dtype=float)
    data[:, 0] = data_dict[target_channel_name].value
    for i in range(1, len(channels_list) + 1):
        data[:, i] = data_dict[channels_list[i - 1]].value

    return data, fs


def upper_envelope(ts):
    """Get upper envelope of a signal.

    Parameters
    ----------
    ts : numpy ndarray
        input signal

    Returns
    -------
    numpy ndarray
        input signal's envelope
    """
    analytic_signal = hilbert(ts)
    upper_env = np.abs(analytic_signal)

    return upper_env


def get_correlation_between(x, y):
    """Get Pearson correlation between `x` and `y`.

    Parameters
    ----------
    x : numpy ndarray
        first series
    y : numpy ndarray
        second series

    Returns
    -------
    float
        correlation value
    """
    r_corr, _ = pearsonr(x, y)
    r_corr = -999.0 if np.isnan(r_corr) else r_corr

    return r_corr
