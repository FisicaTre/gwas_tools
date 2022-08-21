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
from gwpy.segments import DataQualityFlag
from scipy.signal import kaiserord, lfilter, firwin, butter, freqz
from scipy.signal import hilbert
from scipy.stats import pearsonr
import pytvfemd
import os
import re
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

    return predictors


def get_imfs(target_channel, fs, norm=True, max_imf=None):
    """Get imfs with pytvfemd.

    Parameters
    ----------
    target_channel : numpy array
        channel from which extract imfs
    fs : float
        `target_channel` sampling frequency
    norm : bool, optional
        normalize imfs (default : True)
    max_imf : int
        maximum number of imfs to be extracted (default : None)

    Returns
    -------
    numpy ndarray
        `target_channel` imfs matrix
    """
    is_max_imf_none = max_imf is None
    if is_max_imf_none:
        max_imf = 1000

    imfs = pytvfemd.tvfemd(target_channel, max_imf=max_imf + 1)
    fs_int = int(fs)
    imfs = imfs[(defines.EXTRA_SECONDS * fs_int):-(defines.EXTRA_SECONDS * fs_int), :]
    if norm:
        imfs = (imfs - np.nanmean(imfs, axis=0)) / np.nanstd(imfs, axis=0)

    n_imfs = imfs.shape[1]
    if not is_max_imf_none and n_imfs == max_imf + 1:
        imfs = imfs[:, :-1].reshape((imfs.shape[0], n_imfs - 1))

    return imfs


def get_combos(ch_str, envelopes, predictors):
    """Get correlation of combos for a given set of channels.

    Parameters
    ----------
    ch_str : list[str]
        channels correlated with each imf (the first channel corresponds to the first imf, and so on)
    envelopes : numpy array
        imf instantaneous amplitudes matrix
    predictors : numpy array
        predictors matrix

    Returns
    -------
    combos_imfs : list[list[int]]
        list of imfs lists belonging to the same channel
    combos_channels : list[str]
        list of culprits for each combo
    combos_corrs : list[float]
        list of correlations for each combo
    """
    combos_imfs = []
    combos_channels = []
    combos_corrs = []
    seen = set()
    uniq = [x for x in ch_str if x not in seen and not seen.add(x)]
    for u in uniq:
        idxs = np.where(np.array(ch_str) == u)[0]
        if len(idxs) != 1:
            sum_envelope = np.sum(envelopes[:, idxs], axis=1)
            c = get_correlation_between(sum_envelope, predictors[:, idxs[0]])
            combos_imfs.append([idx + 1 for idx in idxs])
            combos_channels.append(u)
            combos_corrs.append(c)

    return combos_imfs, combos_channels, combos_corrs


def get_ifo_of_channel(channel):
    """Get interferometer label for the input channel.

    Parameters
    ----------
    channel : str
        channel name

    Returns
    -------
    ifo : str
        interferometer label
    """
    return channel.split(":")[0]


def get_channel_wo_ifo(channel):
    """Get channel name without ifo.

    Parameters
    ----------
    channel : str
        channel name

    Returns
    -------
    name : str
        channel name without ifo
    """
    return channel.split(":")[1]


def get_instrument_lock_data(lock_channel, gps_start, gps_end):
    """Get data for instrument lock channel.

    Parameters
    ----------
    lock_channel : str
        name of the instrument lock channel
    gps_start : int
        starting GPS
    gps_end : int
        ending GPS

    Returns
    -------
    lock_data : gwpy.Segment, numpy ndarray
        for L1 : segments of instrument active periods in [`gps_start`, `gps_end`]
        for H1 : segments of instrument active periods in [`gps_start`, `gps_end`]
        for V1 : numpy array of lock channel values in [`gps_start`, `gps_end`]
    """
    lock_data = []
    ifo = get_ifo_of_channel(lock_channel)
    if ifo.startswith("L") or ifo.startswith("H"):
        lock_data = DataQualityFlag.query(lock_channel, gps_start, gps_end)
        lock_data = lock_data.active
    elif ifo.startswith("V"):
        lock_data = TimeSeriesDict.get([lock_channel], gps_start, gps_end)
        lock_data = lock_data[lock_channel].value

    return lock_data


def __get_gwf_files_of_interest(gwf_path, start_gps, end_gps, sep, start_gps_pos, n_gps_pos):
    gwf_files = []
    with open(gwf_path, "r") as gwf:
        for line in gwf:
            l = line.rstrip()
            gwf_files.append(re.split(r"\s+", l)[0])

    gwf_to_read = []
    if len(gwf_files) > 0:
        for gwf_file in gwf_files:
            flds = os.path.split(gwf_file)[1].split(".")[0].split(sep)
            if start_gps < int(flds[start_gps_pos]) + int(flds[n_gps_pos]) and end_gps >= int(flds[start_gps_pos]):
                gwf_to_read.append(gwf_file)

    return sorted(gwf_to_read)


def get_data_from_gwf_files(gwf_path, sep, start_gps_pos, n_gps_pos,
                            target_channel, channels, start_gps, end_gps,
                            samp_freq=None, **kwargs):
    """Get data from `gwf` files.

    Parameters
    ----------
    gwf_path : str
        path to `gwf` files
    sep : str
        separator character in a `gwf` file name
    start_gps_pos : int
        index (starting from 0) of the starting gps value in a `gwf` file name
    n_gps_pos : int
        index (starting from 0) of the number of seconds value in a `gwf` file name
    target_channel : str
        name of the target channel
    channels : list[str]
        list of auxiliary channels names
    start_gps : int
        starting GPS
    end_gps : int
        ending GPS
    samp_freq : float, optional
        desired sampling frequency for the channels (default : None)
    kwargs : dict{bool}
        gwpy.TimeSeriesDict keys

    Returns
    -------
    numpy ndarray
        matrix with channels values, first column
        corresponds to the target channel
    float
        common sampling frequency of the channels
        in the matrix
    """
    framelib = False
    try:
        from virgotools.everything import getChannel
        framelib = True
    except:
        pass

    if not framelib:
        if samp_freq is None:
            samp_freq = np.inf

    channels_list = [target_channel] + channels
    start_gps = start_gps - defines.EXTRA_SECONDS
    end_gps = end_gps + defines.EXTRA_SECONDS
    gwf_to_read = __get_gwf_files_of_interest(gwf_path, start_gps, end_gps, sep, start_gps_pos, n_gps_pos)

    data = {}  # dict[list[float]]
    for i, f in enumerate(gwf_to_read):
        flds = os.path.split(f)[1].split(".")[0].split(sep)
        s = max(start_gps, int(flds[start_gps_pos]))
        e = min(end_gps, int(flds[start_gps_pos]) + int(flds[n_gps_pos]))
        if not framelib:
            d = TimeSeriesDict.read(f, channels_list, start=s, end=e, **kwargs)

            if i == 0:
                dict_fs = np.min([d[ch_name].channel.sample_rate.value for ch_name in channels_list])
                if dict_fs < samp_freq:
                    samp_freq = dict_fs
            d.resample(samp_freq)
        else:
            d = {}
            for ch in channels_list:
                ts = getChannel(f, ch, s, e - s, mask=False, fill_value=np.nan)
                d[ch] = ts.data

        if i == 0:
            for k in d.keys():
                if not framelib:
                    data[k] = d[k].value
                else:
                    data[k] = d[k]
        else:
            for k in data.keys():
                if not framelib:
                    data[k] = np.concatenate((data[k], d[k].value), axis=None)
                else:
                    data[k] = np.concatenate((data[k], d[k]), axis=None)

    if framelib and samp_freq is not None:
        for k in data.keys():
            data[k] = data[k][::(len(data[k]) // (samp_freq * (end_gps - start_gps)))]

    data_mtx = np.zeros((data[target_channel].shape[0], len(channels_list)), dtype=float)
    data_mtx[:, 0] = data[target_channel]
    for i in range(1, len(channels_list)):
        data_mtx[:, i] = data[channels_list[i]]

    return data_mtx, samp_freq


def get_data_from_time_series_dict(target_channel_name, channels_list, gps_start, gps_end, fs, verbose=False):
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
    fs : float
        desired sampling frequency for the channels
    verbose : bool, optional
        verbosity option of TimeSeriesDict (default : False)

    Returns
    -------
    data : numpy ndarray
        matrix with channels values, first column
        corresponds to the target channel
    """
    from gwpy.io import datafind as io_datafind

    frametype_tc = io_datafind.find_frametype(target_channel_name, gpstime=gps_start, allow_tape=True)

    if frametype_tc is None:
        tc_dict = TimeSeriesDict.get([target_channel_name],
                                     gps_start - defines.EXTRA_SECONDS,
                                     gps_end + defines.EXTRA_SECONDS,
                                     verbose=verbose)
        tc_dict.resample(fs)
    else:
        tc_dict = TimeSeriesDict.get([target_channel_name],
                                     gps_start - defines.EXTRA_SECONDS,
                                     gps_end + defines.EXTRA_SECONDS,
                                     verbose=verbose, frametype=frametype_tc, resample=fs)

    data_dict = TimeSeriesDict.get(channels_list,
                                   gps_start - defines.EXTRA_SECONDS,
                                   gps_end + defines.EXTRA_SECONDS,
                                   verbose=verbose)
    data_dict.resample(fs)

    data = np.zeros((tc_dict[target_channel_name].value.shape[0], len(channels_list) + 1), dtype=float)
    data[:, 0] = tc_dict[target_channel_name].value
    for i in range(1, len(channels_list) + 1):
        data[:, i] = data_dict[channels_list[i - 1]].value

    return data


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


def get_gps_interval_extremes(gps, duration, event_type):
    """Get the extremes of the gps interval, based on interval duration and event type.

    Parameters
    ----------
    gps : int
        GPS of the event
    duration : int
        duration of the interval
    event_type : str
        event type

    Returns
    -------
    start : int
        left extreme of the interval
    end : int
        right extreme of the interval
    """
    if event_type not in defines.EVENT_LOCATION:
        raise ValueError("Event time can only be: {}".format(", ".join(defines.EVENT_LOCATION)))

    start = gps - duration // 2
    if event_type == "start":
        start = gps
    elif event_type == "end":
        start = gps - duration
    end = start + duration

    return start, end


def get_lock_channel_name_for_ifo(ifo):
    """Get the name of the lock channel for the specific interferometer.

    Parameters
    ----------
    ifo : str
        interferometer id

    Returns
    -------
    lock_ch_name : str
        name of the lock channel for the specific interferometer
    """
    if ifo.startswith("V"):
        return defines.LCK_CH_VIRGO
    elif ifo.startswith("L"):
        return defines.LCK_CH_LIGO
    elif ifo.startswith("H"):
        return defines.LCK_CH_HANFORD
    else:
        return None
