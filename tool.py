#  tool.py - this file is part of the gwasr package.
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


import os
import numpy as np
from gwpy.timeseries import TimeSeriesDict
from scipy.signal import hilbert
from scipy.stats import pearsonr
import pytvfemd
import pickle
from .utils import signal_utils
from .common import defines


LAMBDA = 1.064


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
    verbose : verbosity option of TimeSeriesDict (default : False)

    Returns
    -------
    numpy ndarray
        matrix with channels values, first column
        corresponds to the target channel
    float
        common sampling frequency of the channels
        in the matrix
    """
    data_dict = TimeSeriesDict.get([target_channel_name] + channels_list,
                                   gps_start - defines.EXTRA_SECONDS,
                                   gps_end + defines.EXTRA_SECONDS,
                                   verbose=verbose)
    dict_fs = np.min([data_dict[ch_name].channel.sample_rate.value for ch_name in channels_list])
    if dict_fs < fs:
        fs = dict_fs
    data_dict.resample(fs)
    data = np.zeros((data_dict[target_channel_name].value.shape[0], len(channels_list) + 1), dtype=float)
    data[:, 0] = data_dict[target_channel_name].value
    for i in range(1, len(channels_list) + 1):
        data[:, i] = data_dict[channels_list[i - 1]].value

    return data, fs


def get_predictors(channels, fs, smooth_win=None, n_scattering=1):
    """Get predictors from channels.

    Parameters
    ----------
    channels : numpy ndarray
        channels matrix (shape (len(channels_values), num_channels)
    fs : float
        channels common sampling frequency
    smooth_win : int
        smoothing window (default : None)
    n_scattering : int
        signal bounces number

    Returns
    -------
    numpy ndarray
        channels predictors
    """
    time = np.arange(0, len(channels) / fs, 1 / fs, dtype=float)
    v_mat = np.diff(channels, axis=0) / np.diff(time)[:, np.newaxis]
    if smooth_win is not None:
        for i in range(len(v_mat[0, :])):
            v_mat[:, i] = signal_utils.smooth(v_mat[:, i], smooth_win)

    predictor = n_scattering * (2 / LAMBDA) * np.abs(v_mat)
    fs_int = int(fs)
    predictor = predictor[(defines.EXTRA_SECONDS * fs_int):-(defines.EXTRA_SECONDS * fs_int), :]
    # for i in range(predictor.shape[1]):
    #     predictor[:, i] = butter_lowpass_filter(predictor[:, i], f_lowpass, fs)

    return predictor


def get_imfs(target_channel, fs, norm=True):
    """Get imfs with pytvfemd.

    Parameters
    ----------
    target_channel : numpy array
        channel from which extract imfs
    fs : float
        `target_channel` sampling frequency
    norm : bool
        normalize imfs (default : True)

    Returns
    -------
    numpy ndarray
        `target_channel` imfs matrix
    """
    imfs = pytvfemd.tvfemd(target_channel, MODES=1)
    fs_int = int(fs)
    imfs = imfs[(defines.EXTRA_SECONDS * fs_int):-(defines.EXTRA_SECONDS * fs_int), :]
    if norm:
        imfs = (imfs - np.nanmean(imfs, axis=0)) / np.nanstd(imfs, axis=0)

    return imfs


def get_correlations(imfs, predictors, smooth_win=None, save_envelopes=True,
                     save_env_path=".", save_env_name="upper"):
    """Get correlations between imfs upper envelopes and predictors.

    Parameters
    ----------
    imfs : numpy ndarray
        imfs matrix
    predictors : numpy ndarray
        predictors matrix
    smooth_win : int
        smoothing window for the upper envelopes (default : None)
    save_envelopes : bool
        whether to save upper envelopes (default : True)
    save_env_path : str
        where to save envelopes (default : `.`, ignored if `save_envelopes` is False)
    save_env_name : str
        envelopes file name (default : upper, ignored if `save_envelopes` is False)

    Returns
    -------
    numpy ndarray
        correlation matrix (row: imf, col: predictor)
    """
    corrs = np.zeros((imfs.shape[1], predictors.shape[1]), dtype=float)
    envelopes = np.zeros((imfs.shape[0] - 1, imfs.shape[1]), dtype=float)
    for k in range(imfs.shape[1]):
        analytic_signal = hilbert(imfs[:, k])
        upper_env = np.abs(analytic_signal)[1:]
        if smooth_win is not None:
            upper_env = signal_utils.smooth(upper_env, smooth_win)
        if save_envelopes:
            envelopes[:, k] = upper_env

        for l in range(predictors.shape[1]):
            r_corr, _ = pearsonr(predictors[:, l], upper_env)
            corrs[k, l] = -999.0 if np.isnan(r_corr) else r_corr

    if save_envelopes:
        f_imfs = open(os.path.join(save_env_path, "envelopes_{}.imf".format(save_env_name)), "wb")
        pickle.dump(envelopes, f_imfs)
        f_imfs.close()

    return corrs


def save_predictors(preds, file_name, out_path):
    """Save predictors to binary file.

    Parameters
    ----------
    preds : numpy ndarray
        predictors
    file_name : str
        name of the file
    out_path : str
        where to save the file
    """
    f_pred = open(os.path.join(out_path, "{}.predictors".format(file_name)), "wb")
    pickle.dump(preds, f_pred)
    f_pred.close()

