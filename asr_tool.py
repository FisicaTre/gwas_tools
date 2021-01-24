#  asr_tool.py - this file is part of the asr pagkage,
#  also known as "adaptive scattering recognizer".
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
import argparse
import numpy as np
from gwpy.timeseries import TimeSeriesDict
from scipy.signal import hilbert
from scipy.stats import pearsonr
import pytvfemd
import pickle
from .utils import file_utils
from .utils import signal_utils
from .common import defines


def run_tool(gps, target_channel_name, channels_file, out_path,
             fs, f_lowpass, n_scattering=1, smooth_win=50):
    # initialize variables
    LAMBDA = 1.064
    ch_f = open(channels_file, "r")
    channels_list = [ch.rstrip() for ch in ch_f.readlines() if ch.strip()]
    ch_f.close()
    data = np.array([], dtype=float)
    odir_name = ""

    # build time series matrix
    start_end = gps.split(",")
    if len(start_end) != 2:
        raise ValueError("GPS start or end time not provided.")

    gps_start = int(start_end[0])
    gps_end = int(start_end[1])
    odir_name = str(gps_start + defines.EXTRA_SECONDS) + "_" + str(gps_end - defines.EXTRA_SECONDS)
    data_dict = TimeSeriesDict.get([target_channel_name] + channels_list, gps_start, gps_end)
    dict_fs = np.min([data_dict[ch_name].channel.sample_rate.value for ch_name in channels_list])
    if dict_fs < fs:
        fs = dict_fs
    data_dict.resample(fs)
    data = np.zeros((data_dict[target_channel_name].value.shape[0], len(channels_list) + 1), dtype=float)
    data[:, 0] = data_dict[target_channel_name].value
    for i in range(1, len(channels_list) + 1):
        data[:, i] = data_dict[channels_list[i - 1]].value

    # create folder for results if it does not exist
    out_path = os.path.join(out_path, odir_name)
    if not os.path.isdir(out_path):
        os.makedirs(out_path, exist_ok=True)

    # predictors
    time = np.arange(0, len(data) / fs, 1 / fs, dtype=float)
    v_mat = np.diff(data[:, 1:], axis=0) / np.diff(time)[:, np.newaxis]
    for i in range(len(v_mat[0, :])):
        v_mat[:, i] = signal_utils.smooth(v_mat[:, i], smooth_win)

    predictor = n_scattering * (2 / LAMBDA) * np.abs(v_mat)
    predictor = predictor[(defines.EXTRA_SECONDS * fs):-(defines.EXTRA_SECONDS * fs), :]
    #for i in range(predictor.shape[1]):
    #    predictor[:, i] = butter_lowpass_filter(predictor[:, i], f_lowpass, fs)

    # target channel
    target_channel = signal_utils.butter_lowpass_filter(data[:, 0], f_lowpass, fs)
    #target_channel = signal_utils.lowpass(data[:, 0], f_lowpass)

    # tvf-emd
    imfs = pytvfemd.tvfemd(target_channel, MODES=1)
    imfs = imfs[(defines.EXTRA_SECONDS * fs):-(defines.EXTRA_SECONDS * fs), :]
    imfs = (imfs - np.nanmean(imfs, axis=0)) / np.nanstd(imfs, axis=0)

    # correlations
    corrs = np.zeros((imfs.shape[1], data.shape[1] - 1), dtype=float)
    envelopes = np.zeros((imfs.shape[0] - 1, imfs.shape[1]), dtype=float)
    for k in range(imfs.shape[1]):
        analytic_signal = hilbert(imfs[:, k])
        upper_env = np.abs(analytic_signal)[1:]
        upper_env = signal_utils.smooth(upper_env, smooth_win)
        envelopes[:, k] = upper_env

        for l in range(len(channels_list)):
            r_corr, _ = pearsonr(predictor[:, l], upper_env)
            corrs[k, l] = -999.0 if np.isnan(r_corr) else r_corr

    # save instantaneous amplitudes
    f_imfs = open(os.path.join(out_path, "envelopes_{}.imf".format(target_channel_name)), "wb")
    pickle.dump(envelopes, f_imfs)
    f_imfs.close()

    # max correlations
    max_vals = np.max(corrs, axis=1)
    max_channels = np.argmax(corrs, axis=1)
    max_channel = int(np.argmax(max_vals))
    try:
        max_ch_str = channels_list[max_channels[max_channel]]
        mean_freq = signal_utils.mean_frequency(data[:, max_channels[max_channel] + 1], fs)
    except:
        max_ch_str = "Not found"
        mean_freq = 0

    # output file
    out_file = open(os.path.join(out_path, "output.yml"), "w")
    out_file.write("{}:\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n\n".format(defines.PARAMS_SECT_KEY, defines.GPS_KEY, gps, defines.TARGET_CH_KEY,
                                         target_channel_name, defines.CH_LIST_KEY, channels_file, defines.OUT_PATH_KEY,
                                         out_path, defines.SAMP_FREQ_KEY, fs, defines.LOWPASS_FREQ_KEY, f_lowpass,
                                         defines.SCATTERING_KEY, n_scattering, defines.SMOOTH_WIN_KEY, smooth_win))
    out_file.write("{}:\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n"
                   "  {}: {}\n\n".format(defines.MAX_CORR_SECT_KEY, defines.IMF_KEY, max_channel + 1,
                                         defines.CHANNEL_KEY, max_ch_str, defines.CORR_KEY,
                                         max_vals[max_channel], defines.MEAN_FREQ_KEY, mean_freq))
    out_file.write("{}:\n".format(defines.CORR_SECT_KEY))
    for n_imf in range(len(max_vals)):
        try:
            ch_str = channels_list[max_channels[n_imf]]
            m_fr = signal_utils.mean_frequency(data[:, max_channels[n_imf] + 1], fs)
        except:
            ch_str = "Not found"
            m_fr = 0

        out_file.write("  - {}: {}\n"
                       "    {}: {}\n"
                       "    {}: {}\n"
                       "    {}: {}\n".format(defines.IMF_KEY, n_imf + 1,
                                             defines.CHANNEL_KEY, ch_str,
                                             defines.CORR_KEY, max_vals[n_imf],
                                             defines.MEAN_FREQ_KEY, m_fr))
    out_file.close()

    # save predictors
    f_pred = open(os.path.join(out_path, "{}.predictors".format(target_channel_name)), "wb")
    pickle.dump(predictor[:, max_channels], f_pred)
    f_pred.close()
