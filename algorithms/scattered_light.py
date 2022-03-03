#  scattered_light.py - this file is part of the gwadaptive_scattering package.
#  Copyright (C) 2020- Stefano Bianchi
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
from ..utils import signal_utils, file_utils
from ..common import defines
from gwpy.timeseries import TimeSeriesDict
from gwpy.io import datafind as io_datafind


def scattered_light(gps, seconds, target_channel_name, channels_file, out_path, f_lowpass,
                    event="center", fs=256, n_scattering=1, smooth_win=50):
    """Analysis for scattered light identification.
    The script outputs a folder named as the input gps,
    with inside three files:
        - target channel instantaneous amplitudes (*.imf)
        - most correlated predictor for each imf (*.predictors)
        - output.yml, a summary of the analysis' results

    Parameters
    ----------
    gps : int
        gps of the event
    seconds : int
        how many seconds to analyze in total
    target_channel_name : str
        target channel name
    channels_file : str
        channels list
    out_path : str
        output path where to save results
    f_lowpass : float or str
        lowpass filter frequency
    event : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    fs : float, optional
        channels resample frequency (default : 256)
    n_scattering : int, optional
        number of signal bounces (default : 1)
    smooth_win : int, optional
        signals smoothing window (default : 50)
    """
    if event not in defines.EVENT_LOCATION:
        raise ValueError("Event time can only be: {}".format(", ".join(defines.EVENT_LOCATION)))
    if not isinstance(f_lowpass, int) or not isinstance(f_lowpass, float) or f_lowpass not in defines.LOWP_FREQ_OPTS:
        raise ValueError("Lowpass frequency must be a float or one of these "
                         "strings : {}".format(", ".join(defines.LOWP_FREQ_OPTS)))

    # initialize variables
    ch_f = open(channels_file, "r")
    channels_list = [ch.rstrip() for ch in ch_f.readlines() if ch.strip()]
    ch_f.close()

    # create folder for results if it does not exist
    odir_name = "{:d}".format(gps)
    out_path = os.path.join(out_path, odir_name)
    if not os.path.isdir(out_path):
        os.makedirs(out_path, exist_ok=True)
        
    # build time series matrix
    gps_start, gps_end = signal_utils.get_gps_interval_extremes(gps, seconds, event)
    lock_channel_name = signal_utils.get_lock_channel_name_for_ifo(signal_utils.get_ifo_of_channel(target_channel_name))
    if lock_channel_name is not None:
        lock_channel_data = signal_utils.get_instrument_lock_data(lock_channel_name, gps_start, gps_end)
        if not np.all(lock_channel_data == 1):
            return None

    data, fs = signal_utils.get_data_from_time_series_dict(target_channel_name, channels_list,
                                                           gps_start, gps_end, fs, verbose=True)

    # predictors
    predictor = signal_utils.get_predictors(data[:, 1:], fs, smooth_win=smooth_win, n_scattering=n_scattering)

    # compute lowpass frequency in case it is a string
    if isinstance(f_lowpass, str):
        if f_lowpass == "average":
            f_lowpass = np.max([np.mean(predictor[:, i]) for i in range(predictor.shape[1])])
        elif f_lowpass == "max":
            f_lowpass = np.max([np.max(predictor[:, i]) for i in range(predictor.shape[1])])

    # target channel
    target_channel = signal_utils.butter_lowpass_filter(data[:, 0], f_lowpass, fs)

    # tvf-emd
    imfs = signal_utils.get_imfs(target_channel, fs, max_imf=1)

    # correlations
    corrs = np.zeros((imfs.shape[1], predictor.shape[1]), dtype=float)
    envelopes = np.zeros((imfs.shape[0] - 1, imfs.shape[1]), dtype=float)
    for k in range(imfs.shape[1]):
        upper_env = signal_utils.upper_envelope(imfs[:, k])[1:]
        upper_env = signal_utils.smooth(upper_env, smooth_win)
        envelopes[:, k] = upper_env
        for l in range(predictor.shape[1]):
            corrs[k, l] = signal_utils.get_correlation_between(predictor[:, l], upper_env)

    file_utils.save_envelopes(envelopes, "_".join(target_channel_name.split(":")), out_path)

    # max correlations
    max_vals = np.max(corrs, axis=1)
    max_channels = np.argmax(corrs, axis=1)
    # if len(channels_list) > 1:
    #    max_vals_2 = [np.max([n for n in corrs[i, :] if n != np.max(corrs[i, :])]) for i in range(corrs.shape[0])]
    #    max_channels_2 = [np.where(corrs[i, :] == max_vals_2[i])[0][0] for i in range(corrs.shape[0])]

    # output file
    out_file = file_utils.YmlFile()
    out_file.write_parameters(gps, seconds, event, target_channel_name, channels_file,
                              out_path, fs, f_lowpass, n_scattering, smooth_win)
    ch_str = []
    ch_corr = []
    ch_m_fr = []
    for n_imf in range(len(max_vals)):
        try:
            ch_str.append(channels_list[max_channels[n_imf]])
        except:
            ch_str.append("Not found")
        try:
            ch_corr.append(max_vals[n_imf])
        except:
            ch_corr.append(-999.0)
        try:
            ch_m_fr.append(signal_utils.mean_frequency(channels_list[max_channels[n_imf]], gps_start,
                                                       gps_end, bandpass_limits=(0.03, 10)))
        except:
            ch_m_fr.append(0.0)
    out_file.write_correlation_section(ch_str, ch_corr, ch_m_fr)

    # if len(channels_list) > 1:
    #    ch_2_str = []
    #    ch_2_corr = []
    #    ch_2_m_fr = []
    #    for n_imf in range(len(max_vals_2)):
    #        try:
    #            ch_2_str.append(channels_list[max_channels_2[n_imf]])
    #            ch_2_corr.append(max_vals_2[n_imf])
    #            ch_2_m_fr.append(signal_utils.mean_frequency(data[:, max_channels_2[n_imf] + 1], fs))
    #        except:
    #            ch_2_str.append("Not found")
    #            ch_2_corr.append(-999.0)
    #            ch_2_m_fr.append(0.0)
    #    out_file.write_2nd_best_correlation_section(ch_2_str, ch_2_corr, ch_2_m_fr)

    frametype = io_datafind.find_frametype("L1:ISI-GND_STS_ETMX_X_BLRMS_100M_300M", gpstime=gps)
    seism_channels = ["L1:ISI-GND_STS_ETMX_X_BLRMS_100M_300M",
                      "L1:ISI-GND_STS_ETMY_Y_BLRMS_100M_300M",
                      "L1:ISI-GND_STS_ETMX_Z_BLRMS_100M_300M",
                      "L1:ISI-GND_STS_ETMX_X_BLRMS_30M_100M",
                      "L1:ISI-GND_STS_ETMX_Y_BLRMS_30M_100M",
                      "L1:ISI-GND_STS_ETMX_Z_BLRMS_30M_100M",
                      "L1:ISI-GND_STS_ETMX_X_DQ",
                      "L1:ISI-GND_STS_ETMX_Y_DQ",
                      "L1:ISI-GND_STS_ETMX_Z_DQ"]
    try:
        seismometers = TimeSeriesDict.get(seism_channels, gps_start, gps_end, verbose=True, frametype=frametype,
                                          resample=3)
    except:
        seismometers = TimeSeriesDict.get(seism_channels, gps_start, gps_end, verbose=True, frametype=frametype)
    seis_dict = {}
    for s in seism_channels:
        seis_dict[s.split(":")[1]] = seismometers[s].value.mean()
    out_file.write_seismic_channels(seis_dict)

    out_file.save(out_path)

    # save predictors
    file_utils.save_predictors(predictor[:, max_channels], "_".join(target_channel_name.split(":")), out_path)
