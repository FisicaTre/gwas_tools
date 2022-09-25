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


import numpy as np
from ..utils import signal_utils, file_utils
from ..common import defines
from gwpy.timeseries import TimeSeriesDict
from gwpy.io import datafind as io_datafind


def scattered_light(gps, seconds, target_channel_name, channels_file, out_path, f_lowpass,
                    event="center", fs=256, n_scattering=1, smooth_win=50, max_imf=1,
                    combos=False, seismic=False, save_data=True, check_lock=False):
    """Analysis for scattered light identification.
    The script outputs a folder named as the input gps,
    with inside three files:
        - target channel instantaneous amplitudes (.imf, optional)
        - most correlated predictor for each imf (.predictors, optional)
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
        Can be `start`, `center`, or `end` (default : `center`)
    fs : float, optional
        channels resample frequency (default : 256)
    n_scattering : int, optional
        number of signal bounces (default : 1)
    smooth_win : int, optional
        signals smoothing window (default : 50)
    max_imf :int, optional
        maximum number of imfs to extract (default : 1)
    combos : bool, optional
        if True, combos are computed (default : False)
    seismic : bool, optional
        if True, seismic channels data are saved to output file (default : False)
    save_data : bool, optional
        if True, instantaneous amplitudes and predictors are saved to file (defaults : True)
    check_lock : bool, optional
        if True, lock channel status is checked, and if it is not always locked, the analysis is not performed (default : False)
    """
    if event not in defines.EVENT_LOCATION:
        raise ValueError("Event time can only be: {}".format(", ".join(defines.EVENT_LOCATION)))
    if f_lowpass not in defines.LOWP_FREQ_OPTS:
        try:
            f_lowpass = float(f_lowpass)
        except:
            raise ValueError("Lowpass frequency must be a float or one of these "
                             "strings : {}".format(", ".join(defines.LOWP_FREQ_OPTS)))

    # initialize variables
    ch_f = open(channels_file, "r")
    channels_list = [ch.rstrip() for ch in ch_f.readlines() if ch.strip()]
    ch_f.close()

    gps_start, gps_end = signal_utils.get_gps_interval_extremes(gps, seconds, event)
    ifo = signal_utils.get_ifo_of_channel(target_channel_name)

    # output file
    out_file = file_utils.YmlFile()

    if check_lock:
        lock_channel_name = signal_utils.get_lock_channel_name_for_ifo(ifo)
        if lock_channel_name is not None:
            lock_channel_data = signal_utils.get_instrument_lock_data(lock_channel_name, gps_start, gps_end)
            if ifo.startswith("L") or ifo.startswith("H"):
                if len(lock_channel_data) != 1 or lock_channel_data[0][0] != gps_start or lock_channel_data[0][-1] != gps_end:
                    out_file.write_lock_info(False)
                else:
                    out_file.write_lock_info(True)
            elif ifo.startswith("V"):
                if not np.all(lock_channel_data >= defines.VIRGO_SCIENCE_MODE_THR):
                    out_file.write_lock_info(False)
                else:
                    out_file.write_lock_info(True)

    # build time series matrix
    data = signal_utils.get_data_from_time_series_dict(target_channel_name, channels_list,
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
    imfs = signal_utils.get_imfs(target_channel, fs, max_imf=max_imf)

    # correlations
    corrs = np.zeros((imfs.shape[1], predictor.shape[1]), dtype=float)
    envelopes = np.zeros((imfs.shape[0] - 1, imfs.shape[1]), dtype=float)
    for k in range(imfs.shape[1]):
        upper_env = signal_utils.upper_envelope(imfs[:, k])[1:]
        upper_env = signal_utils.smooth(upper_env, smooth_win)
        envelopes[:, k] = upper_env
        for l in range(predictor.shape[1]):
            corrs[k, l] = signal_utils.get_correlation_between(predictor[:, l], upper_env)

    if save_data:
        # save instantaneous amplitudes
        file_utils.save_envelopes(envelopes, "_".join(target_channel_name.split(":")), out_path)

    # max correlations
    max_vals = np.max(corrs, axis=1)
    max_channels = np.argmax(corrs, axis=1)
    # if len(channels_list) > 1:
    #    max_vals_2 = [np.max([n for n in corrs[i, :] if n != np.max(corrs[i, :])]) for i in range(corrs.shape[0])]
    #    max_channels_2 = [np.where(corrs[i, :] == max_vals_2[i])[0][0] for i in range(corrs.shape[0])]

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

    selected_predictors = predictor[:, max_channels]
    if combos:
        combos_imfs, combos_channels, combos_corrs = signal_utils.get_combos(ch_str, envelopes, selected_predictors)
        out_file.write_combo_section(combos_imfs, combos_channels, combos_corrs)

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

    if seismic:
        if ifo.startswith("L"):
            frametype = io_datafind.find_frametype(defines.LIGO_SEISMIC_CHANNELS[0], gpstime=gps)
            if frametype is not None:
                seismometers = TimeSeriesDict.get(defines.LIGO_SEISMIC_CHANNELS, gps_start, gps_end, verbose=True,
                                                  frametype=frametype, resample=3)
            else:
                seismometers = TimeSeriesDict.get(defines.LIGO_SEISMIC_CHANNELS, gps_start, gps_end, verbose=True)
                seismometers.resample(3)
            seis_dict = {}
            for s in defines.LIGO_SEISMIC_CHANNELS:
                seis_dict[s.split(":")[1]] = seismometers[s].value.mean()
            out_file.write_seismic_channels(seis_dict)
        elif ifo.startswith("H"):
            frametype = io_datafind.find_frametype(defines.HANFORD_SEISMIC_CHANNELS[0], gpstime=gps)
            if frametype is not None:
                seismometers = TimeSeriesDict.get(defines.HANFORD_SEISMIC_CHANNELS, gps_start, gps_end, verbose=True,
                                                  frametype=frametype, resample=3)
            else:
                seismometers = TimeSeriesDict.get(defines.HANFORD_SEISMIC_CHANNELS, gps_start, gps_end, verbose=True)
                seismometers.resample(3)
            seis_dict = {}
            for s in defines.HANFORD_SEISMIC_CHANNELS:
                seis_dict[s.split(":")[1]] = seismometers[s].value.mean()
            out_file.write_seismic_channels(seis_dict)
        elif ifo.startswith("V"):
            from gwdama.io import GwDataManager as gwdm
            seismometers = gwdm.read_from_virgo(defines.VIRGO_SEISMIC_CHANNELS, gps_start,
                                                gps_end, ffl_spec="V1trend100")
            # seismometers.resample(3)
            seis_dict = {}
            for s in defines.VIRGO_SEISMIC_CHANNELS:
                seis_dict[s.split(":")[1]] = seismometers[s].value.mean()
            out_file.write_seismic_channels(seis_dict)

    out_file.save(out_path)

    if save_data:
        # save predictors
        file_utils.save_predictors(selected_predictors, "_".join(target_channel_name.split(":")), out_path)
