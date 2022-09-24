#  plot_utils.py - this file is part of the gwadaptive_scattering package.
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
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from astropy.time import Time
from gwpy.timeseries import TimeSeries
from gwpy.segments import Segment
from ..utils import file_utils
from ..common import defines


def __plot_imfs_arg(arg):
    arg_ok = True
    if arg is not None:
        if isinstance(arg, list):
            for a in arg:
                if not isinstance(a, int):
                    arg_ok = False
                    break
        else:
            arg_ok = False

    return arg_ok


def __corr_thr(arg):
    arg = float(arg)
    if -1.0 <= arg <= 1.0:
        return arg
    else:
        if arg < 0.0:
            return -1.0
        else:
            return 1.0


def plot_imf(pred, pred_name, imf_ia, imf_ia_name, gps, samp_freq,
             title, plot_name, event_time="center", figsize=None):
    """Plot of the instantaneous amplitude and predictor.
    
    Parameters
    ----------
    pred : numpy array
        predictor
    pred_name : str
        predictor name
    imf_ia : numpy array
        imf instantaneous amplitude
    imf_ia_name : str
        imf instantaneous amplitude name
    gps : int
        gps of the event
    samp_freq : float
        sampling frequency
    title : str
        plot title
    plot_name : str
        plot name (full path)
    event_time : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    figsize : tuple[int], optional
        figure size (default : None)
    """
    if event_time not in defines.EVENT_LOCATION:
        raise ValueError("Event time can only be: {}".format(", ".join(defines.EVENT_LOCATION)))

    t1 = Time(gps, format="gps")
    t2 = Time(t1, format="iso", scale="utc")
    gps_date = "Time [seconds] from {} UTC ({:.1f})\n".format(t2.value.split(".")[0], t1.value)

    x = np.arange(0, len(pred) / samp_freq, 1 / samp_freq, dtype=float)
    if event_time == "center":
        x = np.arange(-len(pred) / (2.0 * samp_freq), len(pred) / (2.0 * samp_freq), 1 / samp_freq, dtype=float)
    elif event_time == "end":
        x = np.arange(-len(pred) / samp_freq, 0, 1 / samp_freq, dtype=float)

    font_size = 16
    fig, ax1 = plt.subplots() if figsize is None else plt.subplots(figsize=figsize)
    ax2 = ax1.twinx()
    l1 = ax1.plot(x[:len(pred)], pred, "r-", label=pred_name)
    l2 = ax2.plot(x[:len(imf_ia)], imf_ia, "b-", label=imf_ia_name)
    ls = l1 + l2
    leg = plt.legend(ls, [l.get_label() for l in ls], ncol=1, loc=2)
    try:
        plt.draw()
        leg_height_ratio = leg.get_window_extent().height / ax1.get_window_extent().height
    except:
        leg_height_ratio = 0.25
    
    ax1.set_ylim(0, np.max(pred) + (np.max(pred) - np.min(pred)) * (leg_height_ratio + 0.1))
    ax2.set_ylim(0, np.max(imf_ia) + (np.max(imf_ia) - np.min(imf_ia)) * (leg_height_ratio + 0.1))
    ax1.set_ylabel("Predictor [$Hz$]", color="r", fontsize=font_size)
    ax1.set_xlabel(gps_date, fontsize=font_size)
    ax2.set_ylabel("IA", color="b", fontsize=font_size)
    ax1.grid(True)
    ax2.grid(False)
    plt.title(title)
    plt.savefig(plot_name, bbox_inches="tight", dpi=300)
    plt.close("all")
            
            
def plot_omegagram_download(pred, pred_name, target_name, gps, seconds, plot_name,
                            norm=False, harmonics=None, event_time="center", figsize=None):
    """Omegagram plot with download of the target channel.
    
    Parameters
    ----------
    pred : numpy array
        predictor
    pred_name : str
        name of the predictor channel
    target_name : str
        name of the channel from which compute the omegagram
    gps : int
        gps of the event
    seconds : int
        how many seconds to plot
    plot_name : str
        plot name (full path)
    norm : bool, optional
        normalize predictor to 10 Hz (default : False)
    harmonics : list[int], optional
        predictor harmonics to be plotted (default : [1, 2, 3, 4, 5])
    event_time : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    figsize : tuple[int], optional
        figure size
    """
    if harmonics is None:
        harmonics = [1, 2, 3, 4, 5]

    if event_time not in defines.EVENT_LOCATION:
        raise ValueError("Event time can only be: {}".format(", ".join(defines.EVENT_LOCATION)))

    gps1 = gps
    gps2 = gps + seconds
    if event_time == "center":
        gps1 = gps - seconds // 2
        gps2 = gps + seconds // 2
    elif event_time == "end":
        gps1 = gps - seconds
        gps2 = gps
    if norm:
        correction_factor = 10.0 / np.max(pred)
    else:
        correction_factor = 1.0
            
    line_width = 1
    colormap = "viridis"
    plot_t_min = 0
    plot_t_max = gps2 - gps1
    if event_time == "center":
        plot_t_min = (gps1 - gps2) // 2
        plot_t_max = (gps2 - gps1) // 2
    elif event_time == "end":
        plot_t_min = gps1 - gps2
        plot_t_max = 0

    t1 = Time(gps, format="gps")
    t2 = Time(t1, format="iso", scale="utc")
    gps_date = "Time [seconds] from {} UTC ({:.1f})\n".format(t2.value.split(".")[0], t1.value)

    ts_l2 = pred * correction_factor
    ts = TimeSeries.get(target_name, gps1, gps2).astype("float64")
    if ":" in target_name and target_name.split(":")[0][0] == "V":
        ts = TimeSeries(ts.value, times=np.arange(gps1, gps2 + ts.dt.value, ts.dt.value, dtype=float), dtype=float)
    ts.times = ts.times.value - gps
    tsq = ts.q_transform(outseg=Segment(plot_t_min, plot_t_max + 0.2), tres=0.2, fres=0.2, whiten=True)

    plot_f_min = tsq.yindex[0].value
    plot_f_max = np.max([tsq.yindex[-1].value, np.max(ts_l2) * np.max(harmonics)])
    if plot_f_max > 50:
        plot_f_max = 50
    
    plot = plt.figure(figsize=figsize)
    ax = plot.gca()
    pcm = ax.imshow(tsq, vmin=0, vmax=15, cmap=colormap)
    ax.set_ylabel("Frequency [Hz]", fontsize=16)
    ax.grid(False)
    for i in harmonics:
        ax.plot(np.linspace(plot_t_min, plot_t_max, len(ts_l2)), ts_l2 * i, color="r", lw=line_width)
    ax.set_xlim(plot_t_min, plot_t_max)
    ax.set_ylim(plot_f_min, plot_f_max)
    ax.set_xlabel(gps_date, fontsize=16)
    ax.set_title("{} | Fringes: {}".format(target_name, pred_name), fontsize=16)
    cbar = ax.colorbar(clim=(0, 15), location="right")
    cbar.set_label("Normalized energy", fontsize=16)
    plt.savefig(plot_name, bbox_inches="tight", dpi=300)
    plt.close("all")
        
        
def plot_omegagram(pred, pred_name, target, target_name, gps, seconds, fs, plot_name,
                   norm=False, harmonics=None, event_time="center", figsize=None):
    """Omegagram plot.
    
    Parameters
    ----------
    pred : numpy array
        predictor
    pred_name : str
        name of the predictor channel
    target : numpy array
        channel from which compute the omegagram
    target_name : str
        name of the target channel
    gps : int
        gps of the event
    seconds : int
        how many seconds to plot
    fs : int
        sampling frequency
    plot_name : str
        plot name (full path)
    norm : bool, optional
        normalize predictor to 10 Hz (default : False)
    harmonics : list[int], optional
        predictor harmonics to be plotted (default : [1, 2, 3, 4, 5])
    event_time : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    figsize : tuple[int], optional
        figure size
    """
    if harmonics is None:
        harmonics = [1, 2, 3, 4, 5]

    if event_time not in defines.EVENT_LOCATION:
        raise ValueError("Event time can only be: {}".format(", ".join(defines.EVENT_LOCATION)))

    gps1 = gps
    gps2 = gps + seconds
    if event_time == "center":
        gps1 = gps - seconds // 2
        gps2 = gps + seconds // 2
    elif event_time == "end":
        gps1 = gps - seconds
        gps2 = gps
    if norm:
        correction_factor = 10.0 / np.max(pred)
    else:
        correction_factor = 1.0
            
    line_width = 1
    colormap = "viridis"
    plot_t_min = 0
    plot_t_max = gps2 - gps1
    if event_time == "center":
        plot_t_min = (gps1 - gps2) // 2
        plot_t_max = (gps2 - gps1) // 2
    elif event_time == "end":
        plot_t_min = gps1 - gps2
        plot_t_max = 0

    t1 = Time(gps, format="gps")
    t2 = Time(t1, format="iso", scale="utc")
    gps_date = "Time [seconds] from {} UTC ({:.1f})\n".format(t2.value.split(".")[0], t1.value)

    ts_l2 = pred * correction_factor
    ts = TimeSeries(target, times=np.arange(gps1, gps1 + len(target) / fs, 1 / fs, dtype=float), dtype=float)
    ts.times = ts.times.value - gps
    tsq = ts.q_transform(outseg=Segment(plot_t_min, plot_t_max + 0.2), tres=0.2, fres=0.2, whiten=True)

    plot_f_min = tsq.yindex[0].value
    plot_f_max = np.max([tsq.yindex[-1].value, np.max(ts_l2) * np.max(harmonics)])
    if plot_f_max > 50:
        plot_f_max = 50
    
    plot = plt.figure(figsize=figsize)
    ax = plot.gca()
    pcm = ax.imshow(tsq, vmin=0, vmax=15, cmap=colormap)
    ax.set_ylabel("Frequency [Hz]", fontsize=16)
    ax.grid(False)
    for i in harmonics:
        ax.plot(np.linspace(plot_t_min, plot_t_max, len(ts_l2)), ts_l2 * i, color="r", lw=line_width)
    ax.set_xlim(plot_t_min, plot_t_max)
    ax.set_ylim(plot_f_min, plot_f_max)
    ax.set_xlabel(gps_date, fontsize=16)
    ax.set_title("{} | Fringes: {}".format(target_name, pred_name), fontsize=16)
    cbar = ax.colorbar(clim=(0, 15), location="right")
    cbar.set_label("Normalized energy", fontsize=16)
    plt.savefig(plot_name, bbox_inches="tight", dpi=300)
    plt.close("all")
        

def plot_imfs(folders, imfs_to_plot=None, imf_thr=-1.0, save_ext="png", figsize=None):
    """Plot imfs in the given folders.

    Parameters
    ----------
    folders : list[str]
        paths to the files needed for the plots
    imfs_to_plot : list[int]
        imfs to consider for combos (default : None, i.e. all the imfs)
    imf_thr : float, optional
        correlation value above which to plot imfs (default : -1.0)
    save_ext : str, optional
        plots extension (default : png)
    figsize : tuple[int], optional
        figure size
    """
    if not __plot_imfs_arg(imfs_to_plot):
        raise ValueError("`imfs` parameter must be a list of integers.")
    imf_thr = __corr_thr(imf_thr)

    for res_folder in folders:
        yf = file_utils.YmlFile(res_folder)
        ia = file_utils.load_imfs(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf.get_target_channel()
        gps_event = yf.get_gps()
        event_time = yf.get_event_position()
        fs = yf.get_sampling_frequency()

        if imfs_to_plot is None:
            for n_imf in range(yf.get_imfs_count()):
                if yf.get_corr_of_imf(n_imf + 1) >= imf_thr:
                    plot_imf(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1), ia[:, n_imf], target_channel,
                             gps_event, fs, "$\\rho$ = {:.2f}".format(yf.get_corr_of_imf(n_imf + 1)),
                             os.path.join(res_folder, file_utils.culprit_plot_name(n_imf + 1, save_ext)),
                             event_time=event_time, figsize=figsize)
        else:
            for n_imf in range(yf.get_imfs_count()):
                if n_imf + 1 in imfs_to_plot:
                    if yf.get_corr_of_imf(n_imf + 1) >= imf_thr:
                        plot_imf(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1), ia[:, n_imf],
                                 target_channel, gps_event, fs,
                                 "$\\rho$ = {:.2f}".format(yf.get_corr_of_imf(n_imf + 1)),
                                 os.path.join(res_folder, file_utils.culprit_plot_name(n_imf + 1, save_ext)),
                                 event_time=event_time, figsize=figsize)


def plot_omegagrams(folders, imfs_to_plot=None, omegagram_thr=-1.0, harmonics=None, save_ext="png", figsize=None):
    """Plot omegagrams in the given folders.

    Parameters
    ----------
    folders : list[str]
        paths to the files needed for the plots
    imfs_to_plot : list[int]
        imfs to consider for combos (default : None, i.e. all the imfs)
    omegagram_thr : float, optional
        correlation value above which to plot omegagrams (default : -1.0)
    harmonics : list[int]
        harmonics for the culprit (default : [1, 2, 3, 4, 5])
    save_ext : str, optional
        plots extension (default : png)
    figsize : tuple[int], optional
        figure size
    """
    if harmonics is None:
        harmonics = [1, 2, 3, 4, 5]

    if not __plot_imfs_arg(imfs_to_plot):
        raise ValueError("`imfs` parameter must be a list of integers.")
    omegagram_thr = __corr_thr(omegagram_thr)

    for res_folder in folders:
        yf = file_utils.YmlFile(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf.get_target_channel()
        gps = yf.get_gps()
        seconds = yf.get_seconds()
        event_time = yf.get_event_position()

        if imfs_to_plot is None:
            for n_imf in range(yf.get_imfs_count()):
                if yf.get_corr_of_imf(n_imf + 1) >= omegagram_thr:
                    plot_omegagram_download(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1),
                                            target_channel, gps, seconds,
                                            os.path.join(res_folder, file_utils.omegagram_plot_name(n_imf + 1, save_ext)),
                                            harmonics=harmonics, event_time=event_time, figsize=figsize)
        else:
            for n_imf in range(yf.get_imfs_count()):
                if n_imf + 1 in imfs_to_plot:
                    if yf.get_corr_of_imf(n_imf + 1) >= omegagram_thr:
                        plot_omegagram_download(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1),
                                                target_channel, gps, seconds,
                                                os.path.join(res_folder, file_utils.omegagram_plot_name(n_imf + 1, save_ext)),
                                                harmonics=harmonics, event_time=event_time, figsize=figsize)


def plot_combinations(folders, imfs_to_plot=None, save_ext="png", imf_thr=-1.0, figsize=None):
    """Plot of the sum of more instantaneous amplitudes and the predictor.

    Parameters
    ----------
    folders : list[str]
        paths to the files needed for the plots
    imfs_to_plot : list[int]
        imfs to consider for combos (default : None, i.e. all the imfs)
    imf_thr : float, optional
        correlation value above which to plot the combo (default : -1.0)
    save_ext : str, optional
        plots extension (default : png)
    figsize : tuple[int], optional
        figure size
    """
    if not __plot_imfs_arg(imfs_to_plot):
        raise ValueError("`imfs` parameter must be a list of integers.")
    imf_thr = __corr_thr(imf_thr)

    for res_folder in folders:
        yf = file_utils.YmlFile(res_folder)
        ia = file_utils.load_imfs(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf.get_target_channel()
        gps_event = yf.get_gps()
        event_time = yf.get_event_position()
        fs = yf.get_sampling_frequency()
        combos = yf.get_combos()

        if imfs_to_plot is None:
            for i in range(len(combos)):
                c = combos[i][defines.CORR_KEY]
                if c >= imf_thr:
                    imf_idxs = [j - 1 for j in combos[i][defines.IMF_KEY]]
                    sum_envelope = np.sum(ia[:, imf_idxs], axis=1)
                    plot_imf(preds[:, imf_idxs[0]], combos[i][defines.CHANNEL_KEY], sum_envelope,
                             target_channel + " (combo)", gps_event, fs, "$\\rho$ = {:.2f}".format(c),
                             os.path.join(res_folder, file_utils.combo_plot_name(combos[i][defines.IMF_KEY], save_ext)),
                             event_time=event_time, figsize=figsize)
        else:
            imfs_list = [imf for imf in imfs_to_plot if imf <= yf.get_imfs_count()]
            if len(imfs_list) > 0:
                for i in range(len(combos)):
                    is_ok = False
                    for imf_i in combos[i][defines.IMF_KEY]:
                        if imf_i in imfs_list:
                            is_ok = True
                            break
                    if is_ok:
                        c = combos[i][defines.CORR_KEY]
                        if c >= imf_thr:
                            imf_idxs = [j - 1 for j in combos[i][defines.IMF_KEY]]
                            sum_envelope = np.sum(ia[:, imf_idxs], axis=1)
                            plot_imf(preds[:, imf_idxs[0]], combos[i][defines.CHANNEL_KEY], sum_envelope,
                                     target_channel + " (combo)", gps_event, fs, "$\\rho$ = {:.2f}".format(c),
                                     os.path.join(res_folder,
                                                  file_utils.combo_plot_name(combos[i][defines.IMF_KEY], save_ext)),
                                     event_time=event_time, figsize=figsize)


def plot_seismic_data(folders_path, ifo, combo=True, save_ext="png"):
    """Plot seismic channels along with the correlations time series.

    Parameters
    ----------
    folders_path : str
        path to the folder containing the analysis folders
    ifo : str
        interferometer id
    combo : bool, optional
        if True, plot correlations of the combo,
        if False plot the correlation of the single imf (default : True)
    save_ext : str, optional
        plots extension (default : png)
    """
    # ***** NOTE : the plots refers (by now) only to the first IMF ***** #
    imf_to_plot = 1
    # ***** ***** #
    if ifo.startswith("L"):
        seism_channels = [sc.split(":")[1] for sc in defines.LIGO_SEISMIC_CHANNELS]
    elif ifo.startswith("H"):
        seism_channels = [sc.split(":")[1] for sc in defines.HANFORD_SEISMIC_CHANNELS]
    elif ifo.startswith("V"):
        seism_channels = [sc.split(":")[1] for sc in defines.VIRGO_SEISMIC_CHANNELS]
    else:
        raise ValueError("No seismic data for the current interferometer.")

    seismometers = {}
    seism_group = 3
    for i in range(len(seism_channels) // seism_group):
        seismometers[i] = {}
        for j in range(seism_group):
            seism_key = seism_channels[i * seism_group + j]
            seismometers[i][seism_key] = []

    corrs = []
    lock_periods = []
    res_folders = file_utils.get_results_folders(folders_path, filter_non_valid=False)
    out_folder = os.path.join(folders_path, defines.COMPARISON_FOLDER)
    if len(res_folders) > 0:
        gps = Time(int(os.path.split(res_folders[0])[1]), format="gps")
        date = Time(gps, format="iso", scale="utc").value.split(" ")[0]

        for res_folder in res_folders:
            if res_folder != out_folder:
                if file_utils.is_valid_folder(res_folder):
                    yf = file_utils.YmlFile(res_folder)
                    if combo:
                        try:
                            corrs.append(yf.get_corr_of_combo_with_imf(imf_to_plot))
                        except:
                            corrs.append(np.nan)
                    else:
                        corrs.append(yf.get_corr_of_imf(imf_to_plot))
                    seis_dict = yf.get_seismic_channels()
                    for k in seismometers.keys():
                        for sk in seis_dict.keys():
                            if sk in seismometers[k].keys():
                                seismometers[k][sk].append(seis_dict[sk])
                    lck = yf.is_locked()
                    if lck:
                        lock_periods.append(np.nan)
                    else:
                        lock_periods.append(1)
                else:
                    corrs.append(np.nan)
                    if file_utils.yml_exists(res_folder):
                        yf = file_utils.YmlFile(res_folder)
                        seis_dict = yf.get_seismic_channels()
                        for k in seismometers.keys():
                            for sk in seis_dict.keys():
                                if sk in seismometers[k].keys():
                                    seismometers[k][sk].append(seis_dict[sk])
                        lck = yf.is_locked()
                        if lck:
                            lock_periods.append(np.nan)
                        else:
                            lock_periods.append(1)
                    else:
                        for k in seismometers.keys():
                            for sk in seismometers[k].keys():
                                seismometers[k][sk].append(np.nan)
                        lock_periods.append(1)

        if len(corrs) > 0:
            maverage = pd.DataFrame({"Moving average": corrs})
            mavg = maverage.rolling(60, min_periods=1).mean().shift(-30)

            font_size = 24
            colors = ["b", "g", "orange"]
            labels = ["{:02d}".format(h) for h in range(24)]
            plt.rcParams.update({'font.size': font_size})
            plt.figure(figsize=(26, 14))

            n_plots = len(seismometers.keys()) + 1
            for i in range(n_plots):
                plt.subplot(n_plots, 1, i + 1)
                if i == 0:
                    plt.plot(range(len(corrs)), corrs, "b", label="$\\rho$")
                    plt.plot(range(len(corrs)), mavg.values, "r", label="moving average", linewidth=2)
                    for j in range(len(lock_periods) - 1):
                        if lock_periods[j] == 1:
                            plt.axvspan(j, j + 1, facecolor="0.2", alpha=0.5)
                    plt.xlim(0, len(corrs))
                    if np.isnan(corrs).all():
                        plt.ylim(0, 1)
                    else:
                        plt.ylim(np.nanmin(corrs) - 0.05, np.nanmax(corrs) + 0.05)
                    plt.ylabel("$\\rho$", fontsize=font_size)
                    plt.grid(True)
                    plt.xticks(np.arange(0, len(corrs), step=60), labels)
                    plt.legend(loc=0)
                else:
                    for j, ch in enumerate(seismometers[i - 1].keys()):
                        ts = seismometers[i - 1][ch]
                        plt.plot(range(len(ts)), ts, colors[j], label="{}:{}".format(ifo, ch))
                    plt.xticks(np.arange(0, len(seismometers[i - 1][seism_channels[(i - 1) * seism_group]]), step=60),
                               labels)
                    plt.xlim(0, len(seismometers[i - 1][seism_channels[(i - 1) * seism_group]]))
                    plt.legend(loc=0)
                    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
                    plt.ylabel("[$ms^{-1}$]", fontsize=font_size)
                    plt.grid(True)

            plt.xlabel("Hours from {} 00:00:00".format(date), fontsize=font_size)
            plt.savefig(os.path.join(out_folder, file_utils.seismic_plot_name(save_ext)), bbox_inches="tight", dpi=300)
            plt.close("all")


def plot_summaries(res_table, imf_thr=-1.0, save_ext="png", figsize=None):
    """Summary plots for the daily analysis.

    Parameters
    ----------
    res_table : str
        path to the summary table
    imf_thr : float, optional
        correlation value above which to consider a line in the summary table (default : -1.0)
    save_ext : str, optional
        plots extension (default : png)
    figsize : tuple[int], optional
        figure size
    """
    if not os.path.exists(res_table):
        return None

    # ***** NOTE : the plots refers (by now) only to the first IMF ***** #
    imf_to_plot = 1
    # ***** ***** #
    imf_thr = __corr_thr(imf_thr)
    table = pd.read_csv(res_table)
    # add range column
    table["range"] = "ELSE"
    freq_bands_names = ["{:.2f} Hz <= f < {:.2f} Hz".format(defines.FREQ_BANDS[i], defines.FREQ_BANDS[i + 1])
                        for i in range(len(defines.FREQ_BANDS) - 1)]
    for i in range(len(defines.FREQ_BANDS) - 1):
        table.range[(table[defines.summary_table_mean_frequency_column_of_imf(imf_to_plot)] >= defines.FREQ_BANDS[i]) &
                    (table[defines.summary_table_mean_frequency_column_of_imf(imf_to_plot)] < defines.FREQ_BANDS[i + 1])] = freq_bands_names[i]
    # add location column
    table["location"] = "OTHER"
    for k, v in defines.CHAMBERS.items():
        for vi in v:
            table.location[table[defines.summary_table_culprit_column_of_imf(imf_to_plot)].str.contains(vi, regex=False, na=False)] = k

    table = table[table[defines.summary_table_correlation_column_of_imf(imf_to_plot)] > imf_thr]

    font_size = 16
    plot_path = os.path.split(res_table)[0]
    # plot by frequency range
    freqs = table["range"].value_counts().sort_index()
    ranges = [">= {}\n{}".format(i.split("f")[0].strip()[:-3], i.split("f")[1].strip()) for i in list(freqs.index)]
    ranges_vals = freqs.values / len(table["range"])
    ranges_i = np.arange(len(ranges_vals))
    fig, ax = plt.subplots() if figsize is None else plt.subplots(figsize=figsize)
    plt.bar(ranges_i, ranges_vals, color="darkblue", edgecolor="black", linewidth=1)
    plt.xlabel("Frequency band", fontsize=font_size)
    plt.ylabel("Percentage", fontsize=font_size)
    plt.title("Percentage by frequency band (correlation > {:.2f})".format(imf_thr), fontsize=font_size)
    plt.xticks(ranges_i, ranges)
    plt.ylim(0, 1)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol="%", is_latex=False))
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, file_utils.summary_freq_plot_name(save_ext)), bbox_inches="tight", dpi=300)

    # plot by chamber
    locs = table["location"].value_counts()
    chambers = list(locs.index)
    chambers_vals = locs.values
    chambers_i = np.arange(len(chambers_vals))
    plt.subplots() if figsize is None else plt.subplots(figsize=figsize)
    plt.barh(chambers_i, chambers_vals, height=0.6, color="darkblue", edgecolor="black", linewidth=1)
    plt.title("Scattering glitches associated to location (correlation > {:.2f})".format(imf_thr), fontsize=font_size)
    plt.xlabel("Counts", fontsize=font_size)
    plt.ylabel("Location", fontsize=font_size)
    plt.yticks(chambers_i, chambers)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_path, file_utils.summary_chamber_plot_name(save_ext)), dpi=300, bbox_inches="tight")

    plt.close("all")
