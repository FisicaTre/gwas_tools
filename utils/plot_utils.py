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

# TODO : replace plot names with functions from file_utils

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from gwpy.timeseries import TimeSeries
from gwpy.segments import Segment
from ..utils import file_utils
from ..common import defines


def __plot_imfs_arg(arg):
    if arg == "all":
        return arg
    elif isinstance(arg, list):
        return arg
    else:
        return None


def __corr_thr(arg):
    arg = float(arg)
    if -1.0 <= arg <= 1.0:
        return arg
    else:
        if arg < 0.0:
            return -1.0
        else:
            return 1.0


def plot_imf(pred, pred_name, imf_ia, imf_ia_name, gps, samp_freq, title,
             plot_name, save_path, event_time="center", save_ext="png", figsize=None):
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
        plot name
    save_path : str
        save path
    event_time : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    save_ext : str, optional
        plot extension (default : png)
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
    plt.savefig(os.path.join(save_path, plot_name + "." + save_ext), bbox_inches="tight", dpi=300)
    plt.close("all")
    
    
def plot_combinations(plot_combos, ias, predictors, target_channel_name, gps1, samp_freq,
                      out_path, event_time="center", save_ext="png", thr=-1.0, figsize=None):
    """Plot sum of more instantaneous amplitudes and the predictor.
    
    Parameters
    ----------
    plot_combos : dict
        combos dict, as retrieved from yml file
    ias : numpy array
        imf instantaneous amplitudes matrix
    predictors : numpy array
        predictors matrix
    target_channel_name : str
        name of the channel from which `ias` are computed
    gps1 : int
        gps of the event
    samp_freq : float
        sampling frequency
    out_path : str
        save path
    event_time : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    save_ext : str, optional
        plot extension (default : png)
    thr : float, optional
        correlation threshold for plots (default : -1.0)
    figsize : tuple[int], optional
        figure size
    """
    # seen = set()
    # uniq = [x for x in plot_channels if x not in seen and not seen.add(x)]
    # for u in uniq:
    #    plot_idxs = np.where(plot_channels == u)[0]
    #    if len(plot_idxs) != 1:

    for idx, imf_list in enumerate(plot_combos[defines.IMF_KEY]):
        imf_idxs = [i - 1 for i in imf_list]
        sum_envelope = np.sum(ias[:, imf_idxs], axis=1)
        c = plot_combos[defines.CORR_KEY][idx]
        if c > thr:
            plot_imf(predictors[:, imf_idxs[0]], plot_combos[defines.CHANNEL_KEY][idx],
                     sum_envelope, target_channel_name + " (combo)",
                     gps1, samp_freq, "$\\rho$ = {:.2f}".format(c),
                     "combo_imf_{}_culprit".format("+".join(imf_list)), out_path,
                     event_time=event_time, save_ext=save_ext, figsize=figsize)
            
            
def plot_omegagram_download(pred, pred_name, target_name, gps, seconds, plot_name, save_path,
                            norm=False, harmonics=None, event_time="center", save_ext="png",
                            figsize=None):
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
        plot name
    save_path : str
        save path
    norm : bool, optional
        normalize predictor to 10 Hz (default : False)
    harmonics : list[int], optional
        predictor harmonics to be plotted (default : [1, 2, 3, 4, 5])
    event_time : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    save_ext : str, optional
        plot extension (default : png)
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
    plt.savefig(os.path.join(save_path, plot_name + "." + save_ext), bbox_inches="tight", dpi=300)
    plt.close("all")
        
        
def plot_omegagram(pred, pred_name, target, target_name, gps, seconds, fs, plot_name, save_path,
                   norm=False, harmonics=None, event_time="center", save_ext="png", figsize=None):
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
        plot name
    save_path : str
        save path
    norm : bool, optional
        normalize predictor to 10 Hz (default : False)
    harmonics : list[int], optional
        predictor harmonics to be plotted (default : [1, 2, 3, 4, 5])
    event_time : str
        position of the event's gps in the analysed period.
        Can be `start`, `center`, or `end` (default : center)
    save_ext : str, optional
        plot extension (default : png)
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
    plt.savefig(os.path.join(save_path, plot_name + "." + save_ext), bbox_inches="tight", dpi=300)
    plt.close("all")
        

# def plot_imfs_summary(culprits, title, plot_name, save_path, dsort=True,
#                       batch=10, mean_freqs=None, save_ext="png"):
#     """Summary histogram of culprits found for a certain imf.
#
#     Parameters
#     ----------
#     culprits : list of str
#         channels names
#     title : str
#         plot title
#     plot_name : str
#         plot name
#     save_path : str
#         save path
#     dsort : bool, optional
#         descending order sort (default : True)
#     batch : int, optional
#         maximum number of points per plot (default : 10)
#     mean_freqs : dict, optional
#         dict {channel_name: [mean_freqs]} for mean frequencies histograms (default : None)
#     save_ext : str, optional
#         plot extension (default : png)
#     """
#     culprits = np.array(culprits)
#     seen = set()
#     uniq = [x for x in culprits if x not in seen and not seen.add(x)]
#     counts = []
#     for u in uniq:
#         counts.append(len(np.where(culprits == u)[0]))
#
#     if dsort:
#         uniq = [x for _, x in sorted(zip(counts, uniq), reverse=True)]
#         counts = sorted(counts, reverse=True)
#
#     if len(counts) <= batch:
#         plt.figure()
#         plt.bar(np.arange(len(counts)), counts, width=0.8)
#         plt.xticks(np.arange(len(counts)), uniq, rotation=45, horizontalalignment="right")
#         plt.ylabel("Channel occurrence")
#         plt.title(title)
#         plt.savefig(os.path.join(save_path, plot_name + "." + save_ext), bbox_inches="tight", dpi=300)
#         plt.close("all")
#     else:
#         n = len(counts) // batch
#         for i in range(n):
#             plt.figure()
#             plt.bar(np.arange(len(counts[i*batch:(i+1)*batch])), counts[i*batch:(i+1)*batch], width=0.8)
#             plt.xticks(np.arange(len(counts[i*batch:(i+1)*batch])), uniq[i*batch:(i+1)*batch],
#                        rotation=45, horizontalalignment="right")
#             plt.ylabel("Channel occurrence")
#             plt.title(title)
#             plt.savefig(os.path.join(save_path, plot_name + "_batch_" + str(i + 1) + "." + save_ext),
#                         bbox_inches="tight", dpi=300)
#             plt.close("all")
#         if len(counts) % batch != 0:
#             plt.figure()
#             plt.bar(np.arange(len(counts[n*batch:])), counts[n*batch:], width=0.8)
#             plt.xticks(np.arange(len(counts[n*batch:])), uniq[n*batch:],
#                        rotation=45, horizontalalignment="right")
#             plt.ylabel("Channel occurrence")
#             plt.title(title)
#             plt.savefig(os.path.join(save_path, plot_name + "_batch_" + str(n + 1) + "." + save_ext),
#                         bbox_inches="tight", dpi=300)
#             plt.close("all")
#
#     if mean_freqs is not None:
#         mf_x = []
#         mf_y = []
#         for i, key in enumerate(uniq[:batch]):
#             mf_x += [i for _ in range(len(mean_freqs[key]))]
#             mf_y += mean_freqs[key]
#         plot_mean_freq_summary(mf_x, mf_y, uniq, title, plot_name + "_mean_freq", save_path, save_ext=save_ext)


# def plot_corr_summary(gps_list, corr_list, title, plot_name, save_path, save_ext="png"):
#     """Summary of correlations for each GPS for a certain imf.
#
#     Parameters
#     ----------
#     gps_list : list
#         gps
#     corr_list : list
#         correlations
#     title : str
#         plot title
#     plot_name : str
#         plot name
#     save_path : str
#         save path
#     save_ext : str, optional
#         plot extension (default : png)
#     """
#     t1 = Time(gps_list[0], format="gps")
#     t2 = Time(t1, format="iso", scale="utc")
#     gps_date = "{:d} ({} UTC)\n".format(gps_list[0], t2)
#     x_values = [x - gps_list[0] for x in gps_list]
#
#     plt.figure()
#     plt.bar(x_values, corr_list, width=0.8)
#     plt.ylim(-1, 1)
#     plt.ylabel("$\\rho$")
#     plt.xlabel("t [s]\nfrom " + gps_date)
#     plt.title(title)
#     plt.savefig(os.path.join(save_path, plot_name + "." + save_ext), bbox_inches="tight", dpi=300)
#     plt.close("all")


# def plot_mean_freq_summary(x_vals, y_vals, x_labels, title, plot_name, save_path, save_ext="png"):
#     """Summary of mean frequencies for a channel and a certain imf.
#
#     Parameters
#     ----------
#     x_vals : list
#         list of integers for channel number
#     y_vals : list
#         mean_frequencies
#     x_labels : list[str]
#         channels names corresponding to `x_vals`
#     title : str
#         plot title
#     plot_name : str
#         plot name
#     save_path : str
#         save path
#     save_ext : str, optional
#         plot extension (default : str)
#     """
#     plt.figure()
#     plt.scatter(x_vals, y_vals)
#     plt.ylabel("Mean frequency [Hz]")
#     plt.xticks(np.arange(len(x_labels)), x_labels, rotation=45, horizontalalignment="right")
#     plt.title(title)
#     plt.savefig(os.path.join(save_path, plot_name + "." + save_ext), bbox_inches="tight", dpi=300)
#     plt.close("all")


def plot_imfs(folders, imfs, imf_thr=-1.0, combos=False, save_ext="png", figsize=None):
    """Plot imfs in the given folders.

    Parameters
    ----------
    folders : list[str]
        paths to the files needed for the plots
    imfs : str, list[int]
        imfs to plot, can be "all", or a list of integers
    imf_thr : float, optional
        correlation value above which to plot imfs (default : -1.0)
    combos : bool, optional
        if True, plot also combinations of imfs with the same channel name (default : False)
    save_ext : str, optional
        plots extension (default : png)
    figsize : tuple[int], optional
        figure size
    """
    imfs_to_plot = __plot_imfs_arg(imfs)
    if imfs_to_plot is None:
        return
    imf_thr = __corr_thr(imf_thr)

    for res_folder in folders:
        yf = file_utils.YmlFile(res_folder)
        ia = file_utils.load_imfs(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf.get_target_channel()
        gps_event = yf.get_gps()
        event_time = yf.get_event_position()
        fs = yf.get_sampling_frequency()

        if imfs_to_plot == "all":
            for n_imf in range(yf.get_imfs_count()):
                if yf.get_corr_of_imf(n_imf + 1) >= imf_thr:
                    plot_imf(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1), ia[:, n_imf], target_channel,
                             gps_event, fs, "$\\rho$ = {:.2f}".format(yf.get_corr_of_imf(n_imf + 1)),
                             "imf_{}_culprit".format(n_imf + 1), res_folder, event_time=event_time,
                             save_ext=save_ext, figsize=figsize)
            if combos:
                # ch_list = np.array([yf.get_channel_of_imf(el + 1) for el in range(yf.get_imfs_count())])
                plot_combinations(yf.get_combos(), ia, preds, target_channel, gps_event,
                                  fs, res_folder, thr=imf_thr, event_time=event_time,
                                  save_ext=save_ext, figsize=figsize)
        elif isinstance(imfs_to_plot, list):
            for n_imf in range(yf.get_imfs_count()):
                if n_imf + 1 in imfs_to_plot:
                    if yf.get_corr_of_imf(n_imf + 1) >= imf_thr:
                        plot_imf(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1), ia[:, n_imf],
                                 target_channel, gps_event, fs,
                                 "$\\rho$ = {:.2f}".format(yf.get_corr_of_imf(n_imf + 1)),
                                 "imf_{}_culprit".format(n_imf + 1), res_folder, event_time=event_time,
                                 save_ext=save_ext, figsize=figsize)
            if combos:
                #ch_list = np.array(
                #    [yf.get_channel_of_imf(el + 1) for el in range(yf.get_imfs_count()) if el + 1 in imfs_to_plot])
                # if len(ch_list) != 0:
                imfs_list = [imf for imf in imfs_to_plot if imf <= yf.get_imfs_count()]
                if len(imfs_list) > 0:
                    combos_dict = yf.get_combos()
                    for i in range(len(combos_dict[defines.IMF_KEY])):
                        is_ok = False
                        for imf_i in combos_dict[defines.IMF_KEY][i]:
                            if imf_i in imfs_list:
                                is_ok = True
                                break
                        if not is_ok:
                            _ = combos_dict[defines.IMF_KEY].pop(i)
                            _ = combos_dict[defines.CHANNEL_KEY].pop(i)
                            _ = combos_dict[defines.CORR_KEY].pop(i)
                    # plot_combinations(ch_list, ia[:, imfs_list], preds[:, imfs_list],
                    plot_combinations(combos_dict, ia, preds,
                                      target_channel, gps_event, fs, res_folder, thr=imf_thr,
                                      event_time=event_time, save_ext=save_ext, figsize=figsize)


def plot_omegagrams(folders, imfs, omegagram_thr=-1.0, harmonics=None,
                    save_ext="png", figsize=None):
    """Plot omegagrams in the given folders.

    Parameters
    ----------
    folders : list[str]
        paths to the files needed for the plots
    imfs : str, list[int]
        imfs to consider, can be "all", or a list of integers
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

    imfs_to_plot = __plot_imfs_arg(imfs)
    if imfs_to_plot is None:
        return
    omegagram_thr = __corr_thr(omegagram_thr)

    for res_folder in folders:
        yf = file_utils.YmlFile(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf.get_target_channel()
        gps = yf.get_gps()
        seconds = yf.get_seconds()
        event_time = yf.get_event_position()

        if imfs_to_plot == "all":
            for n_imf in range(yf.get_imfs_count()):
                if yf.get_corr_of_imf(n_imf + 1) >= omegagram_thr:
                    plot_omegagram_download(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1),
                                            target_channel, gps, seconds,
                                            "imf_{}_omegagram".format(n_imf + 1), res_folder,
                                            harmonics=harmonics, event_time=event_time,
                                            save_ext=save_ext, figsize=figsize)
        elif isinstance(imfs_to_plot, list):
            for n_imf in range(yf.get_imfs_count()):
                if n_imf + 1 in imfs_to_plot:
                    if yf.get_corr_of_imf(n_imf + 1) >= omegagram_thr:
                        plot_omegagram_download(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1),
                                                target_channel, gps, seconds,
                                                "imf_{}_omegagram".format(n_imf + 1), res_folder,
                                                harmonics=harmonics, event_time=event_time,
                                                save_ext=save_ext, figsize=figsize)
