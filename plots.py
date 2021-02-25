#  plots.py - this file is part of the asr package,
#  also known as "adaptive scattering recognizer".
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
import os
from .utils import file_utils
from .utils import plot_utils
from .common import defines


def plot_imfs_arg(arg):
    """Check on imfs to plot.
    
    Parameters
    ----------                                            
    arg : str, list
        imfs to plot arg
    """
    if arg in ["all", "max_corr"]:
        return arg
    elif isinstance(arg, list):
        return arg
    else:
        return None


def comparison_arg(arg):
    """Check on comparison.
    
    Parameters
    ----------                                            
    arg : str
        comparison arg
    """
    if arg == "max_corr":
        return arg
    elif isinstance(arg, list):
        return arg
    else:
        return None
        
        
def corr_thr(arg):
    """Check on correlation threshold.
    
    Parameters
    ----------                                            
    arg : float
        correlation threshold
    """
    arg = float(arg)
    if arg >= -1.0 and arg <= 1.0:
        return arg
    else:
        if arg < 0.0:
            return -1.0
        else:
            return 1.0


def plot_imfs(folders, imfs, imf_thr=-1.0, combos=False, save_ext="png"):
    """Plot imfs in the given folders.

    Parameters
    ----------
    folders : list[str]
        paths to the file needed for the plots
    imfs : str, list[int]
        imfs to plot, can be "all", "max_corr", or list of integers
    imf_thr : float
        correlation value above which to plot imfs (default : -1.0)
    combos : bool
        if True, plot also combinations of imfs with the same channel name (default : False)
    save_ext : str
        plots extension
    """
    imfs_to_plot = plot_imfs_arg(imfs)
    if imfs_to_plot is None:
        return
    imf_thr = corr_thr(imf_thr)

    for res_folder in folders:
        yf = file_utils.YmlFile(res_folder)
        ia = file_utils.load_imfs(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf.get_target_channel()
        gps_start, gps_end = yf.get_gps()
        fs = yf.get_sampling_frequency()
        gps_event = (gps_start + gps_end) // 2

        if imfs == "max_corr":
            n_imf = yf.get_max_corr_imf() - 1
            predictor_name = yf.get_max_corr_channel()
            if yf.get_max_corr() >= imf_thr:
                plot_utils.plot_imf(preds[:, n_imf], predictor_name, ia[:, n_imf], target_channel,
                                    gps_event, fs, "$\\rho$ = {:.4f}".format(yf.get_max_corr()),
                                    "max_corr_culprit", res_folder, save_ext=save_ext)
        elif imfs_to_plot == "all":
            for n_imf in range(yf.get_imfs_count()):
                if yf.get_corr_of_imf(n_imf + 1) >= imf_thr:
                    plot_utils.plot_imf(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1), ia[:, n_imf], target_channel,
                                        gps_event, fs, "$\\rho$ = {:.4f}".format(yf.get_corr_of_imf(n_imf + 1)),
                                        "imf_{}_culprit".format(n_imf + 1), res_folder, save_ext=save_ext)
            if combos:
                ch_list = np.array([yf.get_channel_of_imf(el + 1) for el in range(yf.get_imfs_count())])
                plot_utils.plot_combinations(ch_list, ia, preds, target_channel, gps_event,
                                             fs, res_folder, thr=imf_thr, save_ext=save_ext)
        elif isinstance(imfs_to_plot, list):
            for n_imf in range(yf.get_imfs_count()):
                if n_imf + 1 in imfs_to_plot:
                    if yf.get_corr_of_imf(n_imf + 1) >= imf_thr:
                        plot_utils.plot_imf(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1), ia[:, n_imf],
                                            target_channel, gps_event, fs,
                                            "$\\rho$ = {:.4f}".format(yf.get_corr_of_imf(n_imf + 1)),
                                            "imf_{}_culprit".format(n_imf + 1), res_folder, save_ext=save_ext)
            if combos:
                ch_list = np.array(
                    [yf.get_channel_of_imf(el + 1) for el in range(yf.get_imfs_count()) if el + 1 in imfs_to_plot])
                if len(ch_list) != 0:
                    imfs_list = [imf - 1 for imf in imfs_to_plot if imf <= yf.get_imfs_count()]
                    plot_utils.plot_combinations(ch_list, ia[:, imfs_list], preds[:, imfs_list],
                                                 target_channel, gps_event, fs, res_folder, thr=imf_thr,
                                                 save_ext=save_ext)


def plot_omegagrams(folders, imfs, omegagram_thr=-1.0, harmonics=[1, 2, 3, 4, 5], save_ext="png"):
    """Plot omegagrams in the given folders.

    Parameters
    ----------
    folders : list[str]
        paths to the file needed for the plots
    imfs : str, list[int]
        imfs to consider, can be "all", "max_corr", or list of integers
    omegagram_thr : float
        correlation value above which to plot omegagrams (default : -1.0)
    harmonics : list[int]
        harmonics for the culprit (default : [1, 2, 3, 4, 5])
    save_ext : str
        plots extension
    """
    imfs_to_plot = plot_imfs_arg(imfs)
    if imfs_to_plot is None:
        return
    omegagram_thr = corr_thr(omegagram_thr)

    for res_folder in folders:
        yf = file_utils.YmlFile(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf.get_target_channel()
        gps_start, gps_end = yf.get_gps()

        if imfs == "max_corr":
            n_imf = yf.get_max_corr_imf() - 1
            predictor_name = yf.get_max_corr_channel()
            if yf.get_max_corr() >= omegagram_thr:
                plot_utils.plot_omegagram_download(preds[:, n_imf], predictor_name, target_channel, gps_start,
                                                   gps_end, "max_corr_omegagram", res_folder, harmonics=harmonics,
                                                   save_ext=save_ext)
        elif imfs_to_plot == "all":
            for n_imf in range(yf.get_imfs_count()):
                if yf.get_corr_of_imf(n_imf + 1) >= omegagram_thr:
                    plot_utils.plot_omegagram_download(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1),
                                                       target_channel, gps_start, gps_end,
                                                       "imf_{}_omegagram".format(n_imf + 1), res_folder,
                                                       harmonics=harmonics, save_ext=save_ext)
        elif isinstance(imfs_to_plot, list):
            for n_imf in range(yf.get_imfs_count()):
                if n_imf + 1 in imfs_to_plot:
                    if yf.get_corr_of_imf(n_imf + 1) >= omegagram_thr:
                        plot_utils.plot_omegagram_download(preds[:, n_imf], yf.get_channel_of_imf(n_imf + 1),
                                                           target_channel, gps_start, gps_end,
                                                           "imf_{}_omegagram".format(n_imf + 1), res_folder,
                                                           harmonics=harmonics, save_ext=save_ext)


def plot_comparisons(folders, comparison, comparison_thr=-1.0, save_ext="png"):
    """Comparison plots of results.

    Parameters
    ----------
    folders : list[str]
        paths to the file needed for the plots
    comparison : str, list[int]
        imfs for comparison, can be "max_corr", or list
    comparison_thr : float
        correlation threshold for comparison (default : -1.0)
    save_ext : str
        plots extension
    """
    ################################################################
    ## to be lightened, create a table with pandas and then do plots
    ## optionally save table to a csv
    ################################################################
    comparison = plot_imfs_arg(comparison)
    if comparison is None:
        return
    comparison_thr = corr_thr(comparison_thr)

    if len(folders) == 0:
        return

    folders_path = os.path.sep.join(folders[0].split(os.path.sep)[:-1])
    cpath = os.path.join(folders_path, "comparison")
    if not os.path.isdir(cpath):
        os.makedirs(cpath, exist_ok=True)

    target_channel = file_utils.YmlFile(folders[0]).get_target_channel()

    if comparison == "max_corr":
        culprits_list = []
        gps_list = []
        corr_list = []
        m_freq_list = []
        f = open(os.path.join(cpath, "max_corr_summary_table_{}.csv".format("_".join(target_channel.split(":")))), "w")
        f.write("start,end,channel,corr,mean_freq\n")
        for res_folder in folders:
            yf = file_utils.YmlFile(res_folder)
            if yf.get_max_corr() >= comparison_thr:
                culprits_list.append(yf.get_max_corr_channel())
                gps_start, gps_end = yf.get_gps()
                gps_list.append((gps_start + gps_end) // 2)
                corr_list.append(yf.get_max_corr())
                m_freq_list.append(yf.get_max_corr_mean_freq())
                f.write("{:d},{:d},{},{:f},{:f}\n".format(gps_start, gps_end, yf.get_max_corr_channel(),
                                                          yf.get_max_corr(), yf.get_max_corr_mean_freq()))
        f.close()
        plot_utils.plot_imfs_summary(culprits_list, "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                     "max_corr_imf_summary_{}".format("_".join(target_channel.split(":"))),
                                     cpath, save_ext=save_ext)
        plot_utils.plot_corr_summary(gps_list, corr_list, "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                     "max_corr_summary_{}".format("_".join(target_channel.split(":"))),
                                     cpath, save_ext=save_ext)
    elif isinstance(comparison, list):
        culprits_dict = {}
        corr_dict = {}
        f_dict = {}
        for i in comparison:
            culprits_dict[i] = []
            corr_dict[i] = {defines.GPS_KEY: [], defines.CORR_KEY: []}
            f_dict[i] = ""

        for res_folder in folders:
            yf = file_utils.YmlFile(res_folder)
            for n_imf in comparison:
                gps_start, gps_end = yf.get_gps()
                gps_mid_time = (gps_start + gps_end) // 2
                corr_dict[n_imf][defines.GPS_KEY].append(gps_mid_time)
                if yf.get_imfs_count() >= n_imf:
                    if yf.get_corr_of_imf(n_imf) >= comparison_thr:
                        culprits_dict[n_imf].append(yf.get_channel_of_imf(n_imf))
                        corr_dict[n_imf][defines.CORR_KEY].append(yf.get_corr_of_imf(n_imf))
                        f_dict[n_imf] += "{:d},{:d},{},{:f},{:f}\n".format(gps_start, gps_end,
                                                                           yf.get_channel_of_imf(n_imf),
                                                                           yf.get_corr_of_imf(n_imf),
                                                                           yf.get_mean_freq_of_imf(n_imf))
                    else:
                        corr_dict[n_imf][defines.CORR_KEY].append(0.0)
                else:
                    culprits_dict[n_imf].append("Not found")
                    corr_dict[n_imf][defines.CORR_KEY].append(0.0)

        mean_freq_dict = {}
        for n_imf in comparison:
            mf_dict = {}
            for chnl in culprits_dict[n_imf]:
                if culprits_dict[n_imf] != "Not found":
                    mf_dict[chnl] = []
            for res_folder in folders:
                yf = file_utils.YmlFile(res_folder)
                if yf.get_imfs_count() >= n_imf:
                    if yf.get_corr_of_imf(n_imf) >= comparison_thr:
                        mf_dict[yf.get_channel_of_imf(n_imf)].append(yf.get_mean_freq_of_imf(n_imf))
            mean_freq_dict[n_imf] = mf_dict

        for n_imf in comparison:
            if f_dict[n_imf] != "":
                f = open(os.path.join(cpath, "imf_{}_summary_table_{}.csv".format(n_imf, "_".join(target_channel.split(":")))), "w")
                f.write("start,end,channel,corr,mean_freq\n")
                f.write(f_dict[n_imf])
                f.close()

            if len(culprits_dict[n_imf]) > 0 and culprits_dict[n_imf].count("Not found") != len(culprits_dict[n_imf]):
                plot_utils.plot_imfs_summary(culprits_dict[n_imf],
                                             "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                             "imf_{:d}_summary_{}".format(n_imf, "_".join(target_channel.split(":"))),
                                             cpath, mean_freqs=mean_freq_dict[n_imf], save_ext=save_ext)

            corr_check = np.sum(np.abs(corr_dict[n_imf][defines.CORR_KEY]))
            if corr_check != 0.0:
                zero_idxs = np.where(np.asarray(corr_dict[n_imf][defines.CORR_KEY]) != 0.0)[0]
                plot_utils.plot_corr_summary(np.asarray(corr_dict[n_imf][defines.GPS_KEY])[zero_idxs],
                                             np.asarray(corr_dict[n_imf][defines.CORR_KEY])[zero_idxs],
                                             "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                             "imf_{:d}_corr_summary_{}".format(n_imf,
                                                                               "_".join(target_channel.split(":"))),
                                             cpath, save_ext=save_ext)

