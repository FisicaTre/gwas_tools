#  asr_plots.py - this file is part of the asr package,
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
import argparse
import os
from .utils import file_utils
from .utils import plot_utils
from .common import defines


def plot_imfs_arg(arg):
    """Check on imfs to plot.
    
    Parameters
    ----------                                            
    arg : str
        imfs to plot arg
    """
    if arg in ["all", "max_corr"]:
        return arg
    else:
        try:
            return [int(s) for s in arg.split(",")]
        except:
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
    else:
        try:
            return [int(s) for s in arg.split(",")]
        except:
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


def plots(ipath, single_folder=False, imfs_to_plot=None, imf_thr=None,
          omegagram_thr=None, comparison=None, comparison_thr=None):
    """Plot analysis results.
    
    Parameters
    ----------
    ipath : str
        path to results folder
    single_folder : bool
        whether or not `ipath` points to a single folder or a group of folders
    imfs_to_plot : str
        imfs to plot, can be "all", "max_corr", or comma-separated numbers
    imf_thr : float
        imfs threshold
    omegagram_thr : float
        omegagram threshold
    comparison : str
        imfs fo comparison, can be "all", "max_corr", or comma-separated numbers
    comparison_thr : float
        comparison imfs threshold
    """
    imfs_to_plot = plot_imfs_arg(imfs_to_plot)
    if imf_thr is not None:
        imf_thr = corr_thr(imf_thr)
    else:
        imf_thr = -1.0
    if omegagram_thr is not None:
        omegagram_thr = corr_thr(omegagram_thr)
    comparison = comparison_arg(comparison)
    if comparison_thr is not None:
        comparison_thr = corr_thr(comparison_thr)
    else:
        comparison_thr = -1.0
    
    res_folders = []
    if single_folder:
        comparison = None
        comparison_threshold = None
        if file_utils.yml_exists(ipath) and file_utils.imfs_exists(ipath) and file_utils.predictors_exists(ipath):
            res_folders.append(ipath)
    else:
        for folder in os.listdir(ipath):
            curr_dir = os.path.join(ipath, folder)
            if os.path.isdir(curr_dir):
                if file_utils.yml_exists(curr_dir) and file_utils.imfs_exists(curr_dir) and file_utils.predictors_exists(curr_dir):
                    res_folders.append(curr_dir)

    for res_folder in res_folders:
        yf = file_utils.load_yml(res_folder)
        ia = file_utils.load_imfs(res_folder)
        preds = file_utils.load_predictors(res_folder)

        target_channel = yf[defines.PARAMS_SECT_KEY][defines.TARGET_CH_KEY]
        gps = yf[defines.PARAMS_SECT_KEY][defines.GPS_KEY]
        fs = yf[defines.PARAMS_SECT_KEY][defines.SAMP_FREQ_KEY]
        gps_start = int(gps.split(",")[0]) + defines.EXTRA_SECONDS
        gps_end = int(gps.split(",")[1]) - defines.EXTRA_SECONDS
        gps_event = (gps_start + gps_end) // 2

        if imfs_to_plot is not None:
            if imfs_to_plot == "max_corr":
                n_imf = yf[defines.MAX_CORR_SECT_KEY][defines.IMF_KEY] - 1
                predictor_name = yf[defines.MAX_CORR_SECT_KEY][defines.CHANNEL_KEY]
                if yf[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY] >= imf_thr:
                    plot_utils.plot_imf(preds[:, n_imf], predictor_name, ia[:, n_imf], target_channel,
                                        gps_event, fs, "$\\rho$ = {:.4f}".format(yf[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY]),
                                        "max_corr_culprit", res_folder)
                if omegagram_thr is not None and yf[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY] >= omegagram_thr:
                    plot_utils.plot_omegagram_download(preds[:, n_imf], target_channel, gps_start, gps_end,
                                                       "max_corr_omegagram", res_folder)
            elif imfs_to_plot == "all":
                for n_imf, imf_pred in enumerate(yf[defines.CORR_SECT_KEY]):
                    if imf_pred[defines.CORR_KEY] >= imf_thr:
                        plot_utils.plot_imf(preds[:, n_imf], imf_pred[defines.CHANNEL_KEY], ia[:, n_imf], target_channel,
                                            gps_event, fs, "$\\rho$ = {:.4f}".format(imf_pred[defines.CORR_KEY]),
                                            "imf_{}_culprit".format(imf_pred[defines.IMF_KEY]), res_folder)
                    if omegagram_thr is not None and imf_pred[defines.CORR_KEY] >= omegagram_thr:
                        plot_utils.plot_omegagram_download(preds[:, n_imf], target_channel, gps_start, gps_end,
                                                           "imf_{}_omegagram".format(imf_pred[defines.IMF_KEY]), res_folder)

                ch_list = np.array([el[defines.CHANNEL_KEY] for el in yf[defines.CORR_SECT_KEY]])
                plot_utils.plot_combinations(ch_list, ia, preds, target_channel, gps_event,
                                             fs, res_folder, imf_thr)
            elif isinstance(imfs_to_plot, list):
                for n_imf, imf_pred in enumerate(yf[defines.CORR_SECT_KEY]):
                    if imf_pred[defines.IMF_KEY] in imfs_to_plot:
                        if imf_pred[defines.CORR_KEY] >= imf_thr:
                            plot_utils.plot_imf(preds[:, n_imf], imf_pred[defines.CHANNEL_KEY], ia[:, n_imf], target_channel,
                                                gps_event, fs, "$\\rho$ = {:.4f}".format(imf_pred[defines.CORR_KEY]),
                                                "imf_{}_culprit".format(imf_pred[defines.IMF_KEY]), res_folder)
                        if omegagram_thr is not None and imf_pred["corr"] >= omegagram_thr:
                            plot_utils.plot_omegagram_download(preds[:, n_imf], target_channel, gps_start, gps_end,
                                                                "imf_{}_omegagram".format(imf_pred[defines.IMF_KEY]), res_folder)

                ch_list = np.array([el[defines.CHANNEL_KEY] for el in yf[defines.CORR_SECT_KEY] if el[defines.IMF_KEY] in imfs_to_plot])
                if len(ch_list) != 0:
                    imfs_list = [imf - 1 for imf in imfs_to_plot if imf <= len(yf[defines.CORR_SECT_KEY])]
                    plot_utils.plot_combinations(ch_list, ia[:, imfs_list], preds[:, imfs_list],
                                                 target_channel, gps_event, fs, res_folder, imf_thr)

    if comparison is not None:
        cpath = os.path.join(ipath, "comparison")
        if not os.path.isdir(cpath):
            os.makedirs(cpath, exist_ok=True)
            
        if comparison == "max_corr":
            culprits_list = []
            gps_list = []
            corr_list = []
            f = open(os.path.join(cpath, "max_corr_summary_list_{}.csv".format(target_channel)), "w")
            f.write("start,end,corr,channel\n")
            res_folders.sort()
            for res_folder in res_folders:
                yf = file_utils.load_yml(res_folder)
                if yf[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY] >= comparison_thr:
                    culprits_list.append(yf[defines.MAX_CORR_SECT_KEY][defines.CHANNEL_KEY])
                    gps_times = yf[defines.PARAMS_SECT_KEY][defines.GPS_KEY].split(",")
                    gps_list.append(str((int(gps_times[0]) + int(gps_times[1])) // 2))
                    corr_list.append(yf[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY])
                    f.write("{:d},{:d},{},{}\n".format(int(yf[defines.PARAMS_SECT_KEY][defines.GPS_KEY].split(",")[0]) + defines.EXTRA_SECONDS,
                                                       int(yf[defines.PARAMS_SECT_KEY][defines.GPS_KEY].split(",")[1]) - defines.EXTRA_SECONDS,
                                                       yf[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY],
                                                       yf[defines.MAX_CORR_SECT_KEY][defines.CHANNEL_KEY]))
              
            f.close()
            plot_utils.plot_imfs_summary(culprits_list, "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                         "max_corr_imf_summary_{}".format(target_channel), cpath)
            plot_utils.plot_corr_summary(gps_list, corr_list, "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                         "max_corr_summary_{}".format(target_channel), cpath, batch=20)
        elif isinstance(comparison, list):
            culprits_dict = {}
            corr_dict = {}
            f_dict = {}
            for i in comparison:
                culprits_dict[i] = []
                corr_dict[i] = {defines.GPS_KEY: [], defines.CORR_KEY: []}
                f_dict[i] = ""
                
            res_folders.sort()
            for res_folder in res_folders:
                yf = file_utils.load_yml(res_folder)
                for n_imf in comparison:
                    gps_times = yf[defines.PARAMS_SECT_KEY][defines.GPS_KEY].split(",")
                    gps_mid_time = str((int(gps_times[0]) + int(gps_times[1])) // 2)
                    corr_dict[n_imf][defines.GPS_KEY].append(gps_mid_time)
                    if len(yf[defines.CORR_SECT_KEY]) >= n_imf:
                        if yf[defines.CORR_SECT_KEY][n_imf - 1][defines.CORR_KEY] >= comparison_thr:
                            culprits_dict[n_imf].append(yf[defines.CORR_SECT_KEY][n_imf - 1][defines.CHANNEL_KEY])
                            corr_dict[n_imf][defines.CORR_KEY].append(yf[defines.CORR_SECT_KEY][n_imf - 1][defines.CORR_KEY])
                            f_dict[n_imf] += "{:d},{:d},{},{}\n".format(int(yf[defines.PARAMS_SECT_KEY][defines.GPS_KEY].split(",")[0]) + defines.EXTRA_SECONDS,
                                                                       int(yf[defines.PARAMS_SECT_KEY][defines.GPS_KEY].split(",")[1]) - defines.EXTRA_SECONDS,
                                                                       yf[defines.CORR_SECT_KEY][n_imf - 1][defines.CORR_KEY],
                                                                       yf[defines.CORR_SECT_KEY][n_imf - 1][defines.CHANNEL_KEY])
                        else:
                            corr_dict[n_imf][defines.CORR_KEY].append(0.0)
                    else:
                        culprits_dict[n_imf].append("Not found")
                        corr_dict[n_imf][defines.CORR_KEY].append(0.0)
            
            for n_imf in comparison:
                if f_dict[n_imf] != "":
                    f = open(os.path.join(cpath, "imf_{}_summary_list_{}.csv".format(n_imf, target_channel)), "w")
                    f.write("start,end,corr,channel\n")
                    f.write(f_dict[n_imf])
                    f.close()
                
                if len(culprits_dict[n_imf]) > 0 and culprits_dict[n_imf].count("Not found") != len(culprits_dict[n_imf]):
                    plot_utils.plot_imfs_summary(culprits_dict[n_imf], "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                                 "imf_{:d}_summary_{}".format(n_imf, target_channel), cpath)
                    
                corr_check = np.sum(np.abs(corr_dict[n_imf][defines.CORR_KEY]))
                if corr_check != 0.0:
                    zero_idxs = np.where(np.asarray(corr_dict[n_imf][defines.CORR_KEY]) != 0.0)[0]
                    plot_utils.plot_corr_summary(np.asarray(corr_dict[n_imf][defines.GPS_KEY])[zero_idxs],
                                                 np.asarray(corr_dict[n_imf][defines.CORR_KEY])[zero_idxs],
                                                 "{} (thr {:.2f})".format(target_channel, comparison_thr),
                                                 "imf_{:d}_corr_summary_{}".format(n_imf, target_channel),
                                                 cpath, batch=10)
