#  file_utils.py - this file is part of the gwadaptive_scattering package.
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


import scipy.io
import yaml
import os
import glob
import pickle
import numpy as np
import pandas as pd
from ..common import defines


class YmlFile(object):
    """Class for the yml output file.

    Parameters
    ----------
    yml_file : str, optional
        path to the yml file (default : None)
    """

    def __init__(self, yml_file=None):
        self.name = "output.yml"
        if yml_file is None:
            self.content = {}
        else:
            with open(os.path.join(yml_file, self.name), "r") as f:
                self.content = yaml.safe_load(f)

    def get_target_channel(self):
        """Target channel.

        Returns
        -------
        str
            target channel
        """
        try:
            return self.content[defines.PARAMS_SECT_KEY][defines.TARGET_CH_KEY]
        except:
            raise ValueError("Value not found.")

    def get_gps(self):
        """GPS time.

        Returns
        -------
        int
            GPS time
        """
        try:
            gps = self.content[defines.PARAMS_SECT_KEY][defines.GPS_KEY]
            return gps
        except:
            raise ValueError("Value not found.")

    def get_seconds(self):
        """Duration of the analyzed period.

        Returns
        -------
        int
            duration of the analyzed period
        """
        try:
            seconds = self.content[defines.PARAMS_SECT_KEY][defines.SECONDS_KEY]
            return seconds
        except:
            raise ValueError("Value not found.")

    def get_event_position(self):
        """Position of the event in the analyzed period.

        Returns
        -------
        str
            position of the event in the analyzed period
        """
        try:
            event_pos = self.content[defines.PARAMS_SECT_KEY][defines.EVENT_KEY]
            return event_pos
        except:
            raise ValueError("Value not found.")

    def get_sampling_frequency(self):
        """Sampling frequency.

        Returns
        -------
        float
            sampling frequency
        """
        try:
            return self.content[defines.PARAMS_SECT_KEY][defines.SAMP_FREQ_KEY]
        except:
            raise ValueError("Value not found.")

    def get_channels_list(self):
        """Channels list.

        Returns
        -------
        str
            path to channels list file
        """
        try:
            return self.content[defines.PARAMS_SECT_KEY][defines.CH_LIST_KEY]
        except:
            raise ValueError("Value not found.")

    def get_output_path(self):
        """Output path.

        Returns
        -------
        str
            path to where output analysis was saved
        """
        try:
            return self.content[defines.PARAMS_SECT_KEY][defines.OUT_PATH_KEY]
        except:
            raise ValueError("Value not found.")

    def get_lowpass_frequency(self):
        """Lowpass frequency.

        Returns
        -------
        float
            lowpass frequency
        """
        try:
            return self.content[defines.PARAMS_SECT_KEY][defines.LOWPASS_FREQ_KEY]
        except:
            raise ValueError("Value not found.")

    def get_scattering_factor(self):
        """Scattering factor.

        Returns
        -------
        int
            scattering factor
        """
        try:
            return self.content[defines.PARAMS_SECT_KEY][defines.SCATTERING_KEY]
        except:
            raise ValueError("Value not found.")

    def get_smoothing_window(self):
        """Smoothing window.

        Returns
        -------
        int
            smoothing window
        """
        try:
            return self.content[defines.PARAMS_SECT_KEY][defines.SMOOTH_WIN_KEY]
        except:
            raise ValueError("Value not found.")

    def is_locked(self):
        """Lock info.

        Returns
        -------
        bool
            lock info
        """
        try:
            lck = self.content[defines.LOCK_INFO_KEY]
            if lck == "yes":
                return True
            elif lck == "no":
                return False
            else:
                return None
        except:
            return None

    def get_imfs_count(self):
        """Get number of found imfs.

        Returns
        -------
        int
            number of found imfs
        """
        if defines.CORR_SECT_KEY in self.content.keys():
            return len(self.content[defines.CORR_SECT_KEY])
        else:
            raise ValueError("{} section not found.".format(defines.CORR_SECT_KEY))

    def get_channel_of_imf(self, imf, second_best=False):
        """Most correlated channel with imf `imf`.

        Parameters
        ----------
        imf : int
            imf number
        second_best : bool
            if True, returns the value from the section for
            the second best culprit

        Returns
        -------
        str
            most correlated channel with imf `imf`
        """
        try:
            if second_best:
                return self.content[defines.CORR_2_SECT_KEY][imf - 1][defines.CHANNEL_KEY]
            else:
                return self.content[defines.CORR_SECT_KEY][imf - 1][defines.CHANNEL_KEY]
        except:
            raise ValueError("Value not found.")

    def get_corr_of_imf(self, imf, second_best=False):
        """Correlation of imf `imf`.

        Parameters
        ----------
        imf : int
            imf number
        second_best : bool
            if True, returns the value from the section for
            the second best culprit

        Returns
        -------
        float
            correlation of imf `imf`
        """
        try:
            if second_best:
                return self.content[defines.CORR_2_SECT_KEY][imf - 1][defines.CORR_KEY]
            else:
                return self.content[defines.CORR_SECT_KEY][imf - 1][defines.CORR_KEY]
        except:
            raise ValueError("Value not found.")

    def get_mean_freq_of_imf(self, imf, second_best=False):
        """Mean frequency of imf `imf`.

        Parameters
        ----------
        imf : int
            imf number
        second_best : bool
            if True, returns the value from the section for
            the second best culprit

        Returns
        -------
        float
            mean frequency of imf `imf`
        """
        try:
            if second_best:
                return self.content[defines.CORR_2_SECT_KEY][imf - 1][defines.MEAN_FREQ_KEY]
            else:
                return self.content[defines.CORR_SECT_KEY][imf - 1][defines.MEAN_FREQ_KEY]
        except:
            raise ValueError("Value not found.")

    def get_corr_of_combo_with_imf(self, imf):
        """Correlation of combo with imf `imf`.

        Parameters
        ----------
        imf : int
            imf number to be present in combo

        Returns
        -------
        float
            correlation of combo with imf `imf`
        """
        combos = self.content[defines.COMBO_SECT_KEY]
        for i in range(len(combos)):
            imfs = [int(imf) for imf in combos[i][defines.IMF_KEY].split(",")]
            if imf in imfs:
                return combos[i][defines.CORR_KEY]

        raise ValueError("Value not found.")

    def get_channel_of_combo_with_imf(self, imf):
        """Most correlated channel with combo with imf `imf`.

        Parameters
        ----------
        imf : int
            imf number to be present in combo

        Returns
        -------
        float
            correlation of combo with imf `imf`
        """
        combos = self.content[defines.COMBO_SECT_KEY]
        for i in range(len(combos)):
            imfs = [int(imf) for imf in combos[i][defines.IMF_KEY].split(",")]
            if imf in imfs:
                return combos[i][defines.CHANNEL_KEY]

        raise ValueError("Value not found.")

    def write_parameters(self, gps, seconds, event, target_channel_name, channels_file, out_path,
                         fs, f_lowpass, n_scattering, smooth_win):
        """Write parameters section to file.

        Parameters
        ----------
        gps : int
            gps time
        seconds : int
            seconds analyzed
        event : str
            position of `gps` in the analyzed period (`start`, `center`, or `end`)
        target_channel_name : str
            target channel name
        channels_file : str
            channels list file
        out_path : str
            output path for results
        fs : float
            sampling frequency
        f_lowpass : float
            lowpass frequency
        n_scattering : int
            scattered light bounces
        smooth_win : int
            smoothing window
        """
        self.content[defines.PARAMS_SECT_KEY] = {}
        self.content[defines.PARAMS_SECT_KEY][defines.GPS_KEY] = int(gps)
        self.content[defines.PARAMS_SECT_KEY][defines.SECONDS_KEY] = int(seconds)
        self.content[defines.PARAMS_SECT_KEY][defines.EVENT_KEY] = event
        self.content[defines.PARAMS_SECT_KEY][defines.TARGET_CH_KEY] = target_channel_name
        self.content[defines.PARAMS_SECT_KEY][defines.CH_LIST_KEY] = channels_file
        self.content[defines.PARAMS_SECT_KEY][defines.OUT_PATH_KEY] = out_path
        self.content[defines.PARAMS_SECT_KEY][defines.SAMP_FREQ_KEY] = float(fs)
        self.content[defines.PARAMS_SECT_KEY][defines.LOWPASS_FREQ_KEY] = float(f_lowpass)
        self.content[defines.PARAMS_SECT_KEY][defines.SCATTERING_KEY] = int(n_scattering)
        self.content[defines.PARAMS_SECT_KEY][defines.SMOOTH_WIN_KEY] = int(smooth_win)

    def write_lock_info(self, locked):
        """Write lock info to file.

        Parameters
        ----------
        locked : bool
            if interferometers is locked
        """
        self.content[defines.LOCK_INFO_KEY] = "yes" if locked else "no"

    def write_correlation_section(self, channels, corrs, mean_freqs):
        """Write correlation section to file.

        Parameters
        ----------
        channels : list[str]
            list of culprits for each imf
        corrs : list[float]
            list of correlations for each imf
        mean_freqs : list[float]
            list of mean frequencies for each imf
        """
        self.content[defines.CORR_SECT_KEY] = []
        for i in range(len(channels)):
            tmp_dict = {defines.IMF_KEY: int(i + 1), defines.CHANNEL_KEY: channels[i],
                        defines.CORR_KEY: float(corrs[i]), defines.MEAN_FREQ_KEY: float(mean_freqs[i])}
            self.content[defines.CORR_SECT_KEY].append(tmp_dict)

    def write_combo_section(self, imfs, channels, corrs):
        """Write combo section to file.

        Parameters
        ----------
        imfs : list[list[int]]
            list of imfs lists belonging to the same channel
        channels : list[str]
            list of culprits for each combo
        corrs : list[float]
            list of correlations for each combo
        """
        self.content[defines.COMBO_SECT_KEY] = []
        for i in range(len(channels)):
            tmp_dict = {defines.IMF_KEY: ",".join([str(imf) for imf in imfs[i]]), defines.CHANNEL_KEY: channels[i],
                        defines.CORR_KEY: float(corrs[i])}
            self.content[defines.COMBO_SECT_KEY].append(tmp_dict)

    def get_combos(self):
        """Get combos section.
        """
        combos = self.content[defines.COMBO_SECT_KEY]
        for i in range(len(combos)):
            imfs = combos[i][defines.IMF_KEY]
            combos[i][defines.IMF_KEY] = [int(imf) for imf in imfs.split(",")]
        return combos

    def write_seismic_channels(self, seismic_dict):
        """Write seismic channels values.

        Parameters
        ----------
        seismic_dict : dict
            seismic channels values (key : channel name, value : channel value)
        """
        self.content[defines.SEISMIC_SECT_KEY] = {}
        for k in seismic_dict.keys():
            self.content[defines.SEISMIC_SECT_KEY][k] = float(seismic_dict[k])

    def get_seismic_channels(self):
        """Get seismic channels section.
        """
        return self.content[defines.SEISMIC_SECT_KEY]

    def write_2nd_best_correlation_section(self, channels, corrs, mean_freqs):
        """Write correlation section corresponding to the second
        best culprit of each imf to file.

        Parameters
        ----------
        channels : list[str]
            list of culprits for each imf
        corrs : list[float]
            list of correlations for each imf
        mean_freqs : list[float]
            list of mean frequencies for each imf
        """
        self.content[defines.CORR_2_SECT_KEY] = []
        for i in range(len(channels)):
            tmp_dict = {defines.IMF_KEY: int(i + 1), defines.CHANNEL_KEY: channels[i],
                        defines.CORR_KEY: float(corrs[i]), defines.MEAN_FREQ_KEY: float(mean_freqs[i])}
            self.content[defines.CORR_2_SECT_KEY].append(tmp_dict)

    def save(self, save_path):
        """Save file.

        Parameters
        ----------
        save_path : str
            path to the output file
        """
        with open(os.path.join(save_path, self.name), "w") as yaml_file:
            yaml.dump(self.content, yaml_file, default_flow_style=False, sort_keys=False)


def from_mat(mat_file, mat_col):
    """Get array from .mat file.
    
    Parameters
    ----------
    mat_file : str
        path to the .mat file
    mat_col : int
        column of the array to get
        
    Returns
    -------
    numpy array
        array from .mat
    """
    mat = scipy.io.loadmat(mat_file)
    mat_data = mat[mat_col]

    return mat_data


def yml_exists(yml_path):
    """Check if yml file exists.
    
    Parameters
    ----------
    yml_path : str
        path to yml file
        
    Returns
    -------
    bool
        yml existence
    """
    yml_file = os.path.join(yml_path, "output.yml")
    if not os.path.exists(yml_file):
        print("File {} not found.".format(yml_file))
        return False

    return True


def imfs_exists(imfs_path):
    """Check if imfs file exists.
    
    Parameters
    ----------
    imfs_path : str
        path to imfs file
        
    Returns
    -------
    bool
        imfs existence
    """
    imfs_file = glob.glob(os.path.join(imfs_path, "*.imf"))
    if len(imfs_file) != 1:
        print("File {} not found or more than one `.imf` found.".format(imfs_file))
        return False

    return True


def predictors_exists(predictors_path):
    """Check if predictors file exists.
    
    Parameters
    ----------
    predictors_path : str
        path to predictors file
        
    Returns
    -------
    bool
        predictors existence
    """
    predictors_file = glob.glob(os.path.join(predictors_path, "*.predictors"))
    if len(predictors_file) != 1:
        print("File {} not found or more than one `.predictors` found.".format(predictors_file))
        return False

    return True


def load_imfs(imfs_path):
    """Load imfs file.
    
    Parameters
    ----------
    imfs_path : str
        path to imfs file
        
    Returns
    -------
    numpy ndarray
        imfs file
    """
    imfs_file = glob.glob(os.path.join(imfs_path, "*.imf"))
    f = open(imfs_file[0], "rb")
    imfs = pickle.load(f)
    f.close()

    return imfs


def load_predictors(predictors_path):
    """Load predictors file.
    
    Parameters 
    ----------
    predictors_path : str
        path to predictors file
        
    Returns
    -------
    numpy ndarray
        predictors file
    """
    predictors_file = glob.glob(os.path.join(predictors_path, "*.predictors"))
    f = open(predictors_file[0], "rb")
    predictors = pickle.load(f)
    f.close()

    return predictors


def get_results_folders(results_path, sort=True, must_include=None, filter_non_valid=True):
    """Get list of folders with results from one or multiple analyses.

    Parameters
    ----------
    results_path : str
        path to results folders
    sort : bool, optional
        sort folders (default : True)
    must_include : list[str], optional
        patterns of files that must be present in a
        folder, otherwise it is discarded. Ignored if
        `filter_non_valid` is False (default : None)
    filter_non_valid : bool, optional
        whether or not exclude non valid folders (default : True)

    Returns
    -------
    list[str]
        list of folders paths
    """
    res_folders = []
    for folder in os.listdir(results_path):
        curr_dir = os.path.join(results_path, folder)
        if os.path.isdir(curr_dir):
            if filter_non_valid:
                if is_valid_folder(curr_dir):
                    if must_include is not None:
                        add_folder = True
                        for pattern in must_include:
                            if len(glob.glob(os.path.join(curr_dir, pattern))) == 0:
                                add_folder = False
                                break
                        if add_folder:
                            res_folders.append(curr_dir)
                    else:
                        res_folders.append(curr_dir)
            else:
                res_folders.append(curr_dir)

    if sort:
        res_folders.sort()

    return res_folders


def is_valid_folder(folder):
    """Check if all required files are present in `folder`.

    Parameters
    ----------
    folder : str
        path to the folder to validate

    Returns
    -------
    bool
        True if folder is valid
    """
    return yml_exists(folder) and imfs_exists(folder) and predictors_exists(folder)


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


def save_envelopes(envelopes, file_name, out_path):
    """Save imfs' instantaneous amplitudes to binary file.

    Parameters
    ----------
    envelopes : numpy ndarray
        imfs' instantaneous amplitudes
    file_name : str
        name of the file
    out_path : str
        where to save the file
    """
    f_imfs = open(os.path.join(out_path, "{}.imf".format(file_name)), "wb")
    pickle.dump(envelopes, f_imfs)
    f_imfs.close()


def summary_table(folders, comparison, table_name):
    """Comparison plots of results.

    Parameters
    ----------
    folders : list[str]
        paths to the files needed for the plots
    comparison : list[int]
        imfs for comparison
    table_name : str
        name of the output csv
    """
    if len(folders) == 0:
        return None

    folders_path = os.path.sep.join(folders[0].split(os.path.sep)[:-1])
    cpath = os.path.join(folders_path, defines.COMPARISON_FOLDER)
    if not os.path.isdir(cpath):
        os.makedirs(cpath, exist_ok=True)

    gps = []
    culprits = {}
    corrs = {}
    m_freqs = {}

    for i in comparison:
        culprits[i] = []
        corrs[i] = []
        m_freqs[i] = []

    for folder in folders:
        if is_valid_folder(folder):
            yf = YmlFile(folder)
            gps.append(yf.get_gps())

            for i in comparison:
                try:
                    culprits[i].append(yf.get_channel_of_imf(i))
                    corrs[i].append(yf.get_corr_of_imf(i))
                    m_freqs[i].append(yf.get_mean_freq_of_imf(i))
                except:
                    culprits[i].append(np.nan)
                    corrs[i].append(np.nan)
                    m_freqs[i].append(np.nan)

    df_dict = {"gps": gps}
    for i in comparison:
        df_dict[defines.summary_table_culprit_column_of_imf(i)] = culprits[i]
        df_dict[defines.summary_table_correlation_column_of_imf(i)] = corrs[i]
        df_dict[defines.summary_table_mean_frequency_column_of_imf(i)] = m_freqs[i]
    df = pd.DataFrame(df_dict)
    df.to_csv(os.path.join(cpath, table_name), index=False)


def omegagram_plot_name(imf, ext):
    """Name of the omegagram plots.

    Parameters
    ----------
    imf : int
        imf to which the omegagram corresponds
    ext : str
        plot extension

    Returns
    -------
    plot_name : str
        plot name
    """
    return "imf_{:d}_omegagram.{}".format(imf, ext)


def culprit_plot_name(imf, ext):
    """Name of the culprit plots.

    Parameters
    ----------
    imf : int
        imf to which the culprit corresponds
    ext : str
        plot extension

    Returns
    -------
    plot_name : str
        plot name
    """
    return "imf_{:d}_culprit.{}".format(imf, ext)


def combo_plot_name(imf_list, ext):
    """Name of the combo plots.

    Parameters
    ----------
    imf_list : list[int], list[str]
        list of imfs to which the culprit corresponds
    ext : str
        plot extension

    Returns
    -------
    plot_name : str
        plot name
    """
    return "combo_imf_{}_culprit.{}".format("+".join([str(i) for i in imf_list]), ext)


def seismic_plot_name(ext):
    """Name of the seismic plots.

    Parameters
    ----------
    ext : str
        plot extension

    Returns
    -------
    plot_name : str
        plot name
    """
    return "seismic.{}".format(ext)


def summary_freq_plot_name(ext):
    """Name of the summary plot by frequency.

        Parameters
        ----------
        ext : str
            plot extension

        Returns
        -------
        plot_name : str
            plot name
        """
    return "summary_freq.{}".format(ext)


def summary_chamber_plot_name(ext):
    """Name of the summary plot by chamber.

        Parameters
        ----------
        ext : str
            plot extension

        Returns
        -------
        plot_name : str
            plot name
        """
    return "summary_chamber.{}".format(ext)


def create_comparison_folder(where):
    """Create folder for comparison files.

        Parameters
        ----------
        where : str
            where to create the folder
    """
    cfold = os.path.join(where, defines.COMPARISON_FOLDER)
    if os.path.exists(cfold):
        os.system("rm -rf {}".format(cfold))
    os.system("mkdir {}".format(cfold))

