#  file_utils.py - this file is part of the gwasr package.
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
from ..common import defines


class YmlFile(object):
    """Class for the yml output file.

    Parameters
    ----------
    yml_file : str
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
        return self.content[defines.PARAMS_SECT_KEY][defines.TARGET_CH_KEY]

    def get_gps(self):
        """GPS times.

        Returns
        -------
        int
            starting GPS
        int
            ending GPS
        """
        gps = self.content[defines.PARAMS_SECT_KEY][defines.GPS_KEY]
        gps_start = int(gps.split(",")[0])
        gps_end = int(gps.split(",")[1])

        return gps_start, gps_end

    def get_sampling_frequency(self):
        """Sampling frequency.

        Returns
        -------
        float
            sampling frequency
        """
        return self.content[defines.PARAMS_SECT_KEY][defines.SAMP_FREQ_KEY]

    def get_channels_list(self):
        """Channels list.

        Returns
        -------
        str
            path to channels list file
        """
        return self.content[defines.PARAMS_SECT_KEY][defines.CH_LIST_KEY]

    def get_output_path(self):
        """Output path.

        Returns
        -------
        str
            path to where output analysis was saved
        """
        return self.content[defines.PARAMS_SECT_KEY][defines.OUT_PATH_KEY]

    def get_lowpass_frequency(self):
        """Lowpass frequency.

        Returns
        -------
        float
            lowpass frequency
        """
        return self.content[defines.PARAMS_SECT_KEY][defines.LOWPASS_FREQ_KEY]

    def get_scattering_factor(self):
        """Scattering factor.

        Returns
        -------
        int
            scattering factor
        """
        return self.content[defines.PARAMS_SECT_KEY][defines.SCATTERING_KEY]

    def get_smoothing_window(self):
        """Smoothing window.

        Returns
        -------
        int
            smoothing window
        """
        return self.content[defines.PARAMS_SECT_KEY][defines.SMOOTH_WIN_KEY]

    def get_max_corr_imf(self):
        """Number of max correlated imf.

        Returns
        -------
        int
            max correlated imf number
        """
        return self.content[defines.MAX_CORR_SECT_KEY][defines.IMF_KEY]

    def get_max_corr_channel(self):
        """Max correlated channel.

        Returns
        -------
        str
            max correlated channel
        """
        return self.content[defines.MAX_CORR_SECT_KEY][defines.CHANNEL_KEY]

    def get_max_corr(self):
        """Max correlation.

        Returns
        -------
        float
            max correlation
        """
        return self.content[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY]

    def get_max_corr_mean_freq(self):
        """Mean frequency of max correlated channel.

        Returns
        -------
        float
            mean frequency of max correlated channel
        """
        return self.content[defines.MAX_CORR_SECT_KEY][defines.MEAN_FREQ_KEY]

    def get_imfs_count(self):
        """Get number of found imfs.

        Returns
        -------
        int
            number of found imfs
        """
        return len(self.content[defines.CORR_SECT_KEY])

    def get_channel_of_imf(self, imf):
        """Most correlated channel with imf `imf`.

        Parameters
        ----------
        imf : int
            imf number

        Returns
        -------
        str
            most correlated channel with imf `imf`
        """
        return self.content[defines.CORR_SECT_KEY][imf - 1][defines.CHANNEL_KEY]

    def get_corr_of_imf(self, imf):
        """Correlation of imf `imf`.

        Parameters
        ----------
        imf : int
            imf number

        Returns
        -------
        float
            correlation of imf `imf`
        """
        return self.content[defines.CORR_SECT_KEY][imf - 1][defines.CORR_KEY]

    def get_mean_freq_of_imf(self, imf):
        """Mean frequency of imf `imf`.

        Parameters
        ----------
        imf : int
            imf number

        Returns
        -------
        float
            mean frequency of imf `imf`
        """
        return self.content[defines.CORR_SECT_KEY][imf - 1][defines.MEAN_FREQ_KEY]

    def write_parameters(self, gps, target_channel_name, channels_file, out_path,
                         fs, f_lowpass, n_scattering, smooth_win):
        """Write parameters section to file.

        Parameters
        ----------
        gps : str
            gps times (comma separated)
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
        self.content[defines.PARAMS_SECT_KEY][defines.GPS_KEY] = gps
        self.content[defines.PARAMS_SECT_KEY][defines.TARGET_CH_KEY] = target_channel_name
        self.content[defines.PARAMS_SECT_KEY][defines.CH_LIST_KEY] = channels_file
        self.content[defines.PARAMS_SECT_KEY][defines.OUT_PATH_KEY] = out_path
        self.content[defines.PARAMS_SECT_KEY][defines.SAMP_FREQ_KEY] = float(fs)
        self.content[defines.PARAMS_SECT_KEY][defines.LOWPASS_FREQ_KEY] = float(f_lowpass)
        self.content[defines.PARAMS_SECT_KEY][defines.SCATTERING_KEY] = int(n_scattering)
        self.content[defines.PARAMS_SECT_KEY][defines.SMOOTH_WIN_KEY] = int(smooth_win)

    def write_max_corr_section(self, n_imf, max_ch_str, max_corr, mean_freq):
        """Write max correlation section to file.

        Parameters
        ----------
        n_imf : int
            imf number
        max_ch_str : str
            channel with max correlation
        max_corr : float
            max correlation value
        mean_freq : float
            mean frequency of the channel with max correlation
        """
        self.content[defines.MAX_CORR_SECT_KEY] = {}
        self.content[defines.MAX_CORR_SECT_KEY][defines.IMF_KEY] = int(n_imf)
        self.content[defines.MAX_CORR_SECT_KEY][defines.CHANNEL_KEY] = max_ch_str
        self.content[defines.MAX_CORR_SECT_KEY][defines.CORR_KEY] = float(max_corr)
        self.content[defines.MAX_CORR_SECT_KEY][defines.MEAN_FREQ_KEY] = float(mean_freq)

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


# def load_yml(yml_path):
#     """Load yml file.
#
#     Parameters
#     ----------
#     yml_path : str
#         path to yml file
#
#     Returns
#     -------
#     dict
#         yml file
#     """
#     yml_file = os.path.join(yml_path, "output.yml")
#     with open(yml_file, "r") as f:
#         yml = yaml.safe_load(f)
#
#     return yml


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


def get_results_folders(results_path, sort=True, must_include=None):
    """Get list of folders with results from one or multiple analyses.

    Parameters
    ----------
    results_path : str
        path to results folders
    sort : bool
        sort folders (default : True)
    must_include : list[str]
        patterns of files that must be present in the
        folder, otherwise is discarded (default : None)

    Returns
    -------
    list[str]
        list of folders paths
    """
    res_folders = []
    for folder in os.listdir(results_path):
        curr_dir = os.path.join(results_path, folder)
        if os.path.isdir(curr_dir):
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
