#  file_utils.py - this file is part of the asr package,
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


import scipy.io
import yaml
import os
import glob
import pickle


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


def load_yml(yml_path):
    """Load yml file.
    
    Parameters
    ----------
    yml_path : str
        path to yml file
        
    Returns
    -------
    dict
        yml file
    """
    yml_file = os.path.join(yml_path, "output.yml")
    with open(yml_file, "r") as f:
        yml = yaml.safe_load(f)

    return yml


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
