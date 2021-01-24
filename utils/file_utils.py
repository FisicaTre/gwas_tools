import scipy.io
import yaml
import os
import glob
import pickle


def from_mat(mat_file, mat_col):
    mat = scipy.io.loadmat(mat_file)
    mat_data = mat[mat_col]
    
    return mat_data


def yml_exists(yml_path):
    yml_file = os.path.join(yml_path, "output.yml")
    if not os.path.exists(yml_file):
        print("File {} not found.".format(yml_file))
        return False

    return True


def imfs_exists(imfs_path):
    imfs_file = glob.glob(os.path.join(imfs_path, "*.imf"))
    if len(imfs_file) != 1:
        print("File {} not found or more than one `.imf` found.".format(imfs_file))
        return False

    return True


def predictors_exists(predictors_path):
    predictors_file = glob.glob(os.path.join(predictors_path, "*.predictors"))
    if len(predictors_file) != 1:
        print("File {} not found or more than one `.predictors` found.".format(predictors_file))
        return False

    return True


def load_yml(yml_path):
    yml_file = os.path.join(yml_path, "output.yml")
    with open(yml_file, "r") as f:
        yml = yaml.safe_load(f)

    return yml


def load_imfs(imfs_path):
    imfs_file = glob.glob(os.path.join(imfs_path, "*.imf"))
    f = open(imfs_file[0], "rb")
    imfs = pickle.load(f)
    f.close()

    return imfs


def load_predictors(predictors_path):
    predictors_file = glob.glob(os.path.join(predictors_path, "*.predictors"))
    f = open(predictors_file[0], "rb")
    predictors = pickle.load(f)
    f.close()

    return predictors
