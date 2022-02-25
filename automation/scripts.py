def __imports(plots=False):
    imports = ["import os",
               "import argparse",
               "from gwadaptive_scattering.tools import scattered_light"]
    if plots:
        imports.append("from gwadaptive_scattering.utils import plot_utils")

    return imports


def __args():
    args = ["ap = argparse.ArgumentParser()",
            "ap.add_argument(\"--gps\", required=True, type=int)",
            "ap.add_argument(\"--seconds\", required=True, type=int)",
            "ap.add_argument(\"--target_channel\", required=True)",
            "ap.add_argument(\"--channels_list\", required=True)",
            "ap.add_argument(\"--opath\", required=True)",
            "ap.add_argument(\"--samp_freq\", required=True, type=float)",
            "ap.add_argument(\"--lowpass_freq\", required=True, type=float)",
            "args = vars(ap.parse_args())"]

    return args


def generate_algo_script(name, env, plots=False, imfs=None,
                         imf_thr=-1.0, omg_thr=-1.0, harmonics=None,
                         fig_ext="png", fig_size=None):
    """Generate the script for the main algorithm.

    Parameters
    ----------
    name : str
        name of the script
    env : str
        python3 environment shebang
    plots : bool, optional
        generate (or not) plots during the analysis (default : False)
    imfs : str, list[int], optional
        imfs to plot, can be "all", or list of integers (default : [1])
    imf_thr : float, optional
        correlation value above which to plot imfs (default : -1.0)
    omg_thr : float, optional
        correlation value above which to plot omegagrams (default : -1.0)
    harmonics : list[int]
        harmonics for the culprit (default : [1, 3, 5])
    fig_ext: str, optional
        plots extension (default : png)
    fig_size : tuple[int], optional
        figure size
    """
    f = open(name, "w")
    lines = "{}\n".format(env)
    lines += "\n".join(__imports(plots=plots))
    lines += "\n".join(__args())
    lines += "\nscattered_light.scattered_light(args[\"gps\"], args[\"seconds\"], args[\"target_channel\"], " \
             "args[\"channels_list\"], args[\"opath\"], args[\"lowpass_freq\"], fs=args[\"samp_freq\"])"

    if plots:
        if harmonics is None:
            harmonics = [1, 3, 5]
        if imfs is None:
            imfs = [1]
        lines += "\nfolder_name = [os.path.join(args[\"opath\"], \"{:d}\".format(args[\"gps\"]))]"
        lines += "\nplot_utils.plot_imfs(folder_name, {}, imf_thr={:f}, " \
                 "save_ext={}, figsize={})".format(imfs, imf_thr, fig_ext, fig_size)
        lines += "\nplot_utils.plot_omegagrams(folder_name, {}, omegagram_thr={:f}, harmonics={}, " \
                 "save_ext={}, figsize={})".format(imfs, omg_thr, harmonics, fig_ext, fig_size)

    f.write(lines)
    f.close()


def generate_comparison_script():
    pass


def generate_scripts_for_analysis():
    pass
