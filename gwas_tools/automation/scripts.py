from ..common import defines


def __imports(mode, plots=False, html=False, corr=False):
    imports = ["\nimport os",
               "import argparse"]
    if mode == "tool":
        imports.append("from gwas_tools.algorithms import scattered_light")
        if plots:
            imports.append("from gwas_tools.utils import plot_utils")
    elif mode == "comparison":
        imports.append("from gwas_tools.utils import file_utils, plot_utils, signal_utils")
        if html:
            if corr:
                imports.append("from gwas_tools.summary_pages import daily_correlations_page")
            else:
                imports.append("from gwas_tools.summary_pages import scattered_light_page")
    imports.append("\n")

    return imports


def __args(mode, corr=False):
    args = None
    if mode == "tool":
        args = ["\nap = argparse.ArgumentParser()",
                "ap.add_argument(\"--gps\", required=True, type=int)",
                "ap.add_argument(\"--seconds\", required=True, type=int)",
                "ap.add_argument(\"--target_channel\", required=True)",
                "ap.add_argument(\"--channels_list\", required=True)",
                "ap.add_argument(\"--opath\", required=True)",
                "ap.add_argument(\"--samp_freq\", required=True, type=float)",
                "ap.add_argument(\"--lowpass_freq\", required=True)",
                "args = vars(ap.parse_args())"]
    elif mode == "comparison":
        if corr:
            args = ["\nap = argparse.ArgumentParser()",
                    "ap.add_argument(\"--ipath\", required=True)",
                    "ap.add_argument(\"--date\", required=True)",
                    "ap.add_argument(\"--target_channel\", required=True)",
                    "ap.add_argument(\"--aux_channel\", required=True)",
                    "args = vars(ap.parse_args())"]
        else:
            args = ["\nap = argparse.ArgumentParser()",
                    "ap.add_argument(\"--ipath\", required=True)",
                    "ap.add_argument(\"--date\", required=True)",
                    "ap.add_argument(\"--target_channel\", required=True)",
                    "ap.add_argument(\"--channels_list\", required=True)",
                    "ap.add_argument(\"--gps_list\", required=True)",
                    "args = vars(ap.parse_args())"]
    args.append("\n")

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
    imfs : list[int], optional
        imfs to plot (default : None, i.e. all the imfs)
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
    mode = "tool"
    f = open(name, "w")
    lines = "{}\n".format(env)
    lines += "\n".join(__imports(mode, plots=plots))
    lines += "\n".join(__args(mode))
    lines += "\nscattered_light.scattered_light(args[\"gps\"], args[\"seconds\"], args[\"target_channel\"], " \
             "args[\"channels_list\"], args[\"opath\"], args[\"lowpass_freq\"], fs=args[\"samp_freq\"], " \
             "check_lock=True)"

    if plots:
        if harmonics is None:
            harmonics = [1, 3, 5]
        if imfs is None:
            imfs = [1]
        lines += "\nfolder_name = [args[\"opath\"]]"
        lines += "\nplot_utils.plot_imfs(folder_name, imfs_to_plot={}, imf_thr={:f}, " \
                 "save_ext=\"{}\", figsize={})".format(imfs, imf_thr, fig_ext, fig_size)
        lines += "\nplot_utils.plot_omegagrams(folder_name, imfs_to_plot={}, omegagram_thr={:f}, harmonics={}, " \
                 "save_ext=\"{}\", figsize={})".format(imfs, omg_thr, harmonics, fig_ext, fig_size)

    f.write(lines)
    f.close()


def generate_algo_corr_script(name, env, plots=False, imfs=None,
                              imf_thr=-1.0, fig_ext="png", fig_size=None):
    """Generate the script for the main correlation algorithm.

    Parameters
    ----------
    name : str
        name of the script
    env : str
        python3 environment shebang
    plots : bool, optional
        generate (or not) plots during the analysis (default : False)
    imfs : list[int], optional
        imfs to plot (default : None, i.e. all the imfs)
    imf_thr : float, optional
        correlation value above which to plot imfs (default : -1.0)
    fig_ext: str, optional
        plots extension (default : png)
    fig_size : tuple[int], optional
        figure size
    """
    mode = "tool"
    f = open(name, "w")
    lines = "{}\n".format(env)
    lines += "\n".join(__imports(mode, plots=plots))
    lines += "\n".join(__args(mode))
    lines += "\nscattered_light.scattered_light(args[\"gps\"], args[\"seconds\"], args[\"target_channel\"], " \
             "args[\"channels_list\"], args[\"opath\"], args[\"lowpass_freq\"], fs=args[\"samp_freq\"], " \
             "max_imf=3, check_lock=True, combos=True, seismic=True)"

    if plots:
        if imfs is None:
            imfs = [1]
        lines += "\nfolder_name = [args[\"opath\"]]"
        lines += "\nplot_utils.plot_combinations(folder_name, imfs_to_plot={}, imf_thr={:f}, " \
                 "save_ext=\"{}\", figsize={})".format(imfs, imf_thr, fig_ext, fig_size)

    f.write(lines)
    f.close()


def generate_comparison_script(name, env, html=False, imfs=None, summary_imfs=1,
                               imf_thr=-1.0, fig_ext="png", fig_size=None, html_abs_path=""):
    """Generate the script for the analysis' results comparison.

    Parameters
    ----------
    name : str
        name of the script
    env : str
        python3 environment shebang
    html : bool, optional
        if True, the summary web page is generated (default : False)
    imfs : list[int], optional
        imfs to print in the summary table (default : None)
    summary_imfs : int, optional
        imfs to show in the summary web page for each day (default : 1)
    imf_thr : float, optional
        correlation value above which to consider a line in the summary table (default : -1.0)
    fig_ext: str, optional
        plots extension (default : png)
    fig_size : tuple[int], optional
        figure size
    html_abs_path : str, optional
        absolute path for files and images in the html page (default : "")
    """
    mode = "comparison"
    f = open(name, "w")
    lines = "{}\n".format(env)
    lines += "\n".join(__imports(mode, html=html))
    lines += "\n".join(__args(mode))
    lines += "\nres_folders = file_utils.get_results_folders(args[\"ipath\"])"
    lines += "\nfile_utils.create_comparison_folder(args[\"ipath\"])"

    if imfs is None:
        imfs = [1, 2]
    lines += "\nfile_utils.summary_table(res_folders, {}, \"{}\")".format(imfs, defines.SUMMARY_NAME)

    if html:
        lines += "\nplot_utils.plot_summaries(os.path.join(args[\"ipath\"], \"{}\", \"{}\"), " \
                 "imf_thr={:f}, save_ext=\"{}\", figsize={})".format(defines.COMPARISON_FOLDER,
                                                                     defines.SUMMARY_NAME,
                                                                     imf_thr, fig_ext, fig_size)
        lines += "\nscattered_light_page.generate_web_page(args[\"ipath\"], args[\"date\"], " \
                 "args[\"target_channel\"], args[\"channels_list\"], args[\"gps_list\"], " \
                 "summary_imfs={:d}, prepend_path=\"{}\")".format(summary_imfs, html_abs_path)

    f.write(lines)
    f.close()


def generate_comparison_corr_script(name, env, html=False, imfs=None, html_abs_path=""):
    """Generate the script for the correlation analysis' results comparison.

    Parameters
    ----------
    name : str
        name of the script
    env : str
        python3 environment shebang
    html : bool, optional
        if True, the summary web page is generated (default : False)
    imfs : list[int], optional
        imfs to print in the summary table (default : None)
    html_abs_path : str, optional
        absolute path for files and images in the html page (default : "")
    """
    mode = "comparison"
    f = open(name, "w")
    lines = "{}\n".format(env)
    lines += "\n".join(__imports(mode, html=html, corr=True))
    lines += "\n".join(__args(mode, corr=True))
    lines += "\nres_folders = file_utils.get_results_folders(args[\"ipath\"])"
    lines += "\nfile_utils.create_comparison_folder(args[\"ipath\"])"

    if imfs is None:
        imfs = [1, 2]
    lines += "\nfile_utils.summary_table(res_folders, {}, \"{}\")".format(imfs, defines.SUMMARY_NAME)

    if html:
        lines += "\nplot_utils.plot_seismic_data(args[\"ipath\"], " \
                 "signal_utils.get_ifo_of_channel(args[\"target_channel\"]))"
        lines += "\ndaily_correlations_page.generate_web_page(args[\"ipath\"], args[\"date\"], " \
                 "args[\"target_channel\"], args[\"aux_channel\"], prepend_path=\"{}\")".format(html_abs_path)

    f.write(lines)
    f.close()


def generate_condor_post_script(name):
    """Generate the condor post script to allow main jobs failure.

    Parameters
    ----------
    name : str
        name of the condor post script
    """
    f = open(name, "w")
    f.write("#!/bin/bash\n\necho \"$1 failed!\"\n\nexit 0")
    f.close()
