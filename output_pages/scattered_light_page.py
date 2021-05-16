#  scattered_light_page.py - this file is part of the gwadaptive_scattering package.
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
import glob
import re
import sys
from datetime import datetime
from astropy.time import Time
from ..html import html_builder as hb
from ..utils import file_utils


#SAVE_PATH = os.path.expandvars("$HOME/public_html/daily/")
SAVE_PLOTS_FOLDER = "plots"
#SCRIPT_NAME = "./scattered_light.py"
#PLOTS_SCRIPT_NAME = "./scattered_light_plots.py"
#COMPARISON_SCRIPT_NAME = "./scattered_light_comparison.py"
#HTML_SCRIPT_NAME = "./scattered_light_page.py"
#PIPELINE_SCRIPT_NAME = "./gwasr_pipeline.py"
#SUMMARY_IMFS = 2
COLOR_THRESHOLD = 0.5
COPY_OR_MOVE = "cp"


def scattered_light_page(res_path, date, tc_name, ch_list_file, gps_file, save_path,
                         summary_imfs=2):
    """Output page for scattered light analysis.

    Parameters
    ----------
    res_path : str
        path containing analyses' results folders
    date : str
        date to be published on the page (year-month-day)
    tc_name : str
        target channel name
    ch_list_file : str
        channels list
    gps_file : str
        gps glitches list
    save_path : str
        path where to store the page files
    summary_imfs : int, optional
        imfs to show in the page
    """
    flags = []

    res_folders = file_utils.get_results_folders(res_path, must_include=["imf_*_culprit.png"])

    # title
    page_date = datetime.strptime(date, "%Y-%m-%d")
    curr_date = page_date.strftime('%Y%m%d')
    curr_folder = os.path.join(save_path, curr_date, "_".join(tc_name.split(":")))
    curr_plots_folder = os.path.join(curr_folder, SAVE_PLOTS_FOLDER)
    os.system("mkdir -p {}".format(curr_plots_folder))
    title = "Scattered light analysis ({})".format(page_date.strftime('%Y-%m-%d'))
    page = hb.HtmlBuilder(title=title, **{"style": "body { background-color: white; }"})

    # pipeline code
    #code = [
    #    PIPELINE_SCRIPT_NAME,
    #    "--target_channel {}".format(tc_name),
    #    "--channels_list {}".format(ch_list_file)
    #]
    #description = "The entire pipeline (analysis, plots, and summary page) " \
    #              "can be reproduced with the following command line:"
    #page.addCommandLineBlock(" ".join(code), description, "pipeline-command-line")

    # html generation code
    #code = [
    #    HTML_SCRIPT_NAME,
    #    "--ipath {}".format(res_path),
    #    "--date {}".format(date),
    #    "--target_channel {}".format(tc_name),
    #    "--channels_list {}".format(ch_list_file),
    #    "--gps_list {}".format(gps_file)
    #]
    #description = "This page can be reproduced with the following command line:"
    #page.addCommandLineBlock(" ".join(code), description, "page-command-line")

    # info
    page.addSection("Info")
    page.openDiv(**{"id_": "info-list"})

    if len(res_folders) > 0:
        os.system("{} {} {}".format(COPY_OR_MOVE, ch_list_file, curr_folder))
        os.system("{} {} {}".format(COPY_OR_MOVE, gps_file, curr_folder))

    info_dict = {
        "Environment": page.getFormattedCode(sys.prefix),
        "Target channel": tc_name,
        "Flags": " | ".join(flags),
        "GPS list": page.getFormattedLink(os.path.basename(gps_file),
                                          **{"href": os.path.basename(gps_file),
                                             "download": os.path.basename(gps_file)}),
        "Channels list": page.getFormattedLink(os.path.basename(ch_list_file),
                                               **{"href": os.path.basename(ch_list_file),
                                                  "download": os.path.basename(ch_list_file)})
    }
    page.addBulletList(info_dict)
    page.closeDiv()

    # results
    page.addSection("Results")
    page.openDiv(**{"id_": "results"})
    for gps_folder in res_folders:
        gps_path = os.path.join(res_path, gps_folder)
        res_file = file_utils.YmlFile(gps_path)
        parameters = [
            ("GPS", ",".join([str(g) for g in res_file.get_gps()])),
            ("Target channel", res_file.get_target_channel()),
            ("Channels list", res_file.get_channels_list()),
            ("Output path", res_file.get_output_path()),
            ("Sampling frequency", res_file.get_sampling_frequency()),
            ("Lowpass frequency", res_file.get_lowpass_frequency()),
            ("Scattering factor", res_file.get_scattering_factor()),
            ("Smoothing window", res_file.get_smoothing_window())
        ]
        imfs_data = {}
        above_thr = False
        for i in range(1, summary_imfs + 1):
            if res_file.get_imfs_count() >= i:
                imfs_data[i] = {}
                imfs_data[i]["Culprit"] = res_file.get_channel_of_imf(i)
                imfs_data[i]["Mean frequency"] = "{:.4f} Hz".format(res_file.get_mean_freq_of_imf(i))
                imfs_data[i]["omegagram"] = os.path.exists(os.path.join(gps_path, "imf_{:d}_omegagram.png".format(i)))
                imfs_data[i]["combo"] = ""

                regex = "[_+]{:d}[_+]".format(i)
                for cf in glob.glob(os.path.join(gps_path, "combo_imf_*_culprit.png")):
                    if re.search(regex, cf):
                        imfs_data[i]["combo"] = cf
                if res_file.get_corr_of_imf(i) >= COLOR_THRESHOLD:
                    above_thr = True

        gps_start, gps_end = res_file.get_gps()
        gps_event = (gps_start + gps_end) // 2
        res_id = "{:d}-{:d}".format(gps_start, gps_end)
        t1 = Time(gps_event, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        gps_date = "{} UTC (GPS: {:d})\n".format(t2, gps_event)

        if above_thr:
            color_key = "danger"
        else:
            color_key = "warning"

        # card
        page.openCard(gps_date, color_key, res_id)

        # parameters table
        page.parametersTable(parameters, int(gps_start), int(gps_end))

        # tool command line
        #code = [
        #    SCRIPT_NAME,
        #    "--gps {}".format(",".join([str(g) for g in res_file.get_gps()])),
        #    "--target_channel {}".format(res_file.get_target_channel()),
        #    "--channels_list {}".format(res_file.get_channels_list()),
        #    "--opath {}".format(os.path.sep.join(res_file.get_output_path().split(os.path.sep)[:-1])),
        #    "--samp_freq {:.2f}".format(res_file.get_sampling_frequency()),
        #    "--lowpass_freq {:.2f}".format(res_file.get_lowpass_frequency()),
        #    "--scattering {:d}".format(res_file.get_scattering_factor()),
        #    "--smooth_win {:d}".format(res_file.get_smoothing_window())
        #]
        #description = "This analysis can be reproduced with the following command line:"
        #page.addCommandLineBlock(" ".join(code), description, "command-line-{}".format(res_id))

        # plots command line
        #code = [
        #    PLOTS_SCRIPT_NAME,
        #    "--ipath {}".format(gps_path),
        #    "--single_folder True"
        #]
        #description = "Plots can be reproduced with the following command line:"
        #page.addCommandLineBlock(" ".join(code), description, "plots-command-line-{}".format(res_id))

        # plots
        for i in range(1, summary_imfs + 1):
            if i in imfs_data.keys():
                page.openDiv(**{"id_": "imf-{}-{}".format(i, res_id)})
                page.addSubsection("Imf {}".format(i))

                omegagram = imfs_data[i].pop("omegagram")
                combo_file = imfs_data[i].pop("combo")

                page.openDiv(**{"id_": "imf-{}-{}-info".format(i, res_id)})
                page.addBulletList(imfs_data[i])
                page.closeDiv()

                imf_plot_name = os.path.join(gps_path, "imf_{}_culprit.png".format(i))
                imf_to_save = os.path.join(curr_plots_folder, "imf-{}-{}.png".format(i, res_id))
                os.system("{} {} {}".format(COPY_OR_MOVE, imf_plot_name, imf_to_save))
                page.openDiv(**{"id_": "imf-{}-{}-plot".format(i, res_id)})
                page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-{}.png".format(i, res_id)),
                             "imf-{}-{}".format(i, res_id))
                page.closeDiv()

                if combo_file != "":
                    combo_to_save = os.path.join(curr_plots_folder, "combo-{}-{}.png".format(i, res_id))
                    os.system("{} {} {}".format(COPY_OR_MOVE, combo_file, combo_to_save))
                    page.openDiv(**{"id_": "combo-{}-{}-plot".format(i, res_id)})
                    page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "combo-{}-{}.png".format(i, res_id)),
                                 "combo-{}-{}".format(i, res_id))
                    page.closeDiv()

                if omegagram:
                    omegagram_plot_name = os.path.join(gps_path, "imf_{}_omegagram.png".format(i))
                    omegagram_to_save = os.path.join(curr_plots_folder, "omegagram-{}-{}.png".format(i, res_id))
                    os.system("{} {} {}".format(COPY_OR_MOVE, omegagram_plot_name, omegagram_to_save))
                    page.openDiv(**{"id_": "omegagram-{}-{}-plot".format(i, res_id)})
                    page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "omegagram-{}-{}.png".format(i, res_id)),
                                 "omegagram-{}-{}".format(i, res_id))
                    page.closeDiv()

                page.closeDiv()

        page.closeCard()

    page.closeDiv()

    # summary plots
    # if os.path.exists(os.path.join(res_path, "comparison")):
    #    page.addSection("Summary")

    #    # summary command line
    #    code = [
    #        COMPARISON_SCRIPT_NAME,
    #        "--ipath {}".format(res_path)
    #    ]
    #    description = "This summary can be reproduced with the following command line:"
    #    page.addCommandLineBlock(" ".join(code), description, "command-line-summary")

    #    # imfs
    #    for i in range(1, SUMMARY_IMFS + 1):
    #        page.openDiv(**{"id_": "imf-{}-summary".format(i)})
    #        page.addSubsection("Imf {}".format(i))

    #        summary_name = os.path.join(res_path, "comparison",
    #                                    "imf_{}_summary_{}.png".format(i, "_".join(tc_name.split(":"))))
    #        summary_name_first_batch = os.path.join(res_path, "comparison",
    #                                                "imf_{}_summary_{}_batch_1.png".format(i, "_".join(tc_name.split(":"))))
    #        mf_summary_name = os.path.join(res_path, "comparison",
    #                                       "imf_{}_summary_{}_mean_freq.png".format(i, "_".join(tc_name.split(":"))))
    #        corr_summary_name = os.path.join(res_path, "comparison",
    #                                         "imf_{}_corr_summary_{}.png".format(i, "_".join(tc_name.split(":"))))
    #        summary_to_save = os.path.join(curr_plots_folder, "imf-{}-summary.png".format(i))
    #        mf_summary_to_save = os.path.join(curr_plots_folder, "imf-{}-mf-summary.png".format(i))
    #        corr_summary_to_save = os.path.join(curr_plots_folder, "imf-{}-corr-summary.png".format(i))

    #        if os.path.exists(summary_name):
    #            os.system("cp {} {}".format(summary_name, summary_to_save))
    #            page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-summary.png".format(i)),
    #                         "imf-{}-summary".format(i))
    #        elif os.path.exists(summary_name_first_batch):
    #            os.system("cp {} {}".format(summary_name_first_batch, summary_to_save))
    #            page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-summary.png".format(i)),
    #                         "imf-{}-summary".format(i))
    #        if os.path.exists(mf_summary_name):
    #            os.system("cp {} {}".format(mf_summary_name, mf_summary_to_save))
    #            page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-mf-summary.png".format(i)),
    #                         "imf-{}-mf-summary".format(i))
    #        if os.path.exists(corr_summary_name):
    #            os.system("cp {} {}".format(corr_summary_name, corr_summary_to_save))
    #            page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-corr-summary.png".format(i)),
    #                         "imf-{}-corr-summary".format(i))
    #        page.closeDiv()

    html_file = os.path.join(curr_folder, "index.html")
    page.savePage(html_file)
