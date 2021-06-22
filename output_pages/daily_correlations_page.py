#  daily_correlations_page.py - this file is part of the gwadaptive_scattering package.
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
import sys
from datetime import datetime
from astropy.time import Time
from ..html import html_builder as hb
from ..utils import file_utils


#SAVE_PATH = os.path.expandvars("/data/dev/web/detchar/daily_corr/")
SAVE_PLOTS_FOLDER = "plots"
SUMMARY_IMFS = 1
COLOR_THRESHOLD_MIN = 0.5
COLOR_THRESHOLD = 0.7
COPY_OR_MOVE = "cp"


def daily_correlations_page(res_path, date, tc_name, aux_ch, save_path):
    """Output page for daily correlations analysis.

    Parameters
    ----------
    res_path : str
        path containing analyses' results folders
    date : str
        date to be published on the page (year-month-day)
    tc_name : str
        target channel name
    aux_ch : str
        auxiliary channel name
    save_path : str
        path where to store the page files
    """
    res_folders = file_utils.get_results_folders(res_path, must_include=["*.png"])

    # title
    page_date = datetime.strptime(date, "%Y-%m-%d")
    curr_date = page_date.strftime('%Y%m%d')
    curr_folder = os.path.join(save_path, curr_date, "_".join(tc_name.split(":")))
    curr_plots_folder = os.path.join(curr_folder, SAVE_PLOTS_FOLDER)
    os.system("mkdir -p {}".format(curr_plots_folder))
    title = "Scattered light daily analysis ({})".format(page_date.strftime('%Y-%m-%d'))
    page = hb.HtmlBuilder(title=title, **{"style": "body { background-color: white; }"})

    # info
    page.addSection("Info")
    page.openDiv(**{"id_": "info-list"})

    info_dict = {
        "Environment": page.getFormattedCode(sys.prefix),
        "Target channel": tc_name,
        "Auxiliary channel": aux_ch
    }
    page.addBulletList(info_dict)
    page.closeDiv()

    # correlation plot
    if os.path.exists(os.path.join(res_path, "comparison")):
        page.addSection("Correlation between {} and {}".format(tc_name, aux_ch))

        for i in range(1, SUMMARY_IMFS + 1):
            ts_corr_name = os.path.join(res_path, "comparison", "imf_{}_ts_corr.png".format(i))
            if os.path.exists(ts_corr_name):
                page.openDiv(**{"id_": "imf-{}-summary".format(i)})
                page.addSubsection("Imf {}".format(i))

                ts_corr_to_save = os.path.join(curr_plots_folder, "imf-{}-ts-corr.png".format(i))

                os.system("{} {} {}".format(COPY_OR_MOVE, ts_corr_name, ts_corr_to_save))
                page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-ts-corr.png".format(i)),
                             "imf-{}-ts-corr".format(i))
                page.closeDiv()

    hour_dict = {}
    for ih in range(24):
        hour_dict[ih] = []

    for gf in res_folders:
        yf = file_utils.YmlFile(gf)
        gps_event = yf.get_gps()
        t1 = Time(gps_event, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        h = int("{}".format(t2).split(" ")[1].split(":")[0])
        hour_dict[h].append(gf)

    # results
    page.addSection("Single GPS results")

    page.openDiv(**{"id_": "info-color"})

    info_color_dict = {
        "Yellow": "Correlation between {} and {}".format(COLOR_THRESHOLD_MIN, COLOR_THRESHOLD),
        "Red": "Correlation greater than {}".format(COLOR_THRESHOLD)
    }
    page.addBulletList(info_color_dict)
    page.closeDiv()

    page.openDiv(**{"id_": "results"})

    for hh in hour_dict.keys():
        if len(hour_dict[hh]) > 0:
            date_h = "{} {:02d}:00:00".format(page_date.strftime('%Y-%m-%d'), hh)

            # main card
            page.openCard(date_h, "primary", "main-card-{:d}".format(hh))

            for gps_path in hour_dict[hh]:
                res_file = file_utils.YmlFile(gps_path)
                parameters = [
                    ("GPS", res_file.get_gps()),
                    ("Seconds", res_file.get_seconds()),
                    ("Event position", res_file.get_event_position()),
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
                for i in range(1, SUMMARY_IMFS + 1):
                    if res_file.get_imfs_count() >= i:
                        imfs_exists = os.path.exists(os.path.join(gps_path, "imf_{}_culprit.png".format(i)))
                        if imfs_exists:
                            imfs_data[i] = {}
                            imfs_data[i]["Culprit"] = res_file.get_channel_of_imf(i)
                            imfs_data[i]["Mean frequency"] = "{:.4f} Hz".format(res_file.get_mean_freq_of_imf(i))
                            imfs_data[i]["omegagram"] = os.path.exists(os.path.join(gps_path,
                                                                                    "imf_{:d}_omegagram.png".format(i)))

                            if res_file.get_corr_of_imf(i) >= COLOR_THRESHOLD:
                                above_thr = True

                gps_event = res_file.get_gps()
                res_id = "{:d}".format(gps_event)
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
                seconds = res_file.get_seconds()
                event_pos = res_file.get_event_position()
                gps_start = gps_event
                if event_pos == "center":
                    gps_start = gps_event - seconds // 2
                elif event_pos == "end":
                    gps_start = gps_event - seconds
                gps_end = gps_start + seconds
                page.parametersTable(parameters, int(gps_start), int(gps_end))

                # plots
                for i in range(1, SUMMARY_IMFS + 1):
                    if i in imfs_data.keys():
                        page.openDiv(**{"id_": "imf-{}-{}".format(i, res_id)})
                        page.addSubsection("Imf {}".format(i))

                        omegagram = imfs_data[i].pop("omegagram")

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

            page.closeCard()

    page.closeDiv()

    html_file = os.path.join(curr_folder, "index.html")
    page.savePage(html_file)
