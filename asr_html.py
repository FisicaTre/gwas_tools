#  asr_html.py - this file is part of the asr package,
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

import os
import sys
import argparse
import glob
import re
from .utils import file_utils
from datetime import datetime
from datetime import timedelta
from astropy.time import Time
from .html import scattered_light_page as slp
from .common import defines


SAVE_PATH = os.path.expandvars("$HOME/public_html/daily/")
SAVE_PLOTS_FOLDER = "plots"
SCRIPT_NAME = "./scattered_light_tool.py"
PLOTS_SCRIPT_NAME = "./scattered_light_plots.py"
HTML_SCRIPT_NAME = "./scattered_light_html.py"
PIPELINE_SCRIPT_NAME = "./asr_pipeline.py"
SUMMARY_IMFS = 2
OMEGAGRAM_THR = 0.5

def generate_html(res_path, tc_name, ch_list_file, gspy_file, flags=[]):
    """Generate HTML page.
    
    Parameters
    ----------                                            
    res_path : str
        path to the analysis folders
    tc_name : str
        target channel name
    ch_list_file : str
        path to the channels list file
    gspy_file : str
        path to the gravity spy output file
    flags : list
        list of flags used for the analysis (default : [])
    """ 
    res_folders = []
    for folder in os.listdir(res_path):
        curr_dir = os.path.join(res_path, folder)
        if os.path.isdir(curr_dir):
            if file_utils.yml_exists(curr_dir):
                res_folders.append(folder)
    res_folders.sort()

    # title
    page_date = datetime.today() - timedelta(days=1)
    curr_date = page_date.strftime('%Y%m%d')
    curr_folder = os.path.join(SAVE_PATH, curr_date, "_".join(tc_name.split(":")))
    curr_plots_folder = os.path.join(curr_folder, SAVE_PLOTS_FOLDER)
    os.system("mkdir -p {}".format(curr_plots_folder))
    title = "Scattered light analysis ({})".format(page_date.strftime('%Y-%m-%d'))
    page = slp.ScatteredLightPage(title=title, **{"style": "body { background-color: white; }"})
    
    # pipeline code
    page.openDiv(id_="pipeline-command-line-section")
    page.addParagraph("The entire pipeline (analysis, plots, and summary page) can be reproduced with the following command line:", class_="mb-2")
    page.openDiv(id_="pipeline-command-line")
    code = [PIPELINE_SCRIPT_NAME]
    code.append("--target_channel {}".format(tc_name))
    code.append("--channels_list {}".format(ch_list_file))
    page.addCommandLine(" ".join(code))
    page.closeDiv()
    page.closeDiv()

    # html generation code
    page.openDiv(id_="page-command-line-section")
    page.addParagraph("This page can be reproduced with the following command line:", class_="mb-2")
    page.openDiv(id_="page-command-line")
    code = [HTML_SCRIPT_NAME]
    code.append("--ipath {}".format(res_path))
    page.addCommandLine(" ".join(code))
    page.closeDiv()
    page.closeDiv()
    
    # info
    page.addSection("Info", id_="info-section")
    page.openDiv(id_="info-list")
    
    if len(res_folders) > 0:
        os.system("cp {} {}".format(ch_list_file, curr_folder))
        os.system("cp {} {}".format(gspy_file, curr_folder))
    
    info_dict = {"Target channel": tc_name,
                 "Flags": " | ".join(flags),
                 "GPS list": page.getFormattedLink(os.path.basename(gspy_file),
                                                   **{"href": os.path.basename(gspy_file),
                                                      "download": os.path.basename(gspy_file)}),
                 "Channels list": page.getFormattedLink(os.path.basename(ch_list_file),
                                                        **{"href": os.path.basename(ch_list_file),
                                                           "download": os.path.basename(ch_list_file)})
                }
    page.addBulletList(info_dict)
    page.closeDiv()

    # results   
    page.addSection("Results", id_="results-section") 
    page.openDiv(id_="results")
    for gps_folder in res_folders:
        gps_path = os.path.join(res_path, gps_folder)
        res_file = file_utils.load_yml(gps_path)
        parameters = res_file[defines.PARAMS_SECT_KEY]
        imfs_data = {}
        above_thr = False
        for i in range(1, SUMMARY_IMFS + 1):
            try:
                imfs_data[i] = {}
                imfs_data[i]["Culprit"] = res_file[defines.CORR_SECT_KEY][i - 1][defines.CHANNEL_KEY]
                #imfs_data[i]["Correlation"] = "{:.4f}".format(res_file[defines.CORR_SECT_KEY][i - 1][defines.CORR_KEY])
                imfs_data[i]["Mean frequency"] = "{:.4f} Hz".format(res_file[defines.CORR_SECT_KEY][i - 1][defines.MEAN_FREQ_KEY])
                imfs_data[i]["omegagram"] = os.path.exists(os.path.join(gps_path, "imf_{:d}_omegagram.png".format(i)))
                imfs_data[i]["combo"] = ""
                
                regex = "[_+]{:d}[_+]".format(i)
                for cf in glob.glob(os.path.join(gps_path, "combo_imf_*_culprit.png")):
                    if re.search(regex, cf):
                        imfs_data[i]["combo"] = cf
                if res_file[defines.CORR_SECT_KEY][i - 1][defines.CORR_KEY] >= 0.5:
                    above_thr = True
            except:
                continue
        
        gps_start = gps_folder.split("_")[0]
        gps_end = gps_folder.split("_")[1]
        gps_event = str((int(gps_start) + int(gps_end)) // 2)
        res_id = "{}-{}".format(gps_start, gps_end)
        t1 = Time(gps_event, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        gps_date = "{} UTC (GPS: {})\n".format(t2, gps_event)
        
        if above_thr:
            color_key = "danger"
        else:
            color_key = "warning"
        
        page.openDiv(class_="card border-{} mb-1 shadow-sm".format(color_key))
        
        # header
        page.openDiv(class_="card-header text-white bg-{}".format(color_key))
        page.addLink("{}".format(gps_date),
                     **{"class": "btn card-link cis-link",
                        "data-toggle": "collapse",
                        "data-target": "#{}".format(res_id)})
        page.closeDiv()
        
        # body
        page.openDiv(id_=res_id, class_="collapse")
        page.openDiv(class_="card-body")
        page.parametersTable(parameters, int(gps_start), int(gps_end))
        
        page.openDiv(id_="command-line-{}".format(res_id))
        page.addParagraph("This analysis can be reproduced with the following command line:", class_="mb-2")
        code = [SCRIPT_NAME]
        for key in parameters.keys():
            code.append("--" + key)
            code.append(str(parameters[key]))
        page.addCommandLine(" ".join(code))
        page.closeDiv()
        
        page.openDiv(id_="plots-command-line-{}".format(res_id))
        page.addParagraph("Plots can be reproduced with the following command line:", class_="mb-2")
        code = [PLOTS_SCRIPT_NAME]
        code.append("--ipath {}".format(gps_path))
        code.append("--imfs_to_plot {}".format(",".join([str(i) for i in range(1, SUMMARY_IMFS + 1)])))
        code.append("--omegagram_thr {:.2f}".format(OMEGAGRAM_THR))
        page.addCommandLine(" ".join(code))
        page.closeDiv()
        
        for i in range(1, SUMMARY_IMFS + 1):
            if i in imfs_data.keys():
                page.openDiv(id_="imf-{}-{}".format(i, res_id))
                page.addSubsection("Imf {}".format(i), id_="imf-{}-{}-section".format(i, res_id))
                
                omegagram = imfs_data[i].pop("omegagram")
                combo_file = imfs_data[i].pop("combo")
                
                page.openDiv(id_="imf-{}-{}-info".format(i, res_id))
                page.addBulletList(imfs_data[i])
                page.closeDiv()
            
                imf_plot_name = os.path.join(gps_path, "imf_{}_culprit.png".format(i))
                imf_to_save = os.path.join(curr_plots_folder, "imf-{}-{}.png".format(i, res_id))
                os.system("cp {} {}".format(imf_plot_name, imf_to_save))
                page.openDiv(id_="imf-{}-{}-plot".format(i, res_id))
                page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-{}.png".format(i, res_id)), "imf-{}-{}".format(i, res_id))
                page.closeDiv()
                
                if combo_file != "":
                    combo_to_save = os.path.join(curr_plots_folder, "combo-{}-{}.png".format(i, res_id))
                    os.system("cp {} {}".format(combo_file, combo_to_save))
                    page.openDiv(id_="combo-{}-{}-plot".format(i, res_id))
                    page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "combo-{}-{}.png".format(i, res_id)), "combo-{}-{}".format(i, res_id))
                    page.closeDiv()
                
                if omegagram:
                    omegagram_plot_name = os.path.join(gps_path, "imf_{}_omegagram.png".format(i))
                    omegagram_to_save = os.path.join(curr_plots_folder, "omegagram-{}-{}.png".format(i, res_id))
                    os.system("cp {} {}".format(omegagram_plot_name, omegagram_to_save))
                    page.openDiv(id_="omegagram-{}-{}-plot".format(i, res_id))
                    page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "omegagram-{}-{}.png".format(i, res_id)), "omegagram-{}-{}".format(i, res_id))
                    page.closeDiv()
            
                page.closeDiv()
            
        page.closeDiv()
        page.closeDiv()
        page.closeDiv()
        
    page.closeDiv()
    
    # summary plots
    if os.path.exists(os.path.join(res_path, "comparison")):
        page.addSection("Summary", id_="summary-section")
        
        # command line
        page.openDiv(id_="command-line-summary")
        page.addParagraph("This summary can be reproduced with the following command line:", class_="mb-2")
        
        code = [PLOTS_SCRIPT_NAME]
        code.append("--ipath {}".format(res_path))
        code.append("--comparison {}".format(",".join([str(i) for i in range(1, SUMMARY_IMFS + 1)])))
        page.addCommandLine(" ".join(code))
        
        page.closeDiv()
        
        # imfs
        for i in range(1, SUMMARY_IMFS + 1):
            page.openDiv(id_="imf-{}-summary".format(i))
            page.addSubsection("Imf {}".format(i), id_="imf-{}-section-summary".format(i))
            summary_name = os.path.join(res_path, "comparison", "imf_{}_summary_{}.png".format(i, tc_name))
            corr_summary_name = os.path.join(res_path, "comparison", "imf_{}_corr_summary_{}.png".format(i, tc_name))
            summary_to_save = os.path.join(curr_plots_folder, "imf-{}-summary-{}.png".format(i, res_id))
            corr_summary_to_save = os.path.join(curr_plots_folder, "imf-{}-corr-summary-{}.png".format(i, res_id))
            os.system("cp {} {}".format(summary_name, summary_to_save))
            os.system("cp {} {}".format(corr_summary_name, corr_summary_to_save))
            page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-summary-{}.png".format(i, res_id)), "imf-{}-summary".format(i))
            page.addPlot(os.path.join(SAVE_PLOTS_FOLDER, "imf-{}-corr-summary-{}.png".format(i, res_id)), "imf-{}-corr-summary".format(i))
            page.closeDiv()
    
    html_file = os.path.join(curr_folder, "index.html")
    page.savePage(html_file)
    
