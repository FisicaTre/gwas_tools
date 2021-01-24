import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from astropy.time import Time
from scipy.stats import pearsonr
from gwpy.timeseries import TimeSeries
from gwpy.segments import Segment
from gwdetchar.scattering import get_fringe_frequency


def normalize(vec):
    norm_vec = (vec - np.nanmean(vec)) / np.nanstd(vec)
    
    return norm_vec


def plot_imf(pred, pred_name, imf_ia, imf_ia_name, gps1, samp_freq, title, plot_name, save_path):
    """
    Parameters:
    -----------
    pred : numpy array
        predictor
    pred_name : str
        predictor name
    imf_ia : numpy array
        imf instantaneous amplitude
    imf_ia_name : str
        imf instantaneous amplitude name
    gps1 : int
        gps start
    samp_freq : int
        sampling frequency
    title : str
        plot title
    plot_name : str
        plot name
    save_path : str
        save path
    """
    if gps1 is not None:
        t1 = Time(gps1, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        gps_date = "GPS: {:d} | t$_0$: {} UTC\n".format(gps1, t2)
    else:
        gps_date = ""
    x = np.arange(0, len(pred) / samp_freq, 1 / samp_freq, dtype=float)
    font_size = 16

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    l1 = ax1.plot(x[:len(pred)], pred, "r-", label=pred_name)
    l2 = ax2.plot(x[:len(imf_ia)], imf_ia, "b-", label=imf_ia_name)
    ls = l1 + l2
    leg = plt.legend(ls, [l.get_label() for l in ls], ncol=1, loc=2)# bbox_to_anchor=[0, 1], loc="lower left")
    #plt.legend(ls, [l.get_label() for l in ls], ncol=1, bbox_to_anchor=[0, 1], loc="lower left")
    plt.draw()
    leg_height_ratio = leg.get_window_extent().height / ax1.get_window_extent().height
    
    ax1.set_ylim(0, np.max(pred) + (np.max(pred) - np.min(pred)) * (leg_height_ratio + 0.1))
    ax2.set_ylim(0, np.max(imf_ia) + (np.max(imf_ia) - np.min(imf_ia)) * (leg_height_ratio + 0.1))
    #ax1.set_ylim(0, np.max(pred) + (np.max(pred) - np.min(pred)) * 0.15)
    #ax2.set_ylim(0, np.max(imf_ia) + (np.max(imf_ia) - np.min(imf_ia)) * 0.15)
    ax1.set_ylabel("Predictor [$Hz$]", color="r", fontsize=font_size)
    ax1.set_xlabel("t - t$_0$ [s]\n" + gps_date, fontsize=font_size)
    ax2.set_ylabel("IA", color="b", fontsize=font_size)
    ax1.grid(True)
    ax2.grid(False)
    #txt_css = {"backgroundcolor": "black", "color":  "white"}
    #at = AnchoredText(s=title, loc=2, prop=txt_css)
    #ax1.add_artist(at)
    plt.title(title)
    plt.savefig(os.path.join(save_path, plot_name + ".png"), bbox_inches="tight", dpi=300)
    plt.close("all")
    
    
def plot_combinations(plot_channels, ias, predictors, target_channel_name, gps1, samp_freq, out_path, thr=-1.0):
    """
    Parameters:
    -----------
    plot_channels : list of str
        channels names
    ias : numpy array
        imf instantaneous amplitudes matrix
    predictors : numpy array
        predictors matrix
    target_channel_name : str
        name of the channel from which `ias` are computed
    gps1 : int
        gps start
    samp_freq : int
        sampling frequency
    out_path : str
        save path
    thr : float
        correlation threshold for plots
    """
    seen = set()
    uniq = [x for x in plot_channels if x not in seen and not seen.add(x)]
    for u in uniq:
        plot_idxs = np.where(plot_channels == u)[0]
        if len(plot_idxs) != 1:
            sum_envelope = np.sum(ias[:, plot_idxs], axis=1)
            c = pearsonr(sum_envelope, predictors[:, plot_idxs[0]])[0]
            if c > thr:
                plot_imf(predictors[:, plot_idxs[0]], u, sum_envelope, target_channel_name + " (combo)",
                         gps1, samp_freq, "$\\rho$ = {:.4f}".format(c),
                         "combo_imf_{}_culprit".format("+".join([str(i + 1) for i in plot_idxs])), out_path)
                         #gps1, samp_freq, "$\\rho$ = {:.4f}".format(c), "combo_{}".format(u), out_path)
            
            
def plot_omegagram_download(pred, target_name, gps1, gps2, plot_name, save_path, norm=False, harmonics=[1, 2, 3, 4, 5]):
    """
    Parameters:
    -----------
    pred : numpy array
        predictor
    target_name : str
        name of the channel from which compute the omegagram
    gps1 : int
        gps start
    gps2 : int
        gps end
    plot_name : str
        plot name
    save_path : str
        save path
    """
    if gps1 is not None and gps2 is not None:
        epoch = (gps1 + gps2) // 2
        if norm:
            correction_factor = 10.0 / np.max(pred)
        else:
            correction_factor = 1.0
            
        line_width = 1
        colormap = "viridis"
        plot_t_min = (gps1 - gps2) // 2
        plot_t_max = (gps2 - gps1) // 2
        #plot_t_step = 5
        #plot_f_min = 6
        #plot_f_max = 50
        #plot_f_step = 10
        
        t1 = Time(epoch, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        gps_date = "GPS: {:d} | t$_0$: {} UTC\n".format(epoch, t2)

        ts_l2 = pred * correction_factor
        ts = TimeSeries.get(target_name, gps1, gps2).astype("float64")
        ts.times = ts.times.value - epoch
        tsq = ts.q_transform(outseg=Segment(plot_t_min, plot_t_max), tres=0.2, fres=0.2, whiten=True)
    
        plot_f_min = tsq.yindex[0].value
        plot_f_max = np.max([tsq.yindex[-1].value, np.max(ts_l2) * np.max(harmonics)])
        if plot_f_max > 50:
            plot_f_max = 50
    
        plot = plt.figure()
        ax = plot.gca()
        pcm = ax.imshow(tsq, vmin=0, vmax=15, cmap=colormap)
        ax.set_ylabel("Frequency [Hz]", fontsize=16)
        ax.set_ylim(plot_f_min, plot_f_max)
        ax.set_xlim(plot_t_min, plot_t_max)
        ax.grid(False)
        for i in harmonics:
            ax.plot(np.linspace(plot_t_min, plot_t_max, len(ts_l2)), ts_l2 * i, color="r", lw=line_width)
        ax.set_xlim(plot_t_min, plot_t_max)
        ax.set_ylim(plot_f_min, plot_f_max)
        #ax.set_xticks(np.arange(plot_t_min, plot_t_max + plot_t_step, plot_t_step))
        #ax.set_xticklabels(["{:d}".format(val) for val in np.arange(plot_t_min, plot_t_max + plot_t_step, plot_t_step)])
        #ax.set_yticks(list(np.arange(plot_f_min, plot_f_max + plot_f_step, plot_f_step)))
        #ax.set_yticklabels(["{:d}".format(val) for val in np.arange(plot_f_min, plot_f_max + plot_f_step, plot_f_step)])
        ax.set_xlabel("t - t$_0$ [s]", fontsize=16)
        ax.set_title(gps_date, fontsize=16)
        cbar = ax.colorbar(clim=(0, 15), location="right")
        cbar.set_label("Normalized energy", fontsize=16)
        plt.savefig(os.path.join(save_path, plot_name + ".png"), bbox_inches="tight", dpi=300)
        plt.close("all")
    else:
        print("GPS coordinates are None, cannot plot omegagram.")
        
        
def plot_omegagram(pred, target, gps1, gps2, fs, plot_name, save_path, norm=False, harmonics=[1, 2, 3, 4, 5]):
    """
    Parameters:
    -----------
    pred : numpy array
        predictor
    target : numpy array
        channel from which compute the omegagram
    gps1 : int
        gps start
    gps2 : int
        gps end
    plot_name : str
        plot name
    save_path : str
        save path
    """
    if gps1 is not None and gps2 is not None:
        epoch = (gps1 + gps2) // 2
        if norm:
            correction_factor = 10.0 / np.max(pred)
        else:
            correction_factor = 1.0
            
        line_width = 1
        colormap = "viridis"
        plot_t_min = (gps1 - gps2) // 2
        plot_t_max = (gps2 - gps1) // 2
        #plot_t_step = 5
        #plot_f_min = 0
        #plot_f_max = 50
        #plot_f_step = 10
        
        t1 = Time(epoch, format="gps")
        t2 = Time(t1, format="iso", scale="utc")
        gps_date = "GPS: {:d} | t$_0$: {} UTC\n".format(epoch, t2)

        ts_l2 = pred * correction_factor
        #ts = TimeSeries(target, t0=gps1, sample_rate=fs, dtype=float)#.astype("float64")
        ts = TimeSeries(target, times=np.arange(gps1, gps1 + len(target) / fs, 1 / fs, dtype=float), dtype=float)
        ts.times = ts.times.value - epoch
        tsq = ts.q_transform(outseg=Segment(plot_t_min, plot_t_max), tres=0.2, fres=0.2, whiten=True)
        
        plot_f_min = tsq.yindex[0].value
        plot_f_max = np.max([tsq.yindex[-1].value, np.max(ts_l2) * np.max(harmonics)])
        if plot_f_max > 50:
            plot_f_max = 50
    
        plot = plt.figure()
        ax = plot.gca()
        pcm = ax.imshow(tsq, vmin=0, vmax=15, cmap=colormap)
        ax.set_ylabel("Frequency [Hz]", fontsize=16)
        ax.set_ylim(plot_f_min, plot_f_max)
        ax.set_xlim(plot_t_min, plot_t_max)
        ax.grid(False)
        for i in harmonics:
            ax.plot(np.linspace(plot_t_min, plot_t_max, len(ts_l2)), ts_l2 * i, color="r", lw=line_width)
        ax.set_xlim(plot_t_min, plot_t_max)
        ax.set_ylim(plot_f_min, plot_f_max)
        #ax.set_xticks(np.arange(plot_t_min, plot_t_max + plot_t_step, plot_t_step))
        #ax.set_xticklabels(["{:d}".format(val) for val in np.arange(plot_t_min, plot_t_max + plot_t_step, plot_t_step)])
        #ax.set_yticks(list(np.arange(plot_f_min, plot_f_max + plot_f_step, plot_f_step)))
        #ax.set_yticklabels(["{:d}".format(val) for val in np.arange(plot_f_min, plot_f_max + plot_f_step, plot_f_step)])
        ax.set_xlabel("t - t$_0$ [s]", fontsize=16)
        ax.set_title(gps_date, fontsize=16)
        cbar = ax.colorbar(clim=(0, 15), location="right")
        cbar.set_label("Normalized energy", fontsize=16)
        plt.savefig(os.path.join(save_path, plot_name + ".png"), bbox_inches="tight", dpi=300)
        plt.close("all")
    else:
        print("GPS coordinates are None, cannot plot omegagram.")
        

def plot_imfs_summary(culprits, title, plot_name, save_path, dsort=True, batch=10):
    """
    Parameters:
    -----------
    culprits : list of str
        channels names
    title : str
        plot title
    plot_name : str
        plot name
    save_path : str
        save path
    """
    culprits = np.array(culprits)
    seen = set()
    uniq = [x for x in culprits if x not in seen and not seen.add(x)]
    counts = []
    for u in uniq:
        counts.append(len(np.where(culprits == u)[0]))
    
    if dsort:
        uniq = [x for _,x in sorted(zip(counts, uniq), reverse=True)]
        counts = sorted(counts, reverse=True)
    
    if len(counts) <= batch:
        plt.figure()
        plt.bar(np.arange(len(counts)), counts, width=0.8)
        plt.xticks(np.arange(len(counts)), uniq, rotation=45, horizontalalignment="right")
        plt.title(title)
        plt.savefig(os.path.join(save_path, plot_name + ".png"), bbox_inches="tight", dpi=300)
        plt.close("all")
    else:
        n = len(counts) // batch
        for i in range(n):
            plt.figure()
            plt.bar(np.arange(len(counts[i*batch:(i+1)*batch])), counts[i*batch:(i+1)*batch], width=0.8)
            plt.xticks(np.arange(len(counts[i*batch:(i+1)*batch])), uniq[i*batch:(i+1)*batch],
                       rotation=45, horizontalalignment="right")
            plt.title(title)
            plt.savefig(os.path.join(save_path, plot_name + "_batch_" + str(i + 1) + ".png"), bbox_inches="tight", dpi=300)
            plt.close("all")
        if len(counts) % batch != 0:
            plt.figure()
            plt.bar(np.arange(len(counts[n*batch:])), counts[n*batch:], width=0.8)
            plt.xticks(np.arange(len(counts[n*batch:])), uniq[n*batch:],
                       rotation=45, horizontalalignment="right")
            plt.title(title)
            plt.savefig(os.path.join(save_path, plot_name + "_batch_" + str(n + 1) + ".png"), bbox_inches="tight", dpi=300)
            plt.close("all")
    
    
def plot_corr_summary(gps_list, corr_list, title, plot_name, save_path, batch=10):
    """
    Parameters:
    -----------
    gps_list : list of str
        gps
    corr_list : list of float
        correlations
    title : str
        plot title
    plot_name : str
        plot name
    save_path : str
        save path
    """
    if len(gps_list) <= batch:
        plt.figure()
        plt.bar(np.arange(len(gps_list)), corr_list, width=0.8)
        plt.ylim(-1, 1)
        plt.xticks(np.arange(len(gps_list)), gps_list, rotation=45, horizontalalignment="right")
        plt.title(title)
        plt.savefig(os.path.join(save_path, plot_name + ".png"), bbox_inches="tight", dpi=300)
        plt.close("all")
    else:
        n = len(gps_list) // batch
        for i in range(n):
            plt.figure()
            plt.bar(np.arange(len(gps_list[i*batch:(i+1)*batch])), corr_list[i*batch:(i+1)*batch], width=0.8)
            plt.ylim(-1, 1)
            plt.xticks(np.arange(len(gps_list[i*batch:(i+1)*batch])), gps_list[i*batch:(i+1)*batch],
                       rotation=45, horizontalalignment="right")
            plt.title(title)
            plt.savefig(os.path.join(save_path, plot_name + "_batch_" + str(i + 1) + ".png"), bbox_inches="tight", dpi=300)
            plt.close("all")
        if len(gps_list) % batch != 0:
            plt.figure()
            plt.bar(np.arange(len(gps_list[n*batch:])), corr_list[n*batch:], width=0.8)
            plt.ylim(-1, 1)
            plt.xticks(np.arange(len(gps_list[n*batch:])), gps_list[n*batch:], rotation=45, horizontalalignment="right")
            plt.title(title)
            plt.savefig(os.path.join(save_path, plot_name + "_batch_" + str(n + 1) + ".png"), bbox_inches="tight", dpi=300)
            plt.close("all")
