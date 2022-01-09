import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from lib.fun_library import t2Run, save_pdf, getZhbbBaseCuts, getZhbbWeight, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel
import config.ana_cff as cfg

#nn = 'withbbvl_NN'
nn = cfg.nn

def get_val_info(sample):
    s_dict = {}
    for y in cfg.Years:
        s_dict[y] = pd.read_pickle(f'{sys.path[1]}/files/{y}/mc_files/{sample}_val.pkl')
    return s_dict


def signal_wrapper(sig):
    for y in sig:
        x = sig[y].copy()
        #sig[y] = x[((x['Zbb']==True) | (x['Hbb']==True))]
        sig[y] = x[x['matchedGen_ZHbb_bb'] == True]
    return sig

def prep_for_roc(sig,bkg):
    sb_dict = {}
    for y in cfg.Years:
        bkg[y]['Y'] = 0
        sig[y]['Y'] = 1
        sb_dict[y] = pd.concat([bkg[y],sig[y]])
    return sb_dict

#@save_pdf("roc_compary.pdf")
def plot_roc_by_year(sb):
    from sklearn import metrics
    fig, ax = plt.subplots()
    plt.plot([0,1],[0,1], 'k--')
    for y in cfg.Years:
        sb[y] = sb[y][getZhbbBaseCuts(sb[y])]
        fpr, tpr, thresh = metrics.roc_curve(sb[y]['Y'], sb[y][nn], sample_weight = sb[y]['weight'])
        plt.plot(fpr,tpr,
                 label=f'{y}, AUC: {metrics.auc(fpr,tpr):.4f}')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.ylim(0,1.05)
    plt.xlim(0,1)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.grid(True)
    plt.title('ROC Curve per Year')
    plt.legend()
    plt.show()
    plt.close()

def plot_roc_by_nn(sb, alt_name='TTBar'):
    from sklearn import metrics
    fig, ax = plt.subplots()
    ax.plot([0,1],[0,1], 'k--')
    # take 2018
    for nn_type, nn_name in zip([nn, 'binary_NN'],['multi-clss','binary']):
        sb['2018'] = sb['2018'][getZhbbBaseCuts(sb['2018'])]
        fpr, tpr, thresh = metrics.roc_curve(sb['2018']['Y'], sb['2018'][nn_type], sample_weight = sb['2018']['weight'])
        ax.plot(fpr,tpr,
                 label=f'{nn_name}, AUC: {metrics.auc(fpr,tpr):.4f}')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_ylim(0,1.05)
    ax.set_xlim(0,1)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True)
    fig.suptitle('ROC Curve, multi-class vs. binary'+f' (Sig vs {alt_name})')
    ax.legend()
    #plt.show()
    #plt.close()

@save_pdf("roc_comapre_nn.pdf")
def main():
    tt = get_val_info('TTBar')
    ttbb = get_val_info('ttbb')
    #bkg = {k:pd.concat([ttbb[k],tt[k]],axis='rows',ignore_index=True) for k in cfg.Years}
    ttz = signal_wrapper(get_val_info('ttZ'))
    tth = signal_wrapper(get_val_info('ttH'))
    sig = {k: pd.concat([ttz[k],tth[k]],axis='rows',ignore_index=True) for k in cfg.Years}
    #plot_roc_by_year(prep_for_roc(sig,tt))
    plot_roc_by_nn(prep_for_roc(sig,tt))
    plot_roc_by_nn(prep_for_roc(sig,ttbb), alt_name='tt+bb')

if __name__ == '__main__':
    main()
    
