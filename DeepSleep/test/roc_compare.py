import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from lib.fun_library import getZhbbBaseCuts as zh_cuts
import config.ana_cff as cfg


def get_val_info(sample):
    s_dict = {}
    for y in cfg.Years:
        s_dict[y] = pd.read_pickle(f'files/{y}/mc_files/{sample}_val.pkl')
    return s_dict


def signal_wrapper(sig):
    for y in sig:
        x = sig[y].copy()
        sig[y] = x[((x['Zbb']==True) | (x['Hbb']==True))]
    return sig

def prep_for_roc(sig,bkg):
    sb_dict = {}
    for y in cfg.Years:
        bkg[y]['Y'] = 0
        sig[y]['Y'] = 1
        sb_dict[y] = pd.concat([bkg[y],sig[y]])
    return sb_dict

def plot_roc(sb):
    from sklearn import metrics
    fig, ax = plt.subplots()
    plt.plot([0,1],[0,1], 'k--')
    for y in cfg.Years:
        sb[y] = sb[y][zh_cuts(sb[y])]
        fpr, tpr, thresh = metrics.roc_curve(sb[y]['Y'], sb[y]['NN'])
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

def main():
    tt = get_val_info('TTBarLep')
    sig = signal_wrapper(get_val_info('TTZH'))
    plot_roc(prep_for_roc(sig,tt))

if __name__ == '__main__':
    main()
    
