import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#
import config.ana_cff as cfg
from lib.fun_library import  getZhbbBaseCuts as tight_sel


def main():
    # need to get metaData and average weights (after baseline "tight") per year
    store_dict = {}
    y= '2016'
    for y in cfg.Years:
        tot_events, df, tight_events = getData(y)
        store_dict[y] = {'tot':tot_events,'tight':tight_events,'mean_w':df}
    def q_dict(key_, extra_key=None):
        if extra_key is None:
            return np.array([store_dict[y][key_] for y in cfg.Years])
        else:
            if extra_key == 'PrefireWeight':
                return np.array([store_dict[y][key_][extra_key] for y in ['2016','2017']]+[0])
            if extra_key == 'HEM_weight':
                return np.array([0,0]+[store_dict[y][key_][extra_key] for y in ['2018']])
            return np.array([store_dict[y][key_][extra_key] for y in ['2016','2017','2018']])
    print("Tot events (unweighted),              2016:{:12}, 2017:{:12}, 2018:{:12}".format(*q_dict('tot')))
    print("Events passing baseline (no tight),   2016:{:12}, 2017:{:12}, 2018:{:12}".format(*q_dict('tight')))
    print("Ratio of base/tot *100,               2016:{:12.5f}, 2017:{:12.5f}, 2018:{:12.5f}".format(*(100*(q_dict('tight')/q_dict('tot')))))
    for key in store_dict['2016']['mean_w'].keys(): # dummy
        print(f"Average for weight: {key:15}"+"   2016:{:12.5f}, 2017:{:12.5f}, 2018:{:12.5f}".format(*q_dict('mean_w',key)))


    

def getData(year):
    md = pd.read_pickle(cfg.postSkim_dir+f'/{year}/TTBar/TTToSemiLeptonic.pkl')['metaData']
    df = pd.read_pickle(cfg.master_file_path+f'{year}/mc_files/TTBar_val.pkl')
    #df = df[tight_sel(df)]
    df.loc[:,'genWeight'] = np.sign(df['genWeight'])
    df = df[(df['tt_type']=='Semi') & (df['tt_B'] == False)].filter(regex=r'\w*((?<!G)sf|eight)$')
    return md['tot_events'] - md['ttbb_events'], df.mean(), len(df)

if __name__ == '__main__':
    main()
