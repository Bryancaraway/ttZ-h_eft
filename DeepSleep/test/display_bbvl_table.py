import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import numpy as np
import json
import config.ana_cff as cfg


def json_to_latex_table(year='2016'):
    jsf_name = cfg.dataDir+f'/deepak8sf/deepak8_bbvL_sf_{year}.json'
    jsf_dict = json.load(open(jsf_name,'r'))
    # 'score_bins', 'pt_bins', 'sdm_bins', 'real', 'fake'
    # ['score_pt_sdm_sf', 'score_pt_sdm_sf_err']
    score_bins, pt_bins, sdm_bins = [jsf_dict[k] for k in ['score_bins', 'pt_bins', 'sdm_bins']]
    fsf, fsf_err = np.array(jsf_dict['fake']['score_pt_sdm_sf']), np.array(jsf_dict['fake']['score_pt_sdm_sf_err'])
    ##
    print('\hline\hline')
    print(r'\pt bin & '+' & '.join(
        [f'${sdm_bins[j]}\ge\sdm<{sdm_bins[j+1]}$' for j in range(len(pt_bins[:-1])-1)]) + r'\\')
    for i in range(1,len(score_bins)-1):
        print('\hline') 
        print('\multicolumn{5}{c}{bbvL score '+f'$({score_bins[i]:.2f},{score_bins[i+1]:.2f}]$'+r'}\\')
        #print('\hline')
        for j in range(len(pt_bins[:-1])-1):
            print(
                (
                    f'$({pt_bins[j]},{pt_bins[j+1]}]$ & ' if j != len(pt_bins[:-1])-2 else f'$[{pt_bins[j]},\infty)$ & '
                )+' & '.join(
                    #[f'${fsf[i,j,k]:.2f}^{{{fsf_err[i,j,k]:.2f}}}_{{{fsf_err[i,j,k]:.2f}}}$' if not (j == 0 and k == 3) else '- ' for k in range(len(pt_bins[:-1])-1)]
                    [f'${fsf[i,j,k]:.2f}^{{{fsf_err[i,j,k]:.2f}}}_{{{fsf_err[i,j,k]:.2f}}}$'  for k in range(len(pt_bins[:-1])-1)]
                ) + r' \\'
            )


                
            

def main():
    for y in cfg.Years:
        print(f'\n{y}\n')
        json_to_latex_table(y)

if __name__ == '__main__':
    main()
