import uproot
import os
import sys
sys.path.insert(1,'/home/bcaraway/ttZh_ana/DeepSleep')
from modules.AnaDict import AnaDict
import numpy as np


years     = ['2016','2017','2018']
el_id_sel = 'Tight'
mu_id_sel = 'Medium'

def main(lep_type):

    for year in years:
        i_files = [f'Electron_SUSYScaleFactors_2017v2ID_Run{year}.root',
                   f'Electron_GT20GeV_RecoSF_2017v2ID_Run{year}.root',
                   f'Electron_GT10GeV_RecoSF_2017v2ID_Run{year}.root',
                   f'Muon_IDScaleFactor_wSys_{year}GH.root',
                   f'Muon_LooseID_MiniIso0p2SF_{year}.root',
                   f'Muon_IDScaleFactor_wSys_{year}.root',
                   f"Muon_{mu_id_sel}ID_MiniIso0p2SF_{'2017' if year == '2018' else year}.root"]
        h_list  = ['EGamma_SF2D',
                   f'Run{year}_CutBased{el_id_sel}NoIso94XV2',
                   f'Run{year}_Mini',
                   f'Run{year}_MVAVLooseTightIP2DMini',
                   f'NUM_{mu_id_sel}ID_DEN_genTracks_eta_pt',
                   'SF',
                   f'NUM_{mu_id_sel}ID_DEN_genTracks_pt_abseta',
                   f'TnP_MC_NUM_MiniIso02Cut_DEN_{mu_id_sel}ID_PAR_pt_eta',
                   f'NUM_{mu_id_sel}ID_DEN_TrackerMuons_pt_abseta']

        sf = AnaDict()
        for roofile in i_files:
            if os.path.exists(roofile) and lep_type[1:] in roofile:
                with uproot.open(roofile) as f_:
                    for h in h_list:
                        if h in f_:
                            print(roofile, h)
                            sf[h] = {}
                            hist = f_[h]
                            if ((lep_type == 'muon' and year != '2016') or h == 'SF'):
                                pt_bins, eta_bins = hist.edges
                                v    = np.array(hist.values).T
                                verr = np.array(hist.variances).T
                            else:    
                                eta_bins, pt_bins = hist.edges
                                v    = hist.values
                                verr = hist.variances

                            for i in range(v.shape[0]):
                                sf[h][f'{eta_bins[i]},{eta_bins[i+1]}'] = {f'{pt_bins[j]},{pt_bins[j+1]}': 
                                                                           {'values': v[i,j],
                                                                            'up': v[i,j]+verr[i,j],
                                                                            'down': v[i,j]-verr[i,j]}
                                                                           for j in range(v.shape[1])}
                            #
                        #
                    #
                #
            #
        #
        sfdir    = f'/home/bcaraway/ttZh_ana/DeepSleep/files/{year}/lepton_sf_files/'
        sf_out_file = f'{lep_type}_sf_{year}.pkl'
        sf.to_pickle(sfdir+sf_out_file)

if __name__ == '__main__':
    #main('electron')
    main('muon')
                        
