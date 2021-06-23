import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import os
from makeDatacard import Systematic, ShapeSystematic



@t2Run
def add_uncorrSystematics(obj):
    # SIG = ttZbb[0,1,2,3] ttHbb[0,1,2,3]
    # BKG = ttX, TTBar, old_tt_bb, VJets, other
    #process_line = np.array([obj.accepted_sig + obj.accepted_bkg for _ in range(obj.dc_bins)]).flatten()
    all_mc = obj.accepted_sig + obj.accepted_bkg
    all_but_ttbb = obj.accepted_sig + ['TTBar','ttX','VJets','other']
    tth_sig = re.findall(r'ttH\d', ' '.join(obj.accepted_sig))
    ttz_sig = re.findall(r'\w*ttZ\d', ' '.join(obj.accepted_sig))
    ttbar_mc = ['TTBar','tt_bb', 'tt_2b']
    #ttbar_mc = ['TTBar','old_tt_bb']
    #new_ttz_sig = ['new_'+z for z in ttz_sig]
    #
    #Systematic.set_dc_processes(obj.dc_dict, process_line)
    obj.write2dc(f'# Rate uncertainties\n')
    Systematic('lumi_2016', 'lnN',  all_mc, 1.025)
    Systematic('lumi_2017', 'lnN',  all_mc, 1.023)
    Systematic('lumi_2018', 'lnN',  all_mc, 1.025)
    #Systematic('tt2bxsec', 'lnN', ['old_tt_bb'], 1.5)
    # Shape Systatics
    obj.write2dc(100*'-'+'\n')
    obj.write2dc('# Shape uncertainties \n')
    #ShapeSystematic.set_df_histos_histfuncs(obj.data_dict, obj.histos)#, obj.hist3d, obj.ptclip)
    # when defining shape, must incluse whether it is a mcsta, scale, or up/down syst
    for y in obj.years:
        #obj.histos = ShapeSystematic(f'btg_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeight_Up','BTagWeight_Down').get_shape()
        # probably do xsec theo rates here
        # first signal
        Systematic(f'tth_ggpdf_{y}', 'lnN', tth_sig, 1.036)
        Systematic(f'ttz_ggpdf_{y}', 'lnN', tth_sig, 1.035)
        Systematic(f'tth_qsc_{y}' ,  'lnN', tth_sig, 1.058,0.908)
        Systematic(f'ttz_qsc_{y}'  , 'lnN', ttz_sig, 1.081,0.907) 
        # now background
        Systematic(f'ggpdf_{y}', 'lnN', ttbar_mc, 1.042)
        Systematic(f'qqpdf_{y}', 'lnN', ['ttX','VJets','other'],
                   [p_norms.rate_unc['pdf']['ttX'],
                    p_norms.rate_unc['pdf']['VJets'],
                    p_norms.rate_unc['pdf']['other']])# combine into one 
        #Systematic('qqpdf', 'lnN', ['VJets'], 1.038)# combine into one
        #Systematic('qqpdf', 'lnN', ['other'], 1.050)# combine into one
        #Systematic('tt_qsc'   , 'lnN', ttbar_mc+['ttX'], [[1.024,0.965] for _ in ttbar_mc]+[1.300])
        Systematic(f'tt_qsc_{y}'   , 'lnN', ttbar_mc, [[1.024,0.965] for _ in ttbar_mc])
        Systematic(f'ttx_qsc_{y}'  ,  'lnN', ['ttX'], p_norms.rate_unc['QCD_scale']['ttX']) 
        Systematic(f'v_qsc_{y}'    , 'lnN', ['VJets'],p_norms.rate_unc['QCD_scale']['VJets'])#1.008, 0.996) # .821/1.24
        Systematic(f'other_qsc_{y}', 'lnN', ['other'],p_norms.rate_unc['QCD_scale']['other']) #1.05, 0.95)   # .898/1.12
        Systematic(f'tt2bxsec_{y}', 'lnN', ['tt_2b'], 1.5)
        obj.histos = ShapeSystematic(f'btglf_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeightLight_Up','BTagWeightLight_Down').get_shape()
        obj.histos = ShapeSystematic(f'btghf_{y}', 'shape', 'up/down', all_mc, 1, 'BTagWeightHeavy_Up','BTagWeightHeavy_Down').get_shape()
        obj.histos = ShapeSystematic(f'lepsf_{y}', 'shape', 'up/down', all_mc, 1, 'lep_sf_up','lep_sf_down').get_shape()
        obj.histos = ShapeSystematic(f'trigeff_{y}', 'shape', 'up/down', all_mc, 1, 'lep_trig_eff_tight_pt_up','lep_trig_eff_tight_pt_down').get_shape()
        obj.histos = ShapeSystematic(f'pu_{y}',      'shape', 'up/down', all_mc, 1, 'puWeight_Up','puWeight_Down').get_shape()
        obj.histos = ShapeSystematic(f'ak4JER_{y}', 'shape', 'qconly', all_mc, 1,    extraQC=True).get_shape()
        obj.histos = ShapeSystematic(f'ak4JES_{y}', 'shape', 'qconly', all_mc, 1,    extraQC=True).get_shape()
        obj.histos = ShapeSystematic(f'ak8JER_{y}', 'shape', 'qconly', all_mc, 1,    extraQC=True).get_shape()
        obj.histos = ShapeSystematic(f'ak8JES_{y}', 'shape', 'qconly', all_mc, 1,    extraQC=True).get_shape() 
        
        #obj.histos = ShapeSystematic(f'toppt_{y}', 'shape', 'up/down', ttbar_mc, 1, 'Stop0l_topptWeight_Up' ,'Stop0l_topptWeight_Down').get_shape()
        obj.histos = ShapeSystematic(f'isr_{y}', 'shape', 'ps', ['TTBar'], 1, 'ISR_Up','ISR_Down').get_shape()
        obj.histos = ShapeSystematic(f'fsr_{y}', 'shape', 'ps', ['TTBar'], 1, 'FSR_Up','FSR_Down', extraQC=False).get_shape()
        #obj.histos = ShapeSystematic(f'mu_r', 'shape', 'normup/down', all_but_ttbb, 1, 'mu_r_Up','mu_r_Down').get_shape()
        #obj.histos = ShapeSystematic(f'mu_f', 'shape', 'normup/down', all_but_ttbb, 1, 'mu_f_Up','mu_f_Down').get_shape()
        obj.histos = ShapeSystematic(f'mu_r_tt_{y}', 'shape', 'normup/down', ['TTBar'], 1, 'mu_r_Up','mu_r_Down').get_shape()
        obj.histos = ShapeSystematic(f'mu_f_tt_{y}', 'shape', 'normup/down', ['TTBar'], 1, 'mu_f_Up','mu_f_Down').get_shape()
        obj.histos = ShapeSystematic(f'mu_r_tth_{y}', 'shape', 'normup/down', tth_sig, 1, 'mu_r_Up','mu_r_Down').get_shape()
        obj.histos = ShapeSystematic(f'mu_f_tth_{y}', 'shape', 'normup/down', tth_sig, 1, 'mu_f_Up','mu_f_Down').get_shape()
        obj.histos = ShapeSystematic(f'mu_r_ttz_{y}', 'shape', 'normup/down', ttz_sig, 1, 'mu_r_Up','mu_r_Down').get_shape()
        obj.histos = ShapeSystematic(f'mu_f_ttz_{y}', 'shape', 'normup/down', ttz_sig, 1, 'mu_f_Up','mu_f_Down').get_shape()
        obj.histos = ShapeSystematic(f'isr_ttbb_{y}', 'shape', 'ps', ['tt_bb', 'tt_2b'], 1, 'ISR_Up','ISR_Down').get_shape()
        obj.histos = ShapeSystematic(f'fsr_ttbb_{y}', 'shape', 'ps', ['tt_bb', 'tt_2b'], 1, 'FSR_Up','FSR_Down', extraQC=False).get_shape()
        obj.histos = ShapeSystematic(f'mu_r_ttbb_{y}', 'shape', 'normup/down', ['tt_bb', 'tt_2b'], 1, 'mu_r_Up','mu_r_Down').get_shape()
        obj.histos = ShapeSystematic(f'mu_f_ttbb_{y}', 'shape', 'normup/down', ['tt_bb', 'tt_2b'], 1, 'mu_f_Up','mu_f_Down').get_shape()
        obj.histos = ShapeSystematic(f'pdf_{y}', 'shape', 'normup/down', tth_sig+ttz_sig+['TTBar'], 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        obj.histos = ShapeSystematic(f'pdf_ttbb_{y}', 'shape', 'normup/down', ['tt_bb','tt_2b'],                    1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
        obj.histos = ShapeSystematic(f'UE_{y}',     'shape', 'qconly', ['TTBar'], 1, extraQC=True).get_shape()
        #obj.histos = ShapeSystematic('erdOn', 'shape',  ['TTBar'], 1) # not working with just one shape at the moment
        obj.histos = ShapeSystematic(f'hdamp_{y}', 'shape', 'qconly', ['TTBar'], 1, extraQC=True).get_shape()
        obj.histos = ShapeSystematic(f'hdamp_ttbb_{y}', 'shape', 'qconly', ['tt_bb', 'tt_2b'], 1, extraQC=True).get_shape()
        Systematic(f'CMS_ttbbnorm_{y}', 'lnN', ['tt_bb','tt_2b'], 2)
    #
    obj.histos = ShapeSystematic(f'pref_2016', 'shape', 'up/down', all_mc, 1, 'PrefireWeight_Up' ,'PrefireWeight_Down').get_shape()
    obj.histos = ShapeSystematic(f'pref_2017', 'shape', 'up/down', all_mc, 1, 'PrefireWeight_Up' ,'PrefireWeight_Down').get_shape()
    # redundant #obj.histos = ShapeSystematic(f'mu_rf', 'shape', 'normup/down', all_mc, 1, 'mu_rf_Up','mu_rf_Down').get_shape()
    #
    #obj.histos = ShapeSystematic(f'pdf_ttz', 'shape', 'up/down', ttz_sig, 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()
    #obj.histos = ShapeSystematic(f'pdf', 'shape', 'normup/down', all_but_ttbb, 1, 'pdfWeight_Up','pdfWeight_Down').get_shape()

    
    # These shapes are already computed, just need to add to datacard

        
    #
    obj.write2dc(100*'-'+'\n')
    obj.write2dc('# Float tt_bb normalization\n') 
    #obj.write2dc('CMS_ttbbnorm rateParam * tt_*b 1 [-10,10]\n')
    #Systematic('CMS_ttbbnorm', 'lnN', ['tt_bb','tt_2b'], 10)
    #
    obj.write2dc(100*'-'+'\n')
    obj.write2dc('# MC Stats uncertainties\n') 
    #obj.histos = ShapeSystematic(f'mcstat','shape','mcstat', all_mc, 1).get_shape()
    obj.write2dc('#* autoMCStats 10 0  1\n') 
        
    #
