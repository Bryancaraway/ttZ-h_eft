## helper function to retrieve and handle pdf weight ##
import config.ana_cff as cfg
from lib.fun_library import t2Run
from modules.AnaDict import AnaDict
#
import lhapdf
import numpy as np
import pandas as pd
from numba import jit
import re

## to be used on the fly by Skim.py

class PDFHelper():

    PDFVARS = ['Generator_scalePDF', 'Generator_x1', 'Generator_x2', 'Generator_id1', 'Generator_id2']
    PDFSET  = 'NNPDF31_nnlo_hessian_pdfas' 
    #MATCH_VARS = ['run','luminosityBlock','event','LHE_HT', 'LHE_HTIncoming']
    MATCH_VARS = ['luminosityBlock','event']

    pdf_sample_dict = {
        "TTbb_2L2Nu":"TTbb_4f_TTTo2l2nu_TuneCP5-Powheg-Openloops-Pythia8",
        "TTbb_SemiLeptonic":"TTbb_4f_TTToSemiLeptonic_TuneCP5-Powheg-Openloops-Pythia8",
        "TTbb_Hadronic":"TTbb_4f_TTToHadronic_TuneCP5-Powheg-Openloops-Pythia8",
        #'ttHTobb':"ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8",
        #'ttHToNonbb':"ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8",
    }

    pdf_sample_otf = ['TTZToQQ','TTZToBB','TTZToLLNuNu']

    pdf_sample_ken = ['TTToSemiLeptonic','TTToHadronic','TTTo2L2Nu','ttHToNonbb','ttHTobb']

    def __init__(self, obj):
        self.pdfweights = None
        self.pdf_func = (lambda _ : None) # default do nothing, which is the current case for EFT
        if obj.jec_sys is None and not obj.isData: 
            if obj.sample not in list(self.pdf_sample_dict.keys())+self.pdf_sample_otf+self.pdf_sample_ken:
                pass
            elif obj.sample in self.pdf_sample_dict.keys():
                self.pdf_func = self.pdfweight_matcher
            elif obj.sample in self.pdf_sample_otf:
                self.pdf_func = self.pdfweight_otf_calc
            elif obj.sample in self.pdf_sample_ken:
                #pass # for now
                self.pdf_func = self.pdfweight_matcher_v2
        self.pdf_func(obj) # computes pdfweights if necessary

    @t2Run
    def pdfweight_otf_calc(self,obj): # wont work at the moment, need to do this inclusively, preserve cuts
        pdf_info = AnaDict({v: obj.tarray(v) for v in self.PDFVARS})
        pset = lhapdf.getPDFSet(self.PDFSET)
        pdfs = pset.mkPDFs()
        # get central weight
        n_events = pdf_info.size
        # following procedure
        x1_vals = np.array([[pdfs[j].xfxQ(pdf_info['Generator_id1'][i], pdf_info['Generator_x1'][i], pdf_info['Generator_scalePDF'][i]) for j in range(pset.size)] for i in range(n_events)])
        x2_vals = np.array([[pdfs[j].xfxQ(pdf_info['Generator_id2'][i], pdf_info['Generator_x2'][i], pdf_info['Generator_scalePDF'][i]) for j in range(pset.size)] for i in range(n_events)])
        #
        pdfweights = (x1_vals * x2_vals) / (x1_vals[:,0] * x2_vals[:,0])[:,None]
        from awkward import JaggedArray as aj
        # need to convert to jagged array
        self.pdfweights = aj.fromcounts(pdfweights.shape[-1]*np.ones(len(pdfweights), dtype=int), pdfweights.flatten())
        #self.pdfweights = ak.fromiter(pdfweights)
    
    @t2Run
    def pdfweight_matcher(self,obj):
        isDittbb = obj.sample == 'TTbb_2L2Nu' and obj.year == '2017'#  ttbb 2L2Nu 2017 needs extra care
        match_vars = self.MATCH_VARS if not isDittbb else self.MATCH_VARS + ['LHE_HT','LHE_HTIncoming'] 
        do_round  = (lambda _df : _df ) if not isDittbb else (lambda _df : _df.round({'LHE_HT':2,'LHE_HTIncoming':1}) )
        main_tree_df = do_round(pd.DataFrame.from_dict({k: obj.tarray(k) for k in match_vars}))
        # need to get pdf tree
        pdf_file = f'{cfg.pdfNtuples_dir}/{obj.year}/{self.pdf_sample_dict[obj.sample]}.pkl'
        pdf_dict = AnaDict.read_pickle(pdf_file)
        pdf_tree_df = do_round(pdf_dict.loc(match_vars).to_df())
        for var in match_vars:
            del pdf_dict[var]
        # merge
        pdf_tree_df['pdf_index'] = pdf_tree_df.index
        pdf_index_key = main_tree_df.merge(pdf_tree_df, how='left', on=match_vars)['pdf_index'].to_numpy()
        del pdf_tree_df , main_tree_df # trying to clear up memory
        #
        sorted_pdf_weights = pdf_dict['LHEPdfWeight'][pdf_index_key,:]
        #print(sorted_pdf_weights)
        self.pdfweights = sorted_pdf_weights

    @t2Run
    def pdfweight_matcher_v2(self,obj): # opens friend trees and extracts weight stored in Kens file area
        pdf_file = obj.roofile.replace('bcaraway','hatake').replace('.root','_fix.root')
        pdf_branch = "LHEPdfWeightf"
        if obj.sample == 'ttHToNonbb' and obj.year == '2016':
            pdf_file = re.sub(r'Skim_\d*_fix.root','friend.root', pdf_file)
            pdf_branch = "LHEPdfWeight"
        import uproot
        with uproot.open(pdf_file) as p_file:
            tree =  p_file.get("Events")
            pdfweights = tree.array(pdf_branch)
            if type(pdfweights) == np.ndarray: # need to convert to jagged array
                from awkward import JaggedArray as aj
                pdfweights = aj.fromcounts(pdfweights.shape[-1]*np.ones(len(pdfweights), dtype=int), pdfweights.flatten())
            self.pdfweights = pdfweights # LHEPdfWeight
            print(self.pdfweights)
            print(type(self.pdfweights))
        
    def apply_cuts(self,obj):
        if self.pdfweights is not None:
            self.pdfweights = self.pdfweights[obj.precut][obj.event_mask][obj.btag_event_mask] # apply event level cuts

    def add_pdfweights_to_events(self,obj):
        if self.pdfweights is not None:
            n_weights = max(self.pdfweights.counts)
            for i in range(n_weights):
                obj.events[f'pdfweight_{i}'] = self.pdfweights[:,i]
