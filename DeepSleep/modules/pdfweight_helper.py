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
    #MATCH_VARS = ['run','luminosityBlock','event','LHE_HT']
    MATCH_VARS = ['run','luminosityBlock','event']

    pdf_sample_dict = {
        "TTbb_2L2Nu":"TTbb_4f_TTTo2l2nu_TuneCP5-Powheg-Openloops-Pythia8",
        "TTbb_SemiLeptonic":"TTbb_4f_TTToSemiLeptonic_TuneCP5-Powheg-Openloops-Pythia8",
        "TTbb_Hadronic":"TTbb_4f_TTToHadronic_TuneCP5-Powheg-Openloops-Pythia8",
        'ttHTobb':"ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8",
        'ttHToNonbb':"ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8",
    }

    pdf_sample_otf = ['TTZToQQ','TTZToBB','TTZToLLNuNu']

    def __init__(self, obj):
        self.pdfweights = None
        self.pdf_func = (lambda _ : None) # default do nothing
        if obj.jec_sys is None: 
            if obj.sample not in self.pdf_sample_dict.keys() and obj.sample not in self.pdf_sample_otf:
                pass
            elif obj.sample in self.pdf_sample_dict.keys():
                self.pdf_func = self.pdfweight_matcher
            else:# must be in pdf_sample_otf
                self.pdf_func = self.pdfweight_otf_calc
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
        main_tree_df = pd.DataFrame.from_dict({k: obj.tarray(k) for k in self.MATCH_VARS})
        # need to get pdf tree
        pdf_file = f'{cfg.pdfNtuples_dir}/{obj.year}/{self.pdf_sample_dict[obj.sample]}.pkl'
        pdf_dict = AnaDict.read_pickle(pdf_file)
        pdf_tree_df = pdf_dict.loc(self.MATCH_VARS).to_df()
        # merge
        pdf_tree_df['pdf_index'] = pdf_tree_df.index
        pdf_index_key = main_tree_df.merge(pdf_tree_df, how='left', on=self.MATCH_VARS)['pdf_index'].to_numpy()
        #print(main_tree_df[main_tree_df['luminosityBlock'] == 43714])
        #print(pdf_tree_df)
        #print(pdf_index_key)
        sorted_pdf_weights = pdf_dict['LHEPdfWeight'][pdf_index_key,:]
        self.pdfweights = sorted_pdf_weights
        
    def apply_cuts(self,obj):
        if self.pdfweights is not None:
            self.pdfweights = self.pdfweights[obj.precut][obj.event_mask][obj.btag_event_mask] # apply event level cuts

    def add_pdfweights_to_events(self,obj):
        if self.pdfweights is not None:
            n_weights = max(self.pdfweights.counts)
            for i in range(n_weights):
                obj.events[f'pdfweight_{i}'] = self.pdfweights[:,i]
