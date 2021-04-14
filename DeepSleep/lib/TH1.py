'''
Cutstom TH1 Class 
that converts dictionary
with sumw and sumw2 keys
and with edges
into uproot_methods TH1 object
'''

from uproot_methods.classes.TH1 import Methods as TH1Methods
import numpy as np

class TH1(TH1Methods, list) :
    pass

class TAxis(object):
    def __init__(self, fNbins, fXmin, fXmax):
        self._fNbins = fNbins
        self._fXmin = fXmin
        self._fXmax = fXmax


def export1d(histo, name, z_to_e=False):
    '''
    must pass a dict with key 'sumw'
    and optionally sumw2
    and must pass edges and name of hist
    '''
    label = name
    if z_to_e: # protect against completely zero yield bins
        if 0 in histo['sumw']:
            histo['sumw']  = np.where(histo['sumw'] == 0., 0.000001,             histo['sumw'])
            histo['sumw2'] = np.where(histo['sumw'] == 0., np.power(0.000001,2), histo['sumw2'])
    sumw = np.clip(np.pad(np.nan_to_num(histo['sumw']), 1, 'constant', constant_values=0), 0.,np.inf)
    if 'sumw2' in histo:
        sumw2 = np.pad(np.nan_to_num(histo['sumw2']), 1, 'constant', constant_values=0).astype(">f8")
        if 'Up' not in name and 'Down' not in name:
            sumw2 = np.where( np.sqrt(sumw2) > sumw , np.power(sumw,2) , sumw2 )
    else:
        sumw2 = sumw.astype(">f8")
    edges = np.linspace(0, len(sumw[1:-1]), len(sumw[1:-1])+1)

    out = TH1.__new__(TH1)
    out._fXaxis = TAxis(len(edges) - 1, edges[0], edges[-1])
    out._fXaxis._fName = name
    out._fXaxis._fTitle = label
    if len(set(edges[1:] - edges[:-1])) > 1: # means bin_w are not all the same
        out._fXaxis._fXbins = edges.astype(">f8")
    
    centers = (edges[1:] + edges[:-1]) / 2.0
    out._fEntries = out._fTsumw = out._fTsumw2 = sumw[1:-1].sum() 
    out._fTsumwx = (sumw[1:-1] * centers).sum()#  might need  [1:-1] for underflow and overflow, I do...
    out._fTsumwx2 = (sumw[1:-1] * centers**2).sum()
    
    out._fName  = "histogram"
    out._fTitle =  label
    
    out._classname = b"TH1D"
    out.extend(sumw.astype(">f8"))
    out._fSumw2 = sumw2

    return out
