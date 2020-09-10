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


def export1d(histo, name):
    '''
    must pass a dict with key 'sumw'
    and optionally sumw2
    and must pass edges and name of hist
    '''
    label = name
    sumw = histo['sumw']
    edges = np.linspace(0, len(sumw), len(sumw)+1)

    out = TH1.__new__(TH1)
    out._fXaxis = TAxis(len(edges) - 1, edges[0], edges[-1])
    out._fXaxis._fName = name
    out._fXaxis._fTitle = label
    if len(set(edges[1:] - edges[:-1])) > 1: # means bin_w are not all the same
        out._fXaxis._fXbins = edges.astype(">f8")
    
    centers = (edges[1:] + edges[:-1]) / 2.0
    out._fEntries = out._fTsumw = out._fTsumw2 = sumw.sum() # no overflow or underflow here so need to do total sum
    out._fTsumwx = (sumw * centers).sum()#  might need  [1:-1]
    out._fTsumwx2 = (sumw * centers**2).sum()
    
    out._fName  = "histogram"
    out._fTitle =  label
    
    out._classname = b"TH1D"
    out.extend(sumw.astype(">f8"))
    if 'sumw2' in histo:
        out._fSumw2 = histo['sumw2'].astype(">f8")
    else:
        out._fSumw2 = sumw.astype(">f8")
    
    return out
