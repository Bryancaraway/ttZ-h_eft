# Main data container for analysis 
# created by: Bryan Caraway
from collections import UserDict
import pickle5 as pickle


class AnaDict (UserDict): # may want to move this to different file
    ''' 
    pythonic dictionary with the ability to apply numpy indexing
    across all key values
    read/write to pickle with the same syntax as pandas --> may need to import speed
    '''
    # methods
    def loc(self,keys):
        return self({k: self.data[k] for k in keys})
    def pad(self, _i):
        # i is the amount to pad
        return self({k: self.data[k].pad(_i) for k in self.data})
    def fillna(self, _i):
        # i is the value to fill na
        return self({k: self.data[k].fillna(_i) for k in self.data})
    def sum(self):
        return sum(self.data.values())
    def values(self):
        return self.data.values()
    def to_pickle(self, out_name):
        with open(out_name, 'wb') as handle:
            pickle.dump(self.data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # magic methods
    def __getitem__(self,key):
        if type(key) is str:
            return self.data[key]
        else:
            return self({k: v[key] for k,v in self.data.items()}) 
    #class methods
    @classmethod
    def read_pickle(cls,input_name):
        print(input_name)
        with open(input_name, 'rb') as handle:
            return cls(pickle.load(handle))
    @classmethod
    def __call__(cls,dict_):
        return cls(dict_)
            
