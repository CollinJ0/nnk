#!/usr/bin/env python
# filename: variants.py

import inspect
import logging
import pandas as pd

import errors
import config

class Variants:
    #class_attribute
    purpose = config.VPURPOSE

    def __init__(self, twist_file, trimmed_ref):
        self.twist_file = twist_file
        self.trimmed_ref = trimmed_ref
        self.unique = None
        self.variants_map_to_ref = None

#         self.index=
        
#     def __iter__(self):
#         return self

#     def __next__(self):
#         if self.index == 0:
#             raise StopIteration
#         self.index = self.index - 1
#         return self.data[self.index]
    
    @property
    def variant_df(self):
        return Variants.twist_to_variant_df(self.twist_file)
    
    @property
    def variant_fasta_str(self):
        df = self.variant_df
        
        if self.unique:
            return ''.join([config.VFS1.format(self.trimmed_ref)]+
                           [config.VFS2.format(df.loc[i][config.VPOSCOL],
                                               df.loc[i][config.VVDNACOL]) for i in df.index])
    
    @property
    def unique(self):        
        return self._unique
            
    @unique.setter
    def unique(self, value):
        __vvdnacol=self.variant_df[config.VVDNACOL]
        __full = len(list(__vvdnacol))
        __unique = len(list(set(list(__vvdnacol))))
        __test=(__full==__unique)
        
        def get_test_str():
            __test_str = config.UNIQUE if __test else config.NUNIQUE
            return __test_str
        
        def get_my_name():
            ans = []
            frame = inspect.currentframe().f_back
            tmp = dict(list(frame.f_globals.items()) + list(frame.f_locals.items()))
            for k, var in tmp.items():
                if isinstance(var, self.__class__):
                    if hash(self) == hash(var):
                        ans.append(k)
            return ans[0]
        
        if __unique > __full:
            raise NameError(errors.EC1)
        
        if value == None:
            self._unique = __test
            
        elif value == __test:
            self._unique = value
        else:
            print((errors.RUNIQUE.format(str(value),
                                         get_test_str(),
                                         get_my_name(),
                                         str(__test))))
            self._unique= __test
            
        
        return self._unique
    
    @property
    def variants_map_to_ref(self):
        return self._variants_map_to_ref
    
    @variants_map_to_ref.setter
    def variants_map_to_ref(self, value):
        if value != None and type(value)==bool:
            self._variants_map_to_ref = value
        elif any(_v in self.trimmed_ref for _v in list(self.variant_df[config.VVDNACOL])):
            self._variants_map_to_ref = True
        else:
            self._variants_map_to_ref = False
        return self._variants_map_to_ref
    
    def variant_fasta_to_file(self, fname):
        with open(config.VFS3.format(fname), config.W) as f:
            f.write(self.variant_fasta_str)
    
    @staticmethod
    def twist_to_variant_df(tf):
        try:
            df = pd.read_excel(tf, sheet_name=config.VSNAME).drop(0)
            Variants.validate_df(df)
        except:
            raise NameError(errors.EC2)
        return df
    
    @staticmethod
    def validate_df(_df):
        correct_columns = config.VCCOL

        assert (list(_df.columns)==correct_columns),errors.EC3
        return (_df)

# TODO
# if variants map to ref, need to make new variants fasta
# if variants aren't unique, need to make variants unique