#!/usr/bin/env python
# filename: variants.py

import inspect
import logging
import pandas as pd

class Variants:
    #class_attribute
    purpose = 'unknown'

    def __init__(self, twist_file, trimmed_ref):
        self.twist_file = twist_file
        self.trimmed_ref = trimmed_ref
        self.unique = None
        self.variants_map_to_ref = None
    
    @property
    def variant_df(self):
        return twist_to_variant_df(self.twist_file)
    
    @property
    def variant_fasta_str(self):
        df = self.variant_df
        
        if self.unique:
            return ''.join(['>Ref\n{}\n'.format(self.trimmed_ref)]+['>{}\n{}\n'.format(df.loc[i]['Position Name (Optional)'], df.loc[i]['Variant DNA']) for i in df.index])
    
    @property
    def unique(self):        
        return self._unique
            
    @unique.setter
    def unique(self, value):
        __full = len(list(self.variant_df['Variant DNA']))
        __unique = len(list(set(list(self.variant_df['Variant DNA']))))
        __test=(__full==__unique)
        
        if __unique > __full:
            raise NameError('Input Twist Bio Excel file is incorrectly formatted: Error code 1')
        
        if value == None:
            self._unique = __test
            
        elif value == __test:
            self._unique = value
        else:
            __test_str = 'Unique' if __test else 'Not Unique'
            print('No, you are trying to assign {} to a variants object that is {}.\nVariants property {}.unique remains {}'.format(str(value), 
                                                                                                                                    __test_str, 
                                                                                                                                    self.get_my_name(), 
                                                                                                                                    str(__test)))
            self._unique= __test
            
        
        return self._unique
    
    @property
    def variants_map_to_ref(self):
        return self._variants_map_to_ref
    
    @variants_map_to_ref.setter
    def variants_map_to_ref(self, value):
        if value != None and type(value)==bool:
            self._variants_map_to_ref = value
        elif any(_v in self.trimmed_ref for _v in list(v1.variant_df['Variant DNA'])):
            self._variants_map_to_ref = True
        else:
            self._variants_map_to_ref = False
        
        return self._variants_map_to_ref
        
    def get_my_name(self):
        ans = []
        frame = inspect.currentframe().f_back
        tmp = dict(list(frame.f_globals.items()) + list(frame.f_locals.items()))
        for k, var in tmp.items():
            if isinstance(var, self.__class__):
                if hash(self) == hash(var):
                    ans.append(k)
        return ans[0]

    def variant_fasta_to_file(self, fname):
        with open('./{}'.format(fname), 'w') as f:
            f.write(self.variant_fasta_str)

def twist_to_variant_df(tf):
    try:
        df = pd.read_excel(tf, sheet_name='Variant Sheet').drop(0)
        validate_df(df)
    except:
        raise NameError('Input Twist Bio Excel file is incorrectly formatted: Error code 2')
    return df

def validate_df(_df):
    correct_columns = ['Group Position Name',
    'Position Name (Optional)',
    'DNA start index',
    'DNA end index',
    'Upstream sequence (at least 10bp)',
    'Downstream sequence (at least 10bp)',
    'WT DNA',
    'Variant DNA']
    
    assert (list(_df.columns)==correct_columns),"Input Twist Bio Excel file is incorrectly formatted"
    return (_df)

# TODO
# if variants map to ref, need to make new variants fasta
# if variants aren't unique, need to make variants unique
