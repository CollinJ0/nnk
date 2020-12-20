#!/usr/bin/env python

import pandas as pd

class variants:
    #class_attribute
    purpose = 'unknown'

    def __init__(self, twist_file, trimmed_ref):
        self.twist_file = twist_file
        self.trimmed_ref = trimmed_ref

    def variant_df(self):
        return twist_to_variant_df(self.twist_file)

    def variant_fasta_str (self):
        return ''.join(['>{}\n{}\n'.format()])


def twist_to_variant_df(tf):
    try:
        df = pd.read_excel(tf, sheet_name='Variant Sheet').drop(0)
        validate_df(df)
    except:
        raise NameError('Input Twist Bio Excel file is incorrectly formatted')
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
