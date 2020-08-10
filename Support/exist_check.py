#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to check that output/input is present



class FileNotFoundError(Exception):
    pass
from support_functions import print_log
import os

def sofia_output(Configuration,Fits_Files, debug = False):

    req_files= ['MOMENT1','MOMENT0','MOMENT2','MASK']
    for file in req_files:
        if os.path.exists(Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files[file]):
            continue
        else:
            log_statement = f"CHECK_SOFIA_OUTPUT: The file {Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files[file]} is not found."
            print_log(log_statement, Configuration['OUTPUTLOG'])
            raise FileNotFoundError(log_statement)

    if not os.path.exists(Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'):
        log_statement = f"CHECK_SOFIA_OUTPUT: The file {Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'} is not found."
        print_log(log_statement, Configuration['OUTPUTLOG'])
        raise FileNotFoundError(log_statement)


sofia_output.__doc__ =f'''
Simple function to make sure all sofia output is present as expeceted
'''
