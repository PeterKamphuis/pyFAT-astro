#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
from inspect import getframeinfo,stack
import os
import traceback



# A simple function to return the line numbers in the stack from where the functions are called
def linenumber():
    for key in stack():
        if key[1] == 'FAT.py':
            return f'({str(key[2])})'
    return f'Failed'
# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
class Proper_Dictionary(OrderedDict):
    def insert(self, existing_key, new_key, key_value):
        done = False
        if new_key in self:
            self[new_key] = key_value
            done = True
        else:
            new_orderded_dict = self.__class__()
            for key, value in self.items():
                new_orderded_dict[key] = value
                if key == existing_key:
                    new_orderded_dict[new_key] = key_value
                    done = True
            if not done:
                new_orderded_dict[new_key] = key_value
                done = True
                print(
                    "----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")

#batch convert types
def convert_type(array, type = 'float'):
    if type =='int':
        print(array)
        return  [int(x) for x in array]
    elif type =='str':
        return  [str(x) for x in array]
    else:
        return  [float(x) for x in array]




def print_log(log_statement,log, screen = False):
    log_statement = f"{linenumber():8s}{log_statement}"
    if screen or not log:
        print(log_statement)
    if log:
        log_file = open(log,'a')
        log_file.write(log_statement)
        log_file.close()

print_log.__doc__ = '''
;+
; NAME:
;       print_log
;
; PURPOSE:
;       Print statements to log if existent and screen if Requested
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       cleanup(Configuration)
;
; INPUTS:
;     Configuration = Structure that has the current fitting directory and expected steps
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;       os
;
; EXAMPLE:
;

'''
