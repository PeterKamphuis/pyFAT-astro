#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to write text files to Disk




def sofia(template,name):
    with open(name,'w') as file:
        for key in template:
            if key[0] == 'E' or key [0] == 'H':
                file.write(template[key])
            else:
                file.write(f"{key} = {template[key]}\n")




sofia.__doc__ = '''

; NAME:
;       sofia
;
; PURPOSE:
;       write a sofia2 dictionary into file
;
; CATEGORY:
;       write
;
;
; INPUTS:
;       name = Name of the sofia file
;       template = a template dictionary that matches the sof read template
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      split, strip, open
;
; EXAMPLE:
;
;
'''
