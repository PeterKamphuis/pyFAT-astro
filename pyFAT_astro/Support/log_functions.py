# -*- coding: future_fstrings -*-
# This module  contains (eventually all) functions that relate to the 
# logging process
import os
import sys
import traceback

from datetime import datetime
from inspect import stack


def enter_recovery_point(Configuration, Fits_Files = None, Tirific_Template= None,
                         message = 'You should have specified a message',point_ID='Random'):
    
    update_statistics(Configuration, message= message)
    print_log(f'''ENTER_RECOVERY_POINT: Creating a recovery point.
{message}              
''',Configuration,case= ['debug_start','main'])
    write_config(f'{Configuration["LOG_DIRECTORY"]}CFG_{point_ID}.txt', Configuration, Fits_Files=Fits_Files,\
                 Tirific_Template=Tirific_Template,message=message )
def get_usage_statistics(Configuration,process):
        # psutil returns bytes
    memory_in_mb = (process.memory_info()[0])/2**20. 
    cpu_percent = process.cpu_percent(interval=1)
    return cpu_percent,memory_in_mb

get_usage_statistics.__doc__ =f'''
 NAME:
    get_usage_statistics
 PURPOSE:
    use psutil to get the current CPU and memory usage of tirific

 CATEGORY:
    log_functions 

 INPUTS:
    Configuration = Standard FAT configuration
    process_id = process id of the tirific currently running.

 OPTIONAL INPUTS:


 OUTPUTS:
    CPU = The current CPU usage
    mem = current memory usage in Mb

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    pyFAT version < 1.0.0 uses top which only works on unix and has an 
    error in the MB calculation.
'''



# A simple function to return the line numbers in the stack from where
# the functions are called.
def linenumber(debug='short'):
    '''get the line number of the print statement in the main.'''
    line = []
    for key in stack():
        if key[1] == 'main.py':
            break
        if key[3] != 'linenumber' and key[3] != 'print_log' and key[3] != '<module>':
            file = key[1].split('/')
            to_add= f"In the function {key[3]} at line {key[2]}"
            if debug == 'long':
                to_add = f"{to_add} in file {file[-1]}."
            else:
                to_add = f"{to_add}."
            line.append(to_add)
    if len(line) > 0:
        if debug == 'long':
            line = ', '.join(line)+f'\n{"":8s}'
        elif debug == 'short':
            line = line[0]+f'\n{"":8s}'
        else:
            line = f'{"":8s}'
    else:
        for key in stack():
            if key[1] == 'main.py':
                line = f"{'('+str(key[2])+')':8s}"
                break
    return line

linenumber.__doc__ =f'''
 NAME:
    linenumber

 PURPOSE:
    get the line number of the print statement in the main. Not sure 
    how well this is currently working.

 CATEGORY:
    log_functions 

 INPUTS:

 OPTIONAL INPUTS:


 OUTPUTS:
    the line number of the print statement

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    If debug = True the full stack of the line print will be given, 
    in principle the first debug message in every function should set 
    this to true and later messages not.
    !!!!Not sure whether currently the linenumber is produced due to 
    the restructuring.
'''

def print_log(log_statement,Configuration, case = None):
   if case is None:
      case=['main']
   
   debugging = False
   debug= 'empty'
   if Configuration['DEBUG']:
      trig=False
      if Configuration['DEBUG_FUNCTION'] == 'ALL':
         trig = True
      else:
         # get the function  
         for key in stack():
            if key[3] != 'linenumber' and key[3] != 'print_log' and key[3] != '<module>': 
               current_function= f"{key[3]}"
               break
         if current_function.lower() in [x.lower() for x in  Configuration['DEBUG_FUNCTION']]:
            trig=True             
      if trig:
         debugging=True    
         if 'debug_start' in case:
            debug = 'long'
         else:
            debug= 'short'
   if Configuration['TIMING']:
      log_statement = f"{linenumber(debug=debug)} {datetime.now()} {log_statement}"
   else:
      log_statement = f"{linenumber(debug=debug)}{log_statement}"
   print_statement = False
   if (debugging and ('debug_start' in case or 'debug_add' in case))\
      or ('verbose' in case and (Configuration['VERBOSE_LOG'] or debugging))\
         or 'main' in case:
            print_statement = True
   if print_statement:
      if Configuration['VERBOSE_SCREEN'] \
         or not Configuration['OUTPUTLOG']  \
         or 'screen' in case \
         or (Configuration['VERBOSE_LOG'] and 'main' in case):
         print(log_statement)
      if Configuration['OUTPUTLOG']:
         with open(Configuration['OUTPUTLOG'],'a') as log_file:
            log_file.write(log_statement)

print_log.__doc__ =f'''
 NAME:
    print_log
 PURPOSE:
    Print statements to log if existent and screen if Requested
 CATEGORY:
    log_functions 

 INPUTS:
    log_statement = statement to be printed
    Configuration = Standard FAT Configuration

 OPTIONAL INPUTS:


    screen = False
    also print the statement to the screen

 OUTPUTS:
    line in the log or on the screen

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    linenumber, .write

 NOTE:
    If the log is None messages are printed to the screen.
    This is useful for testing functions.
'''

def update_statistics(Configuration,process= None,message = None ):
    if Configuration['TIMING']:
        function = traceback.format_stack()[-2].split('\n')
        function = function[0].split()[-1].strip()
        if not process:
            process = Configuration['FAT_PSUPROCESS']
            program = 'pyFAT'
        else:
            program = 'TiRiFiC'

        CPU,mem = get_usage_statistics(Configuration,process)
        with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'a') as file:
            if message:
                file.write(f"# {function.upper()}: {message} at {datetime.now()} \n")
            # We cannot copy the process so initialize in the configuration
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb for {program} \n")
update_statistics.__doc__ =f'''
 NAME:
    update_statistics
 PURPOSE:
    Update the files that track the usage statistics 
 CATEGORY:
    log_functions

 INPUTS:
    Configuration = Standard FAT Configuration

 OPTIONAL INPUTS:
    process = the process ID to track, if not sset the FAT_PSUPROCESS 
              in the Configuration will tracked.
    
    message = Message that will be printed at line before the usages 
              statistics are printed
    

 OUTPUTS:
    update statistic file   
 
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
   
 NOTE:
    If the message is not provided  the CPU and Memory will be dumped 
    without an indication of the process that has been tracked.
'''


def write_config(file,Configuration, Fits_Files = None, Tirific_Template = None,message=None ):
    #be clear we are pickle dumping
    tmp = os.path.splitext(file)
    file = f'{tmp[0]}_RP_{Configuration["RP_COUNTER"]}.pkl'
    print_log(f'''WRITE_CONFIG: writing the configuration to {file}
''',Configuration,case=['debug_start'])
    # Separate the keyword names
    #Proper dictionaries are not pickable

    Pick_Configuration = {}
    for key in Configuration:
        Pick_Configuration[key] = Configuration[key]
        if key == 'FAT_PSUPROCESS':
            Pick_Configuration[key] = None
    
    if Fits_Files != None:
         Pick_Configuration['Fits_Files'] = Fits_Files
    if Tirific_Template != None:
         Pick_Configuration['Tirific_Template'] = Tirific_Template
    if message != None:
         Pick_Configuration['Message'] = message
    
    pythonv = sys.version_info
    if pythonv[0] > 3 or (pythonv[0] == 3 and pythonv[1] > 6.): 
        import pickle
        with open(file,'wb') as tmp:
            pickle.dump(Pick_Configuration,tmp) 
            Configuration['RP_COUNTER'] +=1
    else:
        print_log(f'''WRITE_CONFIG: This function is unfortunately not surported for python < 3.6 
''',Configuration,case=['main'])
  
    


write_config.__doc__ =f'''
 NAME:
    write_config

 PURPOSE:
    Write a config file to the fitting directory.

 CATEGORY:
    write_functions

 INPUTS:
    file = name of the file to write to
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


 OUTPUTS:
    A FAT config file.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: This doesn't work in python 3.6
This Could be used to pickup the fitting from where we left of

'''
