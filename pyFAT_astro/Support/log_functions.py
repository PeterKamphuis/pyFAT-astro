# -*- coding: future_fstrings -*-
# This module  contains (eventually all) functions that relate to the 
# logging process
import numpy as np
import os
import psutil as psu
import sys
import time
import traceback
import warnings

from datetime import datetime
from inspect import stack
from pyFAT_astro.Support.fat_errors import ProgramError

with warnings.catch_warnings():
   warnings.simplefilter("ignore")
   import matplotlib
   matplotlib.use('pdf')
   import matplotlib.pyplot as plt
   import matplotlib.font_manager as mpl_fm


class full_system_tracking:
   def __init__(self,Configuration):
      self.stop = False
      self.pid = os.getpid()
      self.main_pyFAT = psu.Process(self.pid)
      self.user = self.main_pyFAT.username()
      self.python = self.main_pyFAT.name()
      self.tirific = Configuration['TIRIFIC']
      self.font_file = Configuration['FONT_FILE']
      self.sofia = Configuration['SOFIA2']
      self.file = f"{Configuration['MAIN_DIRECTORY']}FAT_Resources_Used.txt"
      self.plot_name= f"{Configuration['MAIN_DIRECTORY']}pyFAT_Resources_Monitor.pdf"
      self.cpus= psu.cpu_count()
      with open(self.file,'w') as resources:
         resources.write("# This file contains an estimate of all resources used for a pyFAT run. \n")
         resources.write(f"# {'Time':20s} {'Sys CPU':>10s} {'Sys RAM':>10s} {'FAT CPU':>10s} {'FAT RAM':>10s} \n")
         resources.write(f"# {'YYYY-MM-DD hh:mm:ss':20s} {'%':>10s} {'Gb':>10s} {'%':>10s} {'Gb':>10s} \n")
      self.interval = 30 # amount of second when to do new monitor

   def start_monitoring(self):
      while not self.stop:
         try:
               self.sys_cpu= psu.cpu_percent(interval=1)
               self.sys_ram= psu.virtual_memory().used/2**30.
               self.CPU = 0.
               self.RAM = 0.
               for proc in psu.process_iter():
                  if proc.username() == self.user \
                     and proc.status() == 'running'\
                     and (proc.name() == self.python or\
                        proc.name() == self.tirific or\
                        proc.name() == self.sofia or\
                        proc.name() == 'python3'):
                     try:
                           self.CPU += proc.cpu_percent(interval=0.5)/self.cpus
                           self.RAM += (proc.memory_info()[0])/2**30.
                     except:
                           pass
               #file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Gb for TiRiFiC \n")
               with open(self.file,'a') as resources:
                  resources.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S'):20s} {self.sys_cpu:>10.1f} {self.sys_ram:>10.2f} {self.CPU:>10.1f} {self.RAM:>10.2f} \n")
         except Exception as e:
               #We do not care if something goes wrong once. We don't want the monitor to crash
               #but we would like to know what went wrong
               traceback.print_exception(type(e),e,e.__traceback__)
               pass
         time.sleep(self.interval)

   def stop_monitoring(self):
      self.interval = 0.1
      self.stop = True
      try:
         loads = {'Time':[]}
         keys=['SCPU','SRAM','FCPU','FRAM']
         for key in keys:
               loads[key]  = []
         with open(self.file) as file:
               lines = file.readlines()
         startdate = 0
         #load data from file into dictionary
         for line in lines:
               line = line.split()
               if line[0] == '#':
                  continue
               else:
                  date = extract_date(f"{line[0]} {line[1]}")
               if startdate == 0:
                  startdate = date
               diff = date - startdate
               time = diff.total_seconds()/(3600.)
               loads['Time'].append(time)
               for i,key in enumerate(keys):
                  loads[key].append(float(line[int(2+i)]))
         #Plot the parameters
         try:
               mpl_fm.fontManager.addfont(self.font_file)
               font_name = mpl_fm.FontProperties(fname=self.font_file).get_name()
         except FileNotFoundError:
               font_name = 'DejaVu Sans'
         
         labelfont = {'family': font_name,
                  'weight': 'normal',
                  'size': 4}
         fig, ax1 = plt.subplots(figsize = (8,6))
         fig.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.3, top = 0.7)
         ax1.plot(loads['Time'],loads['SRAM'],'b--',lw=0.5,alpha=0.5, label='System RAM')
         ax1.plot(loads['Time'],loads['FRAM'],'b-',lw=0.5,alpha=1.0, label='pyFAT RAM')
         ax1.set_ylim(0,np.max(np.array(loads['SRAM']+loads['FRAM'],dtype=float))*1.1)
         ax1.set_ylabel('RAM (Gb) ', color='b')
         ax1.tick_params(axis='y', labelcolor='b')
         ax1.set_xlabel('Run Duration (h)', color='k',zorder=5)
         ax2 = ax1.twinx()
         ax2.plot(loads['Time'],loads['SCPU'],'r--',lw=0.5,alpha=0.5, label='System CPU')
         ax2.plot(loads['Time'],loads['FCPU'],'r-',lw=0.5,alpha=1.0, label='pyFAT CPU')
         ax2.set_ylim(0,np.max(np.array(loads['SCPU']+loads['FCPU'],dtype=float))*1.1)
         ax2.set_ylabel('CPUs (%)',color='r')
         ax2.tick_params(axis='y', labelcolor='r')
         fig.savefig(self.plot_name)
         plt.close()
      except Exception as e:
         traceback.print_exception(type(e),e,e.__traceback__)
         print(f'We failed to write the statistics plot.')
         pass
full_system_tracking.__doc__ =f'''
NAME:
   full_system_tracking

PURPOSE:
   A class that can be started in a thread and run in the background to track 
   all instances of tirific sofia 2 and the main python id   

 CATEGORY:
    log_functions 

 INPUTS:
   Configuration = Standard FAT configuration
   

 OPTIONAL INPUTS:


 OUTPUTS:


 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 
'''


def enter_recovery_point(Configuration, Fits_Files = None, Tirific_Template= None,
                         message = 'You should have specified a message',point_ID='Random'):
    
    update_statistics(Configuration, message= message)
    print_log(f'''ENTER_RECOVERY_POINT: Creating a recovery point.
{message}              
''',Configuration,case= ['debug_start','main'])
    write_config(f'{Configuration["LOG_DIRECTORY"]}CFG_{point_ID}.txt', Configuration, Fits_Files=Fits_Files,\
                 Tirific_Template=Tirific_Template,message=message )
enter_recovery_point.__doc__ =f'''
 NAME:
   enter_recovery_point
 PURPOSE:
   use pickle to enter a recovery point where the current configuration and other 
   tracking is written.

 CATEGORY:
    log_functions 

 INPUTS:
   Configuration = Standard FAT configuration
   Fits_Files = standard FAT Fits Files
   Tirific_Template = 
   message = message for the point
   point_ID = Name of Recovery point

 OPTIONAL INPUTS:


 OUTPUTS:


 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 
'''

def extract_date(string):
    tmp = string.split(' ')
    tmp2 = tmp[0].split('-')
    if len(tmp2) == 3:
        try:
            date =  datetime.strptime(f"{tmp[0]} {tmp[1]}", '%Y-%m-%d %H:%M:%S.%f')
        except ValueError:
            date =  datetime.strptime(f"{tmp[0]} {tmp[1]}", '%Y-%m-%d %H:%M:%S')
    else:
        raise ProgramError("There is no date in the provided string.")
    return date
extract_date.__doc__ =f'''
 NAME:
    extract_date

 PURPOSE:
    convert a string into a date object

 CATEGORY:
    write_functions

 INPUTS:
    string = string

 OPTIONAL INPUTS:


 OUTPUTS:
    date = the date object

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
    

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
      if 'ALL' in Configuration['DEBUG_FUNCTION']:
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
