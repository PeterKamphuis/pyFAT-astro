
from omegaconf import OmegaConf
import sys
from pyFAT_astro.Support.fat_errors import InputError
from pyFAT_astro.config.defaults import defaults
import psutil
import os
from importlib.resources import files as import_pack_files

def process_input(argv):
    if '-v' in argv or '--version' in argv:
        #The version is always printed so only need to exit()
        sys.exit()

    if '-h' in argv or '--help' in argv:
        from pyFAT_astro.Support.log_functions import return_help_message
        print(return_help_message())
        sys.exit()
    #First import the defaults
    cfg = OmegaConf.structured(defaults)
    #if ncpu is the total available remove 1
    if cfg.ncpu == psutil.cpu_count():
        cfg.ncpu -= 1

    # read command line arguments anything list input should be set in brackets '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
    inputconf = OmegaConf.from_cli(argv)
    cfg_input = OmegaConf.merge(cfg,inputconf)
    # Print examples if requested
    if cfg_input.print_examples:
        no_cube = OmegaConf.masked_copy(cfg, ['ncpu','input','output',\
            'fitting','advanced'])
        with open('FAT_defaults.yml','w') as default_write:
            default_write.write(OmegaConf.to_yaml(no_cube))
        my_resources = import_pack_files('pyFAT_astro.config')
        data = (my_resources / 'FAT_Input_Catalogue.txt').read_bytes()
        with open('FAT_Example_Catalogue.txt','w+b') as default_write:
            default_write.write(data)
        print(f'''We have printed the file FAT_defaults.yml FAT_Input_Catalogue.txt in {os.getcwd()}.
''')
        sys.exit()
    #if a configuration file is provided read it
    if not cfg_input.configuration_file is None:
        succes = False
        while not succes:
            try:
                yaml_config = OmegaConf.load(cfg_input.configuration_file)
        #merge yml file with defaults
                cfg = OmegaConf.merge(cfg,yaml_config)
                succes = True
            except FileNotFoundError:
                cfg_input.configuration_file = input(f'''
You have provided a config file ({cfg_input.configuration_file}) but it can't be found.
If you want to provide a config file please give the correct name.
Else press CTRL-C to abort.
configuration_file = ''')



    cfg = OmegaConf.merge(cfg,inputconf)
    
    if not any([cfg.cube_name, cfg.configuration_file, cfg.installation_check\
                ,cfg.print_examples,cfg.input.catalogue]):
        from pyFAT_astro.Support.fat_errors import InputError
        raise InputError(f'''You have not provided a cube_name, configuration_file or input_catalogue
nor requested the installation_check or print_examples.
We are crashes because that means you want to do nothing.''')
    
    return cfg