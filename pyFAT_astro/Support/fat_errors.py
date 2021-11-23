# -*- coding: future_fstrings -*-
#Define some errors
class BadCatalogueError(Exception):
    pass
class BadConfigurationError(Exception):
    pass
class BadSourceError(Exception):
    pass
class CfluxError(Exception):
    pass
class DefFileError(Exception):
    pass
class FileNotFoundError(Exception):
    pass
class FunctionCallError(Exception):
    pass
class InclinationRunError(Exception):
    pass
class InitializeError(Exception):
    pass
class InputError(Exception):
    pass
class ProgramError(Exception):
    pass
class NoConfigFile(Exception):
    pass
class SmallSourceError(Exception):
    pass
class SofiaFaintError(Exception):
    pass
class SofiaRunError(Exception):
    pass
class SupportRunError(Exception):
    pass
class TirificKillError(Exception):
    pass
