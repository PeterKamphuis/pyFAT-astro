# -*- coding: future_fstrings -*-
#Define some errors
class BadCatalogueError(Exception):
    pass
class BadConfigurationError(Exception):
    pass
class BadCubeError(Exception):
    pass
class BadMaskError(Exception):
    pass
class BadSourceError(Exception):
    pass
class BadHeaderError(Exception):
    pass
class CfluxError(Exception):
    pass
class DefFileError(Exception):
    pass
class FaintSourceError(Exception):
    pass
class FileNotFoundError(Exception):
    pass
class FittingError(Exception):
    pass
class FunctionCallError(Exception):
    pass
class InclinationRunError(Exception):
    pass
class InitializeError(Exception):
    pass
class InputError(Exception):
    pass
class NoConfigFile(Exception):
    pass
class MissingProgramError(Exception):
    pass
class ProgramError(Exception):
    pass
class SmallSourceError(Exception):
    pass
class SofiaMissingError(Exception):
    pass
class SofiaRunError(Exception):
    pass
class SupportRunError(Exception):
    pass
class TirificOutputError(Exception):
    pass
class TirificKillError(Exception):
    pass
