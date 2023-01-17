#!/usr/bin/env python3
#This is to check the version of all imported packages

from load import *
import sys
import subprocess
import types
import re
import importlib

def get_pkgV():
    """
    # get the module and version list 
    pkgV = subprocess.check_output(["pip list"],shell=True)#.split("\n")
    pkgV = str(pkgV).split(r'\n')
    pkgV = [ re.sub(' +', ' ', string) for string in pkgV][2:]
    # the pkgDict [module name] = version number
    pkgDict = {}
    for string in pkgV:
        try:
            name, version = string.split(" ")
            pkgDict[name] = version
        except:
            pass
    """
    # get full name of global variables
    global_set = set()
    dict_global = globals()
    for k,v in dict_global.items():
        if isinstance(v, types.ModuleType):
            # get the module name in globals ; and in case of PIL.DrawImage
            module_ = v.__name__.split(".")[0]
            global_set.add(module_)
            

    module_set = set(sys.modules)
    module_set = set([i.split(".")[0] for i in module_set])
    modulenames =  module_set & global_set
    # get the version
    output_l = []#the output list for the install_requires
    for mn in modulenames:
        try:
            module2 = importlib.import_module(mn)
            version = module2.__version__
            output_l.append( module2.__name__ + ">=" + version  )
        except:
            print(mn,"is a built-in function")

    """
    allmodules = modulenames#[sys.modules[name].__name__ for name in modulenames]
    output_l = []#the output list for the install_requires
    for modulename in allmodules:
        try:
            version = pkgDict[modulename]
            output_l.append( modulename + ">=" + version  )
        except:
            pass
    
    """
    print(output_l)
    return output_l
