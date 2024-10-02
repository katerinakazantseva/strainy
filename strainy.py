#!/usr/bin/env python3

"""
This script sets up environment paths
and invokes Strainy without installation.
"""

import os
import sys


#Setting executable paths and additional import dirs
strainy_root = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, strainy_root)

flye_root = os.path.join(strainy_root, "submodules", "Flye")
sys.path.insert(0, flye_root)

bin_absolute = os.path.join(flye_root, "bin")   #for flye binaries
os.environ["PATH"] = bin_absolute + os.pathsep + os.environ["PATH"]
###

def main():
    #Strainy entry point
    import strainy.main
    sys.exit(strainy.main.main())




if __name__ == "__main__":
    main()
