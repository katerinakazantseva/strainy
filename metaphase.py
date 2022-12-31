#!/usr/bin/env python3

import sys
import os
from metaphase.phase import phase_main
from metaphase.transform import transform_main

def main():
    #Setting executable paths
    metaphase_root = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, metaphase_root)

    #sys.exit(phase_main())
    sys.exit(transform_main())


if __name__ == "__main__":
    main()

