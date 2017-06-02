#!/usr/bin/env python

from __future__ import print_function
import sys

from heat_conduction.main import main as heat_main

if __name__ == '__main__':
    if sys.argv[1] == 'heat':
        heat_main()
