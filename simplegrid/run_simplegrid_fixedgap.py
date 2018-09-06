#!/usr/bin/env python2.7

"""
Script used to execute SimpleGrid package from command line.

"""


import sys

from simplegrid_fixedgap import main

if __name__ == '__main__':
    
    try:
        main()
    except KeyboardInterrupt as err:
        print(err)
        sys.exit()
        
        
