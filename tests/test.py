#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os, sys
this_path, this_file = os.path.split(os.path.abspath(__file__))

def main(args):
    if args.foobar:
        print (args.foobar)

if __name__ == '__main__':
    import argparse, traceback
    try:
        parser = argparse.ArgumentParser(prog=this_file)
        parser.add_argument('--foobar')
        args = parser.parse_args()
        main(args)
        sys.exit(0)
    except KeyboardInterrupt as e:
        raise e
    except SystemExit as e:
        raise e
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)