#!/usr/bin/env python
# encoding: utf-8

import sys
from func_tdmt import tdmt
from myParser import parseMyLine
from myParser import checkConsistency,parseMyLine



###################################################
# ---  Load arguments
args=parseMyLine()
ll = sys.argv[1:]
if not ll:
       print "Use -h or --help option for Help"
       sys.exit(0)


###################################################
# ---  Check  arguments
checkConsistency(args)


###################################################
# --- Call timedoman MT inversion 
tdmt(args)
