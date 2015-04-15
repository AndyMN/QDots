#!/usr/bin/python

import sys
import pstats

p = pstats.Stats(sys.argv[1])
p.sort_stats('time').print_stats(10)


