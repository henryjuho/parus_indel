#!/usr/bin/env python

import anavar_utils as an
import sys

out = sys.argv[1]
target_files = [x.rstrip() for x in sys.stdin]
an.merge_results(target_files, out)
