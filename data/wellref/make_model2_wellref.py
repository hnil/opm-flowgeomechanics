#!/usr/bin/env python3
"""Write a WELLREF variant of a model2 deck: one 5x5x1 ring around INJ1.

usage: make_model2_wellref.py IN.DATA OUTDIR
Copies IN.DATA's INCLUDE files next to the output.
"""
import os
import re
import shutil
import sys

src, outdir = sys.argv[1], sys.argv[2]
text = open(src).read()

# Refined wells get one connection per refined cell; fractures add more.
text, n = re.subn(r"(WELLDIMS\s*\n\s*\d+\s+)\d+", r"\g<1>800", text, count=1)
assert n == 1, "no WELLDIMS"

wellref = """WELLREF
-- well   NX NY NZ  RING1 RING2 RING3  LAYERS  KLAYERS
-- INJ1 is completed in (2,13), layers 13-16. One ring, one cell wide
-- (I 1-3, J 12-14: fault F2 is on its east face), 3 layers above and below.
 'INJ1'    5  5  1    1     0     0     PERF     3 /
/

"""
text, n = re.subn(r"^EDIT\s*$", wellref + "EDIT", text, count=1, flags=re.M)
assert n == 1, "no EDIT section"

os.makedirs(outdir, exist_ok=True)
name = os.path.splitext(os.path.basename(src))[0] + "_WELLREF.DATA"
open(os.path.join(outdir, name), "w").write(text)
for inc in re.findall(r"^INCLUDE\s*\n\s*'([^']+)'", text, flags=re.M):
    shutil.copy(os.path.join(os.path.dirname(src), inc), outdir)
print(os.path.join(outdir, name))
