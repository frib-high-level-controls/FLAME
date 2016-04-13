#!/usr/bin/env python

import sys
"[includes]"
from uscsi import Machine
"[includes]"

"[Parse lattice file]"
with open(sys.argv[1], "r") as F:
    mymachine = Machine(F)
"[Parse lattice file]"

"[Allocate State]"
thestate = mymachine.allocState({})
"[Allocate State]"

print "Initial state", thestate

"[Run simulation]"
mymachine.propagate(thestate)
"[Run simulation]"

print "Final state", thestate
