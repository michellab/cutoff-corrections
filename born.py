#!/usr/bin/python
#
# Born ion solvation free energy. For debugging
#

from math import *


EPSFAC = 138.93569
eps = 78.40
q = 1.0
radius = 0.138
DELTAG= 0.5*EPSFAC*( -q**2 ) * (1 - 1 / eps) * ( 1 / radius )

print ("DELTAG = %8.5f kJ.mol-1 (%8.5f kcal.mol-1)" % (DELTAG, DELTAG/4.184))
