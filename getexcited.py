#/usr/bin/python

'''
 _________________________________________________________
|                                                         |
| The path below must be in string format, meaning        |
| PATHTOPATH = '<PATH>', where <PATH> is in single        |
| quotations.                                             |
|_________________________________________________________|

'''

PATHTOPACK = 'INSERT PATH TO getexcited_package HERE, DO NOT INCLUDE "getexcited_package" IN PATH'

'''
         _______________________________________
        |#######################################|
        |##|                                 |##|
        |##|        USE WITH CAUTION!        |##|
        |##|                                 |##|
        |##| Any questions or suggestions    |##|
        |##| regarding this script, feel     |##|
        |##| free to contact sifain@usc.edu. |##|
        |##|_________________________________|##|
        |#######################################|
 _________________________________________________________
|                                                         |
| Currently, there are 13 main functions:                 |
|                                                         |
| (1) prepare inputs for single-point calculations        |
| (2) generate a combined optical spectrum from           |
| single-point calculations                               |
| (3) prepare inputs for non-adiabatic excited-state      |
| molecular dynamics (NEXMD)                              |
| (4) prepare input files for an adiabatic simulation,    |
| with geometries taken from NEXMD                        |
| (5) collect populations as a function of time from hops |
| and  quantum coefficients                               |
| (6) collect PESs and NACTs from NEXMD                   |
| (7) prepare restart input files for NEXMD               |
| (8) clean the directories of unfinished trajectories    |
| (9) access options for geometry analysis                |
| (10) access options for dipole analysis                 |
| (11) access options for transition density analysis     |
| (12) access options for pump-push-probe spectroscopy    |
| (13) access code testing tools                          |
|                                                         |
| NOTE: Ground-state trajectory must be completed first.  |
| Coordinates and velocities are selected from a single   |
| ground-state trajectory. The files that contain         |
| coordinates and velocities are coords.xyz and           |
| velocity.out.                                           |
|                                                         |
| NOTE: Main output file from NEXMD program must be       |
| called 'md.out'.  Currently, this is how it is defined  |
| in the main submission script.  If this output file is  |
| renamed, 'md.out' in the 'collectceo.sh' script must    |
| be changed accordingly in order to generate an optical  |
| spectrum from single-point calculations.  Also,         |
| 'md.out' in 'optspec.py' must be changed.               |
|                                                         |
| NOTE: This script can generate a new set of random      |
| seeds or use a list provided by the user.  The latter   |
| option may be important for code testing or             |
| benchmarking purposes.  The user-defined list must      |
| strictly be a list of random seeds, with no header or   |
| footer, and the number of random seeds must be equal to |
| or greater than the number of trajectories requested.   |
|_________________________________________________________|

'''

import sys
import os
if not os.path.exists(PATHTOPACK):
    print 'You must provide the path to getexcited_package in getexcited.py (PATHTOPACK).'
    sys.exit()
sys.dont_write_bytecode = True
sys.path.append('%s' % (PATHTOPACK))
from getexcited_package.spcalc import SPCALC
from getexcited_package.optspec import OPTSPEC
from getexcited_package.nexmd import NEXMD
from getexcited_package.population import POPULATION
from getexcited_package.pesnact import PESNACT
from getexcited_package.restart import RESTART
from getexcited_package.newsim import NEWSIM
from getexcited_package.cleandir import CLEANDIR
from getexcited_package.dihedral import DIHEDRAL
from getexcited_package.bondlength import BONDLENGTH
from getexcited_package.bla import BLA
from getexcited_package.timing import TIMING
from getexcited_package.permdipole import PERMDIPOLE
from getexcited_package.tdiagonal import TDIAGONAL

FUNQ = input('\nSelect a task from the following list:\n\n[1] Prepare input files for single-point calculations\n[2] Generate an optical spectrum from single-point calculations\n[3] Prepare input files for NEXMD\n[4] Prepare input files for adiabatic dynamics with geometries from NEXMD\n[5] Collect populations from NEXMD\n[6] Collect PESs and NACTs from NEXMD\n[7] Prepare restart input files for NEXMD\n[8] Clean out the directories of NEXMD trajectories that are incomplete\n[9] Access options for geometry analysis\n[10] Access options for dipole analysis\n[11] Access options for transition density analysis\n[12] Access options for pump-push-probe spectroscopy (*** UNDER DEVELOPMENT, DO NOT USE ***)\n[13] Access code testing tools\n\nEnter the number corresponding to the desired task: ')
if FUNQ not in [1,2,3,4,5,6,7,8,9,10,11,12,13]:
    print 'Answer must be 1 through 13.'
    sys.exit()
if FUNQ == 1:
    SPCALC()
if FUNQ == 2:
    OPTSPEC(PATHTOPACK)
if FUNQ == 3:
    NEXMD()
if FUNQ == 4:
    NEWSIM()
if FUNQ == 5:
    POPULATION()
if FUNQ == 6:
    PESNACT()
if FUNQ == 7:
    RESTART(PATHTOPACK)
if FUNQ == 8:
    CLEANDIR()
if FUNQ == 9:
    ADVQ = input('\nSelect a task from the following list:\n\n[1] Calculate a dihedral angle\n[2] Calculate a bond length\n[3] Calculate a bond length alternation\n\nEnter the number corresponding to the desired task: ')
    if ADVQ not in [1,2,3]:
        print 'Answer must be 1 through 3.'
        sys.exit()
    if ADVQ == 1:
        DIHEDRAL()
    if ADVQ == 2:
        BONDLENGTH()
    if ADVQ == 3:
        BLA()
if FUNQ == 10:
    ADVQ = input('\nSelect a task from the following list:\n\n[1] Collect excited-state permanent dipole moment\n\nEnter the number corresponding to the desired task: ')
    if ADVQ != 1:
        print 'Answer must be 1.'
        sys.exit()
    if ADVQ == 1:
        PERMDIPOLE(PATHTOPACK)
if FUNQ == 11:
    ADVQ = input('\nSelect a task from the following list:\n\n[1] Analyze occupancy according to diagonal elements of the transition density matrix\n\nEnter the number corresponding to the desired task: ')
    if ADVQ not in [1]:
        print 'Answer must be 1.'
        sys.exit()
    if ADVQ == 1:
        TDIAGONAL()
if FUNQ == 12:
    sys.exit()
    ADVQ = input('\nSelect a task from the following list:\n\n[1] Prepare input files for single-point calculations after pump-push delay time\n[2] Generate optical spectrum from single-point calculations after pump-push delay time\n[3] Prepare input files for NEXMD after push pulse\n\nEnter the number corresponding to the desired task: ')
    if ADVQ not in [1,2,3]:
        print 'Answer must be 1 through 3.'
        sys.exit()
    if ADVQ == 1:
        SPCALC_PUSH()
    if ADVQ == 2:
        OPTSPEC_PUSH()
    if ADVQ == 3:
        NEXMD_PUSH()
if FUNQ == 13:
    ADVQ = input('Select a task from the following list:\n\n[1] Collect timing data from trajectories\n\nEnter the number corresponding to the desired task: ')
    if ADVQ != 1:
        print 'Answer must be 1.'
        sys.exit()
    if ADVQ == 1:
        TIMING(PATHTOPACK)
