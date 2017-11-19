#/usr/bin/python

'''

This function calculates excited-state permanent dipole
moment in time.

The dipole moment at every time-step is determined by the
dipole moment of the occupied state according to the
adiabatic or nonadiabatic dynamics and may be tracked
along a single trajectory or an ensemble of trajectories.

Type of calculation:

[1] Single Trajectory
< collection time (fs)
> time (fs), dipole moment (Debye)

[2] Ensemble of Trajectories

[2a] Mean
< collection time (fs)
> time (fs), dipole moment (Debye), standard deviation (Debye)

[2b] All time-steps
> trajectory directory, dipole moment (Debye) at all times and 
trajectories

Other options include:

[1] Relative direction of dipole
< line numbers of two atoms in the system, as presented
in 'input.ceon', which will be used to construct a vector and 
dotted with the dipole moment
> for option [1]: angle (degrees)
> for option [2a]: angle (degrees), standard deviation (degrees)
> for options [2b]: angle (degrees)

Output Files:
- permdipole_[type].out, where [type] = single, mean_ensemble, raw_ensemble

Error Files:
- permdipole_[type].err, where [type] = single, mean_ensemble, raw_ensemble

'''

import numpy as np
import os
import sys
import glob
import subprocess
import shlex
import fileinput

CWD = os.getcwd()

def PERMDIPOLE(PATHTODIPOLE):

    print 'Calculating excited-state permanent dipole moment as a function of time.'

    ## TYPE OF CALCULATION AND DIRECTORY CHECK ##
    DYNQ = input('Calculate dipole moment along one trajectory or an ensemble of trajectories?\nAnswer ONE [1] or ENSEMBLE [0]: ')
    if DYNQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if DYNQ == 0: ## ENSEMBLE
        NEXMDIR = raw_input('Ensemble directory [e.g. nexmd]: ')
        if not os.path.exists(NEXMDIR):
            print 'Path %s does not exist.' % (NEXMDIR)
            sys.exit()
        ## CHECK IF NEXMD FOLDERS EXIST ##
        NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
        NEXMDS.sort()
        if len(NEXMDS) == 0:
            print 'There are no NEXMD folders in %s.' % (NEXMDIR)
            sys.exit()
        ## DETERMINE MEAN OR ALL ##
        TYPEQ = input('Output mean dipole in time or output dipoles at all time-steps and trajectories?\nAnswer MEAN [0] or ALL [1]: ')
        if TYPEQ not in [0,1]:
            print 'Answer must be 0 or 1.'
            sys.exit()
    if DYNQ == 1: ## SINGLE TRAJECTORY
        TYPEQ = 0
        NEXMDIR = raw_input('Single trajectory directory: ')
        if not os.path.exists(NEXMDIR):
            print 'Path %s does not exist.' % (NEXMDIR)
            sys.exit()

    ## USER-DEFINED LENGTH OF ANALYSIS AND INITIALIZE ARRAYS ##
    if DYNQ == 0: ## ENSEMBLE
        if not os.path.exists('%s/header' % (NEXMDIR)):
            print 'Path %s/header does not exist.' % (NEXMDIR)
            sys.exit()
        HEADER = open('%s/header' % (NEXMDIR),'r')
        HEADER = HEADER.readlines()
    if DYNQ == 1: ## SINGLE TRAJECTORY
        if not os.path.exists('%s/input.ceon' % (NEXMDIR)):
            print 'Path %s/input.ceon does not exist.' % (NEXMDIR)
            sys.exit()
        HEADER = open('%s/input.ceon' % (NEXMDIR),'r')
        HEADER = HEADER.readlines()
    STATEINIT = None
    for LINE in HEADER:
        if 'bo_dynamics_flag' in LINE:
            BOFLAG = np.int(LINE.split()[0][len('bo_dynamics_flag='):-1])
        if 'exc_state_init=' in LINE:
            STATEINIT = np.int(LINE.split()[0][len('exc_state_init='):-1])
        if 'n_exc_states_propagate' in LINE:
            NSTATES = np.int(LINE.split()[0][len('n_exc_states_propagate='):-1])
        if 'time_init' in LINE:
            TINITH = np.float(LINE.split()[0][len('time_init='):-1])
        if 'time_step' in LINE:
            DT = np.float(LINE.split()[0][len('time_step='):-1])
        if 'n_class_steps' in LINE:
            TSMAX = np.int(LINE.split()[0][len('n_class_steps='):-1]) + 1
        if 'out_data_steps' in LINE:
            ODATA = np.int(LINE.split()[0][len('out_data_steps='):-1])
        if 'out_coords_steps' in LINE:
            CDATA = np.int(LINE.split()[0][len('out_coords_steps='):-1])
        if 'natoms' in LINE:
            NATOMS = np.int(LINE.split()[0][len('natoms='):-1])
    if BOFLAG == 1 and STATEINIT == None:
        print 'Dynamics are set to Born-Oppenheimer, but the initial state is not set.\nPlease check bo_dynamics_flag and exc_state_init in header.'
        sys.exit()
    
    ## COLLECTION TIME ##
    if TYPEQ == 0: ## MEAN DIPOLE
        if DYNQ == 0: ## ENSEMBLE
            TCOLL = input('Calculate dipole up to what time in femtoseconds?\nNote that averaged results will only include trajectories that are complete up to this time: ')
        if DYNQ == 1: ## SINGLE TRAJECTORY
            TCOLL = input('Calculate dipole up to what time in femtoseconds? ')
        if isinstance(TCOLL, int) == False and isinstance(TCOLL, float) == False:
            print 'Time must be integer or float.'
            sys.exit()
        if TCOLL < 0:
            print 'Time must be integer or float greater than zero.'
            sys.exit()
        TCOLL = np.float(TCOLL)
        if TCOLL > (TSMAX - 1)*DT:
            TCOLL = (TSMAX - 1)*DT
    if TYPEQ == 1: ## ALL DIPOLES
        TCOLL = (TSMAX - 1)*DT

    ## DETERMINE DIRECTION OF DIPOLE ##
    DOTQ = input('Find the angle between the dipole and a user-defined vector on the molecule?\nAnswer YES [1] or NO [0]: ')
    if DOTQ not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if DOTQ == 0:
        if NATOMS < 2:
            print 'Number of atoms set under natoms is less than zero.\nPlease check header (for ensemble) or input.ceon (for single trajectory).'
            sys.exit()
        else:
            LINES = [0,1]
    if DOTQ == 1:
        LINES = input('Input an array of the form [atom1, atom2], where atom# = line number of atom (0 is the first line).\nThese two atoms will be used to construct a vector: ')
        if isinstance(LINES, list) == False:
            print 'Input must be an array of the form [atom 1, atom2], where atom# = line number of atom (0 is the first line).'
            sys.exit()
        if len(LINES) != 2:
            print 'Input must be an array with two elements labeling the line numbers of two atoms.'
            sys.exit()
        INDEX = 0
        for i in LINES:
            if isinstance(i, int) == False:
                print 'Element number %d of input array must be integer.\nUser inputted [%s, %s], which is not allowed.' % (INDEX + 1, LINES[0], LINES[1])
                sys.exit()
            if i < 0:
                print 'Element number %d of input array must be a positive integer.\nUser inputted [%s, %s], which is not allowed.' % (INDEX + 1, LINES[0], LINES[1])
                sys.exit()
            if i > NATOMS - 1:
                print 'Element number %d of input array must be less than the max number of atoms (-1).\nUser inputted [%s, %s], which is not allowed.' % (INDEX + 1, LINES[0], LINES[1])
                sys.exit()
            INDEX += 1
        if len(np.unique(LINES)) != 2:
            print 'All elements of input array must be unique.\nUser inputted [%s, %s], which is not allowed.' % (LINES[0], LINES[1])
            sys.exit()

    ## NUMBER OF CLASSICAL TIME-STEPS ##
    TSCOL = 0
    while TSCOL*DT*ODATA <= TCOLL:
        TSCOL += 1
    CCOLL = 0
    
    ## NUMBER OF TIME-STEPS FOR COORDINATES ##
    NUM = 0
    while CCOLL <= TCOLL:
        CCOLL += DT*ODATA*CDATA
        NUM += 1
    
    ## NUMBER OF DIPOLES ##
    EDIPOLES = CCOLL

    ## COLLECTION TIME ARRAY ##
    TIMES = np.linspace(TINITH, CCOLL - DT*ODATA*CDATA, NUM)

    ## GREP PERMANENT DIPOLE MOMENTS AND STATES FROM ENSEMBLES ##
    if DYNQ == 0: ## ENSEMBLE
        print 'Checking permanent dipole moments and states.  Please wait ...'
        ## CHECKS TO MAKE SURE SCRIPTS ARE AVAILABLE ##
        if not os.path.exists('%s/getexcited_package/collectdipline.sh' % (PATHTODIPOLE)):
            print 'The script, collectdipline.sh, must be in the getexcited_package.'
            sys.exit()
        if not os.path.exists('%s/getexcited_package/collectpermdipole.sh' % (PATHTODIPOLE)):
            print 'The script, collectpermdipole.sh, must be in the getexcited_package.'
            sys.exit()
        ## GENERATION OF ERROR FILE ##
        ERROR = open('%s/dipole_collection_ensemble.err' % (CWD),'w')
        ERRFLAG = 0
        for NEXMD in NEXMDS:
            ## CHECK AND OPEN LIST OF DIRECTORIES ##
            if not os.path.exists('%s/%s/dirlist1' % (CWD,NEXMD)):
                print 'Path %s/%s/dirlist1 does not exist.' % (CWD,NEXMD)
                sys.exit()
            DIRLIST1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (CWD,NEXMD)))
            if isinstance(DIRLIST1,int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                ## CHECK IF DIRECTORY EXISTS ##
                if not os.path.exists('%s/%s/%04d' % (CWD,NEXMD,DIR)):
                    print >> ERROR, '%s%04d' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    continue
                ## GO TO DIRECTORY ##
                os.chdir('%s/%s/%04d' % (CWD,NEXMD,DIR))
                ## CHECK IF STANDARD OUTPUT EXISTS ##
                if not os.path.exists('%s/%s/%04d/md.out' % (CWD,NEXMD,DIR)):
                    print >> ERROR, '%s%04d/md.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    continue
                ## GREP LINE NUMBER OF CLASSICAL STEP AND CLASSICAL STEP ##
                subprocess.call(shlex.split('sh %s/getexcited_package/collectdipline.sh %d' % (PATHTODIPOLE, ODATA*CDATA)))
                if not os.path.exists('%s/%s/%04d/dipline.out' % (CWD,NEXMD,DIR)):
                    print >> ERROR, '%s%04d/dipline.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    continue
                print '%s%04d' % (NEXMD,DIR), 'dipole lines in md.out found'
                ## DATA = [LINE NUMBER OF CLASSICAL STEP, CLASSICAL STEP] ##
                DATA = np.genfromtxt('%s/%s/%04d/dipline.out' % (CWD,NEXMD,DIR))
                TDIPOLES = len(DATA)
                ## CHECK TO ENSURE DIPOLE CALCULATION ##
                if np.array_equal(np.around(DATA[1:EDIPOLES:1,1]*DT, decimals = 3), TIMES[1:EDIPOLES:1]) == False:
                    print >> ERROR, 'There is an inconsistency in time-step in %s%04d/dipline.out' % (NEXMD,DIR)
                    ERRFLAG = 1
                    continue
                ## DELETE PREVIOUS PERMDIPOLE.OUT IF EXISTS ##
                if os.path.exists('%s/%s/%04d/permdipole.out' % (CWD,NEXMD,DIR)):
                    os.remove('%s/%s/%04d/permdipole.out' % (CWD,NEXMD,DIR))
                ## GREP DIPOLES FROM STANDARD OUTPUT ##
                subprocess.call(shlex.split('sh %s/getexcited_package/collectpermdipole.sh %d' % (PATHTODIPOLE, NSTATES + 2)))
                if not os.path.exists('%s/%s/%04d/permdipole.out' % (CWD,NEXMD,DIR)):
                    print >> ERROR, '%s%04d/permdipole.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    continue
                print '%s%04d' % (NEXMD,DIR), 'dipoles in md.out extracted'
                ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
                with open('%s/%s/%04d/permdipole.out' % (CWD,NEXMD,DIR),'r') as DATA:
                    if len(DATA.readlines()) != TDIPOLES*(NSTATES + 3):
                        print >> ERROR, '%s%04d/permdipole.out' % (NEXMD,DIR), 'is incomplete'
                        ERRFLAG = 1
                        os.remove('%s/%s/%04d/permdipole.out' % (CWD,NEXMD,DIR))
                        continue
                ## DELETE PREVIOUS POP.OUT IF EXISTS ##
                if os.path.exists('%s/%s/%04d/pop.out' % (CWD,NEXMD,DIR)):
                    os.remove('%s/%s/%04d/pop.out' % (CWD,NEXMD,DIR))
                ## GET STATES ##
                if BOFLAG == 1: ## ADIABATIC
                    STATES = np.int_(np.ones(EDIPOLES)*STATINIT)
                if BOFLAG == 0: ## NONADIABATIC
                    ## CHECK COEFFICIENT FILE EXISTS ##
                    if not os.path.exists('%s/%s/%04d/coeff-n.out' % (CWD,NEXMD,DIR)):
                        print >> ERROR, '%s%04d/coeff-n.out' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        continue
                    DATA = open('%s/%s/%04d/coeff-n.out' % (CWD,NEXMD,DIR),'r')
                    DATA = DATA.readlines()
                    STATES = np.zeros(EDIPOLES)
                    INDEX = 0
                    for LINE in DATA[0:TSCOL:ODATA*CDATA]:
                        VAL = LINE.split()
                        PES = np.int(VAL[0])
                        TIME = np.around(np.float(VAL[1]), decimals = 3)
                        ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
                        if TIME != TIMES[INDEX]:
                            print >> ERROR, 'There is an inconsistency in time-step in %s%04d/coeff-n.out at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                            ERRFLAG = 1
                            break
                        STATES[INDEX] = PES
                        INDEX += 1
                    ## IF SIMULATION BECOMES ADIABATIC AFTER NONADIABATIC ##
                    if len(DATA[0:TSCOL:ODATA*CDATA]) < EDIPOLES:
                        STATES[INDEX::] = PES
                ## SAVE POPULATIONS TO FILE ##
                np.savetxt('pop.out', np.transpose([TIMES,STATES]), fmt=['%10.5e','%d'])
        if ERRFLAG == 1:
            print 'One or more trajectories have experienced an error, check dipole_collection_ensemble.err.'
            CONTQ = input('Continue? Answer YES [1] or NO [0]: ')
            if CONTQ not in [1,0]:
                print 'Answer must to be 1 or 0.'
                sys.exit()
            if CONTQ == 0:
                sys.exit()
        else:
            os.remove('%s/dipole_collection_ensemble.err' % (CWD))

    ## GREP PERMANENT DIPOLE MOMENTS AND STATES FROM SINGLE TRAJECTORY ##
    if DYNQ == 1: ## SINGLE TRAJECTORY
        print 'Checking permanent dipole moments and states.  Please wait ...'
        ## CHECKS TO MAKE SURE SCRIPTS ARE AVAILABLE ##
        if not os.path.exists('%s/getexcited_package/collectdipline.sh' % (PATHTODIPOLE)):
            print 'The script, collectdipline.sh, must be in the getexcited_package.'
            sys.exit()
        if not os.path.exists('%s/getexcited_package/collectpermdipole.sh' % (PATHTODIPOLE)):
            print 'The script, collectpermdipole.sh, must be in the getexcited_package.'
            sys.exit()
        ## GO TO DIRECTORY ##
        os.chdir('%s/%s' % (CWD,NEXMDIR))
        ## CHECK IF STANDARD OUTPUT EXISTS ##
        if not os.path.exists('%s/%s/md.out' % (CWD,NEXMDIR)):
            print '%s/md.out' % (NEXMDIR), 'does not exist'
            sys.exit()
        ## GREP LINE NUMBER OF CLASSICAL STEP AND CLASSICAL STEP ##
        subprocess.call(shlex.split('sh %s/getexcited_package/collectdipline.sh %d' % (PATHTODIPOLE, ODATA*CDATA)))
        if not os.path.exists('%s/%s/dipline.out' % (CWD,NEXMDIR)):
            print '%s/dipline.out' % (NEXMDIR), 'does not exist'
            sys.exit()
        print '%s' % (NEXMDIR), 'dipole lines in md.out found'
        ## DATA = [LINE NUMBER OF CLASSICAL STEP, CLASSICAL STEP] ##
        DATA = np.genfromtxt('%s/%s/dipline.out' % (CWD,NEXMDIR))
        TDIPOLES = len(DATA)
        ## CHECK TO ENSURE DIPOLE CALCULATION ##
        if np.array_equal(np.around(DATA[1:EDIPOLES:1,1]*DT, decimals = 3), TIMES[1:EDIPOLES:1]) == False:
            print 'There is an inconsistency in time-step in %s/dipline.out' % (NEXMDIR)
            sys.exit()
        ## DELETE PREVIOUS PERMDIPOLE.OUT IF EXISTS ##
        if os.path.exists('%s/%s/permdipole.out' % (CWD,NEXMDIR)):
            os.remove('%s/%s/permdipole.out' % (CWD,NEXMDIR))
        ## GREP DIPOLES FROM STANDARD OUTPUT ##
        subprocess.call(shlex.split('sh %s/getexcited_package/collectpermdipole.sh %d' % (PATHTODIPOLE, NSTATES + 2)))
        if not os.path.exists('%s/%s/permdipole.out' % (CWD,NEXMDIR)):
            print '%s/permdipole.out' % (NEXMDIR), 'does not exist'
            sys.exit()
        print '%s' % (NEXMDIR), 'dipoles in md.out extracted'
        ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
        with open('%s/%s/permdipole.out' % (CWD,NEXMDIR),'r') as DATA:
            if len(DATA.readlines()) != TDIPOLES*(NSTATES + 3):
                print '%s/permdipole.out' % (NEXMDIR), 'is incomplete'
                sys.exit()
        ## DELETE PREVIOUS POP.OUT IF EXISTS ##
        if os.path.exists('%s/%s/pop.out' % (CWD,NEXMDIR)):
            os.remove('%s/%s/pop.out' % (CWD,NEXMDIR))
        ## GET STATES ##
        if BOFLAG == 1: ## ADIABATIC
            STATES = np.int_(np.ones(EDIPOLES)*STATINIT)
        if BOFLAG == 0: ## NONADIABATIC
            ## CHECK COEFFICIENT FILE EXISTS ##
            if not os.path.exists('%s/%s/coeff-n.out' % (CWD,NEXMDIR)):
                print '%s/coeff-n.out' % (NEXMDIR), 'does not exist'
                sys.exit()
            DATA = open('%s/%s/coeff-n.out' % (CWD,NEXMDIR),'r')
            DATA = DATA.readlines()
            STATES = np.zeros(EDIPOLES)
            INDEX = 0
            for LINE in DATA[0:TSCOL:ODATA*CDATA]:
                VAL = LINE.split()
                PES = np.int(VAL[0])
                TIME = np.around(np.float(VAL[1]), decimals = 3)
                ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
                if TIME != TIMES[INDEX]:
                    print 'There is an inconsistency in time-step in %s/coeff-n.out at %.3f fs' % (NEXMDIR,TIMES[INDEX])
                    sys.exit()
                STATES[INDEX] = PES
                INDEX += 1
            ## IF SIMULATION BECOMES ADIABATIC AFTER NONADIABATIC ##
            if len(DATA[0:TSCOL:ODATA*CDATA]) < EDIPOLES:
                STATES[INDEX::] = PES
        ## SAVE POPULATIONS TO FILE ##
        np.savetxt('pop.out', np.transpose([TIMES,STATES]), fmt=['%10.5e','%d'])
            
    ## COLLECT USER-DEFINED VECTOR IN THE MOLECULE TO DETERMINE DIPOLE DIRECTION FOR ENSEMBLES ##
    if DYNQ == 0: ## ENSEMBLE
        os.chdir('%s' % (CWD))
        print 'Checking coordinate files for vector generation.  Please wait ...'
        ## GENERATION OF ERROR FILE ##
        ERROR = open('%s/uservec_collection.err' % (CWD),'w')
        ERRFLAG = 0
        for NEXMD in NEXMDS:
            ## CHECK AND OPEN LIST OF DIRECTORIES ##
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %s/dirlist1 does not exist.' % (NEXMD)
                sys.exit()
            DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(DIRLIST1,int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                ## CHECK IF TRAJECTORY DIRECTORY EXISTS ##
                if not os.path.exists('%s/%04d' % (NEXMD,DIR)):
                    print >> ERROR, '%s%04d' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    continue
                ## CHECK IF COORDINATE FILE EXISTS ##
                if not os.path.exists('%s/%04d/coords.xyz' % (NEXMD,DIR)):
                    print >> ERROR, '%s%04d/coords.xyz' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    continue
                ## FIND GEOMETRIES ##
                DATA = open('%s/%04d/coords.xyz' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                LENC = len(DATA)
                NCOORDS = 0
                INDEX = 0
                TFLAG1 = 0
                TFLAG2 = 0
                TFLAG3 = 0
                ARRAY = np.array([])
                for LINE in DATA:
                    if 'time' in LINE:
                        if NCOORDS == 0:
                            TINIT = np.float(LINE.split()[-1])
                            if TINIT != TINITH:
                                TFLAG1 = 1
                                break
                        else:
                            TIME = np.around(np.float(LINE.split()[-1]), decimals = 3)
                            if TIME > TCOLL:
                                TFLAG3 = 1
                                break
                            if TIME != TIMES[NCOORDS]:
                                TFLAG2 = 1
                                break
                        NCOORDS += 1
                        ARRAY = np.append(ARRAY,INDEX)
                    INDEX += 1
                if TFLAG1 == 1:
                    print >> ERROR, 'Initial time in %s%04d/coords.xyz does not match time_init in %s%04d/input.ceon.' % (NEXMD,DIR)
                    ERRFLAG = 1
                    continue
                if TFLAG2 == 1:
                    print >> ERROR, 'There is an inconsistency in time-step in %s%04d/coords.xyz.' % (NEXMD,DIR)
                    ERRFLAG = 1
                    continue
                if TFLAG3 == 1:
                    ARRAY = np.append(ARRAY,INDEX)
                else:
                    ARRAY = np.append(ARRAY, LENC + 1)
                ARRAY = np.int_(ARRAY)
                ## CHECKS TO ENSURE USER-DEFINED VECTOR CALCULATION ##
                if NCOORDS == 0:
                    print >> ERROR, 'No coordinates were found in %s%04d.' % (NEXMD,DIR)
                    ERRFLAG = 1
                    continue
                if NCOORDS == 1:
                    print >> ERROR, 'Only initial coordinates, at %.2f fs, were found in %s%04d.' % (TINIT,NEXMD,DIR)
                    ERRFLAG = 1
                    continue
                ## GENERATION OF USER-DEFINED VECTOR FILE ##
                OUTPUT = open('%s/%04d/uservec.out' % (NEXMD,DIR),'w')
                ## PRINT USER-DEFINED VECTOR TO FILE ##
                for NCOORD in np.arange(NCOORDS):
                    COORDS = DATA[ARRAY[NCOORD] + 1: ARRAY[NCOORD + 1] - 1:1]
                    VEC0 = np.float_(COORDS[LINES[0]].split()[1:])
                    VEC1 = np.float_(COORDS[LINES[1]].split()[1:])
                    UVEC = np.subtract(VEC1, VEC0)
                    OUTPUT.write('{:>12}  {:>12}  {:>12}  {:>12}\n'.format(TIMES[NCOORD],UVEC[0],UVEC[1],UVEC[2]))
                print '%s%04d' % (NEXMD,DIR), 'user-defined vector from coords.xyz extracted'
        if ERRFLAG == 1:
            print 'One or more trajectories have experienced an error, check uservec_collection.err.'
            CONTQ = input('Continue? Answer YES [1] or NO [0]: ')
            if CONTQ not in [1,0]:
                print 'Answer must to be 1 or 0.'
                sys.exit()
            if CONTQ == 0:
                sys.exit()
        else:
            os.remove('%s/uservec_collection.err' % (CWD))

    ## COLLECT USER-DEFINED VECTOR IN THE MOLECULE TO DETERMINE DIPOLE DIRECTION FOR SINGLE TRAJECTORY ##
    if DYNQ == 1: ## SINGLE TRAJECTORY
        os.chdir('%s' % (CWD))
        print 'Checking coordinate files for vector generation.  Please wait ...'
        ## CHECK IF TRAJECTORY DIRECTORY EXISTS ##
        if not os.path.exists('%s' % (NEXMDIR)):
            print '%s' % (NEXMDIR), 'does not exist'
            sys.exit()
        ## CHECK IF COORDINATE FILE EXISTS ##
        if not os.path.exists('%s/coords.xyz' % (NEXMDIR)):
            print '%s/coords.xyz' % (NEXMDIR), 'does not exist'
            sys.exit()
        ## FIND GEOMETRIES ##
        DATA = open('%s/coords.xyz' % (NEXMDIR),'r')
        DATA = DATA.readlines()
        LENC = len(DATA)
        NCOORDS = 0
        INDEX = 0
        TFLAG1 = 0
        TFLAG2 = 0
        TFLAG3 = 0
        ARRAY = np.array([])
        for LINE in DATA:
            if 'time' in LINE:
                if NCOORDS == 0:
                    TINIT = np.float(LINE.split()[-1])
                    if TINIT != TINITH:
                        TFLAG1 = 1
                        break
                else:
                    TIME = np.around(np.float(LINE.split()[-1]), decimals = 3)
                    if TIME > TCOLL:
                        TFLAG3 = 1
                        break
                    if TIME != TIMES[NCOORDS]:
                        TFLAG2 = 1
                        break
                NCOORDS += 1
                ARRAY = np.append(ARRAY,INDEX)
            INDEX += 1
        if TFLAG1 == 1:
            print >> ERROR, 'Initial time in %s/coords.xyz does not match time_init in %s/input.ceon.' % (NEXMDIR)
            sys.exit()
        if TFLAG2 == 1:
            print >> ERROR, 'There is an inconsistency in time-step in %s/coords.xyz.' % (NEXMDIR)
            sys.exit()
        if TFLAG3 == 1:
            ARRAY = np.append(ARRAY,INDEX)
        else:
            ARRAY = np.append(ARRAY, LENC + 1)
        ARRAY = np.int_(ARRAY)
        ## CHECKS TO ENSURE USER-DEFINED VECTOR CALCULATION ##
        if NCOORDS == 0:
            print >> ERROR, 'No coordinates were found in %s/coords.xyz.' % (NEXMDIR)
            sys.exit()
        if NCOORDS == 1:
            print >> ERROR, 'Only initial coordinates, at %.2f fs, were found in %s/coords.xyz.' % (TINIT,NEXMDIR)
            sys.exit()
        ## GENERATION OF USER-DEFINED VECTOR FILE ##
        OUTPUT = open('%s/uservec.out' % (NEXMDIR),'w')
        ## PRINT USER-DEFINED VECTOR TO FILE ##
        for NCOORD in np.arange(NCOORDS):
            COORDS = DATA[ARRAY[NCOORD] + 1: ARRAY[NCOORD + 1] - 1:1]
            VEC0 = np.float_(COORDS[LINES[0]].split()[1:])
            VEC1 = np.float_(COORDS[LINES[1]].split()[1:])
            UVEC = np.subtract(VEC1, VEC0)
            OUTPUT.write('{:>12}  {:>12}  {:>12}  {:>12}\n'.format(TIMES[NCOORD],UVEC[0],UVEC[1],UVEC[2]))
        print '%s' % (NEXMDIR), 'user-defined vector from coords.xyz extracted'
            
    ## COLLECT DIPOLE MOMENT ALONG A SINGLE TRAJECTORY ##
    if DYNQ == 1: ## SINGLE TRAJECTORY
        os.chdir('%s' % (CWD))
        print 'Collecting permanent dipole moment along a single trajectory.  Please wait ...'
        ## GENERATE OUTPUT FILE ##
        OUTPUT = open('%s/permdipole_single.out' % (CWD),'w')
        TTRAJ = 0
        CTRAJ = 0
        ETRAJ = 0
        ## DETERMINE COMPLETED NUMBER OF TIME-STEPS ##
        if not os.path.exists('%s/energy-ev.out' % (NEXMDIR)):
            print 'Path %s/energy-ev.out does not exist.' % (NEXMDIR)
            sys.exit()
        DATA = open('%s/energy-ev.out' % (NEXMDIR),'r')
        DATA = DATA.readlines()
        TSTEPS = len(DATA) - 1
        ## GENERATE ARRAY WITH INDICES OF THE DIPOLE BLOCKS ALONG TRAJECTORY ##
        if TSTEPS >= TSCOL:
            if not os.path.exists('%s/permdipole.out' % (NEXMDIR)):
                print 'Path %s/permdipole.out does not exist.' % (NEXMDIR)
                sys.exit()
            DATA = open('%s/permdipole.out' % (NEXMDIR),'r')
            DATA = DATA.readlines()
            LEND = len(DATA)
            NDIPOLES = 0
            TFLAG = 0
            INDEX = 0
            ARRAY = np.array([])
            for LINE in DATA:
                if 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)' in LINE:
                    if NDIPOLES == 0:
                        TIME = TINITH
                    else:
                        TIME += DT*ODATA*CDATA
                        if TIME > TCOLL:
                            TFLAG = 1
                            break
                    NDIPOLES += 1
                    ARRAY = np.append(ARRAY,INDEX)
                INDEX += 1
            ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
            if NDIPOLES != EDIPOLES:
                print 'Number of dipoles detected in %s/permdipole.out, %d, does not match the expected %d.' % (NEXMDIR,NDIPOLES,EDIPOLES)
                sys.exit()
            ## APPEND LINES FOR LAST DIPOLE SET ##
            if TFLAG == 1:
                ARRAY = np.append(ARRAY, INDEX)
            else:
                ARRAY = np.append(ARRAY, LEND + 1)
            ARRAY = np.int_(ARRAY)
            ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
            if NDIPOLES == 0:
                print 'No dipoles were found in %s/permdipole.out' % (NEXMDIR)
                sys.exit()
            ## OPEN THE USER-DEFINED VECTOR FILE = [TIME, Vx, Vy, Vz] ##
            if not os.path.exists('%s/uservec.out' % (NEXMDIR)):
                print 'Path %s/uservec.out' % (NEXMDIR), 'does not exist'
                sys.exit()
            USERVEC = np.genfromtxt('%s/uservec.out' % (NEXMDIR))
            ## OPEN POPULATION DATA = [TIME, STATE] ##
            STATES = np.genfromtxt('%s/pop.out' % (NEXMDIR))
            ## COLLECT DIPOLE ALONG A SINGLE TRAJECTORY ##
            SDIPOLE = np.zeros((NDIPOLES,2))
            for NDIPOLE in np.arange(NDIPOLES):
                DIPOLES = DATA[ARRAY[NDIPOLE] + 1:ARRAY[NDIPOLE + 1]:1]
                VDIPOLE = np.float_(DIPOLES[np.int(STATES[NDIPOLE, 1])].split()[2:5:1]) ## EXTRACTS DIPOLE VECTOR
                MDIPOLE = np.float(DIPOLES[np.int(STATES[NDIPOLE, 1])].split()[5]) ## EXTRACTS DIPOLE MAGNITUDE
                COSINE = np.arccos(np.dot(VDIPOLE, USERVEC[NDIPOLE,1:4:1])/(np.linalg.norm(VDIPOLE)*np.linalg.norm(USERVEC[NDIPOLE,1:4:1]))) if DOTQ == 1 else 0 ## INVERSE COSINE OF THE DOT PRODUCT
                SDIPOLE[NDIPOLE] = np.array([MDIPOLE, COSINE])
            ## DELETE EXTRANEOUS DATA ##
            os.remove('%s/permdipole.out' % (NEXMDIR))
            os.remove('%s/pop.out' % (NEXMDIR))
            os.remove('%s/uservec.out' % (NEXMDIR))
            print '%s' % (NEXMDIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
            CTRAJ += 1
            if TSTEPS == TSMAX:
                ETRAJ += 1
        else:
            print '%s' % (NEXMDIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
        TTRAJ += 1
        ## SUMMARY OF RESULTS ##
        if CTRAJ == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(TSMAX)), TCOLL)
        else:
            print 'Total Trajectories:', '%04d' % (TTRAJ)
            print 'Completed Trajectories:', '%04d' % (CTRAJ)
            print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            print >> OUTPUT, 'Total Trajectories: ', '%04d' % (TTRAJ)
            print >> OUTPUT, 'Completed Trajectories: ', '%04d' % (CTRAJ)
            print >> OUTPUT, 'Excellent Trajectories: ', '%04d' % (ETRAJ)
            for NDIPOLE in np.arange(NDIPOLES):
                print >> OUTPUT, '%0*.2f' % (len(str((TSMAX))) + 2, DT*ODATA*CDATA*NDIPOLE), '%d' % (STATES[NDIPOLE, 1]), '%03.6f' % (SDIPOLE[NDIPOLE,0]), '%.6f' % (np.degrees(SDIPOLE[NDIPOLE,1]))

    ## CALCULATE MEAN DIPOLE FROM ENSEMBLE OF TRAJECTORIES ##
    if DYNQ == 0 and TYPEQ == 0: ## MEAN FROM ENSEMBLE
        os.chdir('%s' % (CWD))
        print 'Collecting mean permanent dipole moment from ensemble.  Please wait ...'
        ## DETERMINE TOTAL NUMBER OF TRAJECTORIES IN ENSEMBLE ##
        with open('%s/totdirlist' % (NEXMDIR),'w') as DATA:
            for NEXMD in NEXMDS:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %sdirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                INPUT = fileinput.input('%s/dirlist1' % (NEXMD))
                DATA.writelines(INPUT)
        DIRLIST1 = np.int_(np.genfromtxt('%s/totdirlist' % (NEXMDIR)))
        if isinstance(DIRLIST1,int) == True:
            DIRLIST1 = np.array([DIRLIST1])
        os.remove('%s/totdirlist' % (NEXMDIR))
        ## GENERATE OUTPUT AND ERROR FILES ##
        OUTPUT = open('%s/permdipole_mean_ensemble.out' % (CWD),'w')
        ERROR = open('%s/permdipole_mean_ensemble.err' % (CWD),'w')
        ## GENERATE DIPOLE ARRAYS FOR FINAL RESULTS ##
        FDIPOLE = np.zeros(len(TIMES)) ## MAGNITUDE
        FCOSINE = np.zeros(len(TIMES)) ## DIRECTION
        EDIPOLE = np.zeros((EDIPOLES, len(DIRLIST1)))
        ECOSINE = np.zeros((EDIPOLES, len(DIRLIST1)))
        TTRAJ = 0
        CTRAJ = 0
        ETRAJ = 0
        ERRFLAG = 0
        for NEXMD in NEXMDS:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %sdirlist1 does not exist.' % (NEXMD)
                sys.exit()
            DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(DIRLIST1, int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                ## DETERMINE COMPLETED NUMBER OF TIME-STEPS ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,DIR)):
                    print >> ERROR, '%s%04d/energy-ev.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                DATA = open('%s/%04d/energy-ev.out' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                TSTEPS = len(DATA) - 1
                ## GENERATE ARRAY WITH INDICES OF THE DIPOLE BLOCKS ALONG A SINGLE TRAJECTORY ##
                if TSTEPS >= TSCOL:
                    if not os.path.exists('%s/%04d/permdipole.out' % (NEXMD,DIR)):
                        print >> ERROR, '%s%04d/permdipole.out' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    DATA = open('%s/%04d/permdipole.out' % (NEXMD,DIR),'r')
                    DATA = DATA.readlines()
                    LEND = len(DATA)
                    NDIPOLES = 0
                    TFLAG = 0
                    INDEX = 0
                    ARRAY = np.array([])
                    for LINE in DATA:
                        if 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)' in LINE:
                            if NDIPOLES == 0:
                                TIME = TINITH
                            else:
                                TIME += DT*ODATA*CDATA
                                if TIME > TCOLL:
                                    TFLAG = 1
                                    break
                            NDIPOLES += 1
                            ARRAY = np.append(ARRAY,INDEX)
                        INDEX += 1
                    ## CHECK TO ENSURE DIPOLE CALCULATION ##
                    if NDIPOLES != EDIPOLES:
                        print >> ERROR, 'Number of dipoles detected in %s%04d/permdipole.out, %d, does not match the expected %d.' % (NEXMD,DIR,NDIPOLES,EDIPOLES)
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    ## APPEND LINES FOR LAST DIPOLE SET ##
                    if TFLAG == 1:
                        ARRAY = np.append(ARRAY,INDEX)
                    else:
                        ARRAY = np.append(ARRAY, LEND + 1)
                    ARRAY = np.int_(ARRAY)
                    ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
                    if NDIPOLES == 0:
                        print >> ERROR, 'No dipoles were found in %s%04d/permdipole.out' % (NEXMD,DIR)
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    ## OPEN THE USER-DEFINED VECTOR FILE = [TIME, Vx, Vy, Vz] ##
                    if not os.path.exists('%s/%04d/uservec.out' % (NEXMD,DIR)):
                        print >> ERROR, 'Path %s/%04d/uservec.out' % (NEXMD,DIR)
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    USERVEC = np.genfromtxt('%s/%04d/uservec.out' % (NEXMD,DIR))
                    ## OPEN POPULATION DATA = [TIME, STATE] ##
                    STATES = np.genfromtxt('%s/%04d/pop.out' % (NEXMD,DIR))
                    ## COLLECT DIPOLE ALONG A SINGLE TRAJECTORY ##
                    SDIPOLE = np.zeros(NDIPOLES)
                    SCOSINE = np.zeros(NDIPOLES)
                    for NDIPOLE in np.arange(NDIPOLES):
                        DIPOLES = DATA[ARRAY[NDIPOLE] + 1:ARRAY[NDIPOLE + 1]:1]
                        VDIPOLE = np.float_(DIPOLES[np.int(STATES[NDIPOLE, 1])].split()[2:5:1]) ## EXTRACTS DIPOLE VECTOR
                        MDIPOLE = np.float(DIPOLES[np.int(STATES[NDIPOLE, 1])].split()[5]) ## EXTRACTS DIPOLE MAGNITUDE
                        COSINE = np.arccos(np.dot(VDIPOLE, USERVEC[NDIPOLE,1:4:1])/(np.linalg.norm(VDIPOLE)*np.linalg.norm(USERVEC[NDIPOLE,1:4:1]))) if DOTQ == 1 else 0 ## INVERSE COSINE OF THE DOT PRODUCT
                        SDIPOLE[NDIPOLE] = MDIPOLE
                        SCOSINE[NDIPOLE] = COSINE
                        EDIPOLE[NDIPOLE,CTRAJ] = SDIPOLE[NDIPOLE]
                        ECOSINE[NDIPOLE,CTRAJ] = SCOSINE[NDIPOLE]
                    FDIPOLE += SDIPOLE
                    FCOSINE += SCOSINE
                    ## DELETE EXTRANEOUS DATA ##
                    os.remove('%s/%04d/permdipole.out' % (NEXMD,DIR))
                    os.remove('%s/%04d/pop.out' % (NEXMD,DIR))
                    os.remove('%s/%04d/uservec.out' % (NEXMD,DIR))
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    CTRAJ += 1
                    if TSTEPS == TSMAX:
                        ETRAJ += 1
                else:
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                TTRAJ += 1
        ## SUMMARY OF RESULTS ##
        if CTRAJ == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(TSMAX)), TCOLL)
        else:
            ## MEAN AND STANDARD DEVIATION FOR DIPOLE MAGNITUDE ##
            EDIPOLE = np.delete(EDIPOLE, np.arange(CTRAJ, TTRAJ), axis = 1)
            EDIPOLE = np.std(EDIPOLE, axis = 1)
            FDIPOLE = FDIPOLE/CTRAJ
            ## MEAN AND STANDARD DEVIATION FOR DIPOLE DIRECTION ##
            ECOSINE = np.delete(ECOSINE, np.arange(CTRAJ, TTRAJ), axis = 1)
            ECOSINE = np.std(ECOSINE, axis = 1)
            FCOSINE = FCOSINE/CTRAJ
            ## PRINT FINAL RESULTS ##
            print 'Total Trajectories:', '%04d' % (TTRAJ)
            print 'Completed Trajectories:', '%04d' % (CTRAJ)
            print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            print >> OUTPUT, 'Total Trajectories:', '%04d' % (TTRAJ)
            print >> OUTPUT, 'Completed Trajectories:', '%04d' % (CTRAJ)
            print >> OUTPUT, 'Excellent Trajectories:', '%04d' % (ETRAJ)
            for NDIPOLE in np.arange(NDIPOLES):
                print >>  OUTPUT, '%0*.2f' % (len(str((TSMAX))) + 2, DT*ODATA*CDATA*(NDIPOLE)), '%03.6f' % (FDIPOLE[NDIPOLE]), '%03.6f' % (EDIPOLE[NDIPOLE]), '%.6f' % (np.degrees(FCOSINE[NDIPOLE])), '%.6f' % (np.degrees(ECOSINE[NDIPOLE]))
        if ERRFLAG == 1:
            print 'One or more trajectories have experienced an error, check permdipole_mean_ensemble.err.'
        else:
            os.remove('%s/permdipole_mean_ensemble.err' % (CWD))
                
    ## COLLECT DIPOLES FROM ENSEMBLE OF TRAJECTORIES AT ALL TIME-STEPS ##
    if DYNQ == 0 and TYPEQ == 1: ## ALL FROM ENSEMBLE
        os.chdir('%s' % (CWD))
        print 'Collecting all permanent dipole moments from ensemble.  Please wait ...'
        ## GENERATE OUTPUT AND ERROR FILES ##
        OUTPUT = open('%s/permdipole_raw_ensemble.out' % (CWD),'w')
        ERROR = open('%s/permdipole_raw_ensemble.err' % (CWD),'w')
        TTRAJ = 0
        ETRAJ = 0
        ERRFLAG = 0
        for NEXMD in NEXMDS:
            if not os.path.exists('%s/dirlist1' % (NEXMD)):
                print 'Path %sdirlist1 does not exist.' % (NEXMD)
                sys.exit()
            DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
            if isinstance(DIRLIST1, int) == True:
                DIRLIST1 = np.array([DIRLIST1])
            for DIR in DIRLIST1:
                ## DETERMINE NUMBER OF TIME-STEPS COMPLETED ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,DIR)):
                    print >> ERROR, '%s%04d/energy-ev.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                DATA = open('%s/%04d/energy-ev.out' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                TSTEPS = len(DATA) - 1
                ## GENERATE ARRAY WITH INDICES OF THE DIPOLE BLOCKS ALONG TRAJECTORY ##
                if not os.path.exists('%s/%04d/permdipole.out' % (NEXMD,DIR)):
                    print >> ERROR, '%s%04d/permdipole.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                DATA = open('%s/%04d/permdipole.out' % (NEXMD,DIR),'r')
                DATA = DATA.readlines()
                LEND = len(DATA)
                NDIPOLES = 0
                TFLAG = 0
                INDEX = 0
                ARRAY = np.array([])
                for LINE in DATA:
                    if 'Frequencies (eV) and Total Molecular Dipole Moments (Debye)' in LINE:
                        if NDIPOLES == 0:
                            TIME = TINITH
                        else:
                            TIME += DT*ODATA*CDATA
                        NDIPOLES += 1
                        ARRAY = np.append(ARRAY,INDEX)
                    INDEX += 1
                ## ANOTHER CHECK TO ENSURE DIPOLE CALCULATION ##
                if NDIPOLES != EDIPOLES:
                    print >> ERROR, 'Number of dipoles detected in %s%04d/permdipole.out, %d, does not match the expected %d.' % (NEXMD,DIR,NDIPOLES,EDIPOLES)
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                ## APPEND LINES FOR LAST DIPOLE SET ##
                ARRAY = np.append(ARRAY, LEND + 1)
                ARRAY = np.int_(ARRAY)
                ## ANOTHER CHECK TO ENSURE THE DIPOLE CALCULATION ##
                if NDIPOLES == 0:
                    print >> ERROR, 'No dipoles were found in %s%04d/permdipole.out' % (NEXMD,DIR)
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                ## OPEN THE USER-DEFINED VECTOR FILE = [TIME, Vx, Vy, Vz] ##
                if not os.path.exists('%s/%04d/uservec.out' % (NEXMD,DIR)):
                    print >> ERROR, 'Path %s/%04d/uservec.out' % (NEXMD,DIR), 'does not exist'
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                USERVEC = np.genfromtxt('%s/%04d/uservec.out' % (NEXMD,DIR))
                ## OPEN POPULATION DATA = [TIME, STATE] ##
                if not os.path.exists('%s/%04d/pop.out' % (NEXMD,DIR)):
                    print >> ERROR, 'Path %s/%04d/pop.out' % (NEXMD,DIR)
                    ERRFLAG = 1
                    TTRAJ += 1
                    continue
                STATES = np.genfromtxt('%s/%04d/pop.out' % (NEXMD,DIR))
                ## COLLECT DIPOLE ALONG A SINGLE TRAJECTORY ##
                for NDIPOLE in np.arange(NDIPOLES):
                    DIPOLES = DATA[ARRAY[NDIPOLE] + 1:ARRAY[NDIPOLE + 1]:1]
                    VDIPOLE = np.float_(DIPOLES[np.int(STATES[NDIPOLE, 1])].split()[2:5:1]) ## EXTRACTS DIPOLE VECTOR
                    MDIPOLE = np.float(DIPOLES[np.int(STATES[NDIPOLE, 1])].split()[5]) ## EXTRACTS DIPOLE MAGNITUDE
                    COSINE = np.arccos(np.dot(VDIPOLE, USERVEC[NDIPOLE,1:4:1])/(np.linalg.norm(VDIPOLE)*np.linalg.norm(USERVEC[NDIPOLE,1:4:1]))) if DOTQ == 1 else 0 ## INVERSE COSINE OF THE DOT PRODUCT
                    print >> OUTPUT,  '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, DT*ODATA*CDATA*NDIPOLE), '%d' % (STATES[NDIPOLE, 1]), '%03.6f' % (MDIPOLE), '%.6f' % (np.degrees(COSINE))
                ## DELETE EXTRANEOUS DATA ##
                os.remove('%s/%04d/permdipole.out' % (NEXMD,DIR))
                os.remove('%s/%04d/pop.out' % (NEXMD,DIR))
                os.remove('%s/%04d/uservec.out' % (NEXMD,DIR))
                print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                if TSTEPS == TSMAX:
                    ETRAJ += 1
                TTRAJ += 1
        ## SUMMARY OF RESULTS ##
        if TTRAJ == 0:
            print 'No trajectories completed with %0*.2f.' % (len(str(TSMAX)), TCOLL)
        else:
            print 'Total Trajectories:', '%04d' % (TTRAJ)
            print 'Excellent Trajectories:', '%04d' % (ETRAJ)
        if ERRFLAG == 1:
            print 'One or more trajectories have experienced an error, check permdipole_raw_ensemble.err.'
        else:
            os.remove('%s/permdipole_raw_ensemble.err' % (CWD))
