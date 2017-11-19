#/usr/bin/python

'''
    
This function collects PESs and NACTs from NEXMD during adiabatic or 
non-adiabatic dynamics.

A maximum total of four output files will be generated in the current 
working directory if collecting PESs and NACTs is requested.  In 
'pes_raw_ensemble.out' and 'nact_raw_ensemble.out', the PESs and NACTs
at all time-steps or up to a time defined by the user are shown, 
respectively.  In 'pes_raw_ensemble.out', columns from left to right are: 
directory of trajectory, current state, new state, followed by all PESs.  
Likewise, 'nact_raw_ensemble.out' is in similar format.  The columns 
showing NACTs are consecutive rows of the NACT matrix, same as that 
shown in 'nact.out', located in the directory of each trajectory.  
In 'pes_hop_ensemble.out' and 'nact_hop_ensemble.out' are same data as 
'...raw_ensemble.out', but only at time-steps where hops occur.  If the 
BO flag (i.e. bo_dynamics_flag) in 'header' is set to '1', the 
simulation is adiabatic and only 'pes_raw_ensemble.out' will be 
generated.  An error file will be generated if certain files do not 
exist or if trajectories did not finish up to the user-defined length of 
analysis.

Type of calculation

[1] Ensemble of Trajectories

[1a] All time-steps
In 'pes_raw_ensemble.out':
> trajectory directory, current state, new state, PESs
In 'nact_raw_ensemble.out':
> trajectory directory, current state, new state, NACTs

[2a] All up to user-defined time
In 'pes_raw_ensemble.out':
> trajectory directory, current state, new state, PESs
In 'nact_raw_ensemble.out':
> trajectory directory, current state, new state, NACTs

[2b] Hops up to user-defined time
In 'pes_hop_ensemble.out':
> trajectory directory, current state, new state, PESs
In 'nact_hop_ensemble.out':
> trajectory directory, current state, new state, NACTs
Note: current state will not equal new state for ..._hops_...
output files

Output files:
- pes_[type].out, where [type] = hop_ensemble, raw_ensemble
- nact_[type].out, where [type] = hop_ensemble, raw_ensemble

Error files:
- pesnact_[type].out, where [type] = hop_ensemble, raw_ensemble

'''

import numpy as np
import os
import sys
import glob

CWD = os.getcwd()

def PESNACT():

    print 'Collecting PESs and/or NACTs.'

    ## TYPE OF CALCULATON AND DIRECTORY ##
    NEXMDIR = raw_input('NEXMD directory: ')
    if not os.path.exists(NEXMDIR):
        print 'Path %s does not exist.' % (NEXMDIR)
        sys.exit()
    TYPEQ = input('Output PESs and/or NACTs at all time-steps and trajectories, up to some user-defined time, or hops only?\nAnswer ALL [0], USER-DEFINED TIME [1], or HOPS [2]: ')
    if TYPEQ not in [0,1,2]:
        print 'Answer must be 0, 1, or 2.'
        sys.exit()

    ## INFORMATION FROM HEADER ##
    if not os.path.exists('%s/header' % (NEXMDIR)):
        print 'Path %s/header does not exist.' % (NEXMDIR)
        sys.exit()
    HEADER = open('%s/header'% (NEXMDIR),'r')
    HEADER = HEADER.readlines()
    NUM = 0
    TLINE = len(HEADER)
    VERB = None
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
        if 'n_quant_steps' in LINE:
            NQSTEP = np.int(LINE.split()[0][len('n_quant_steps='):-1])
            if NQSTEP == 0:
                NQSTEP = 1
        if '&moldyn' in LINE:
            TLINE = NUM
        if 'verbosity' in LINE and NUM > TLINE and VERB is None:
            VERB = np.int(LINE.split()[0][len('verbosity='):-1])
        if 'out_data_steps' in LINE:
            ODATA = np.int(LINE.split()[0][len('out_data_steps='):-1])
            if ODATA == 0:
                print 'No data has been printed to files because out_data_steps = 0 in header.'
                sys.exit()
        NUM += 1
    if BOFLAG == 1 and STATEINIT == None:
        print 'Dynamics are set to Born-Oppenheimer, but the initial state is not set.\nPlease check bo_dynamics_flag and exc_state_init in header.'
        sys.exit()
    if BOFLAG == 1 and TYPEQ == 2:
        print 'Dynamics are set to Born-Oppenheimer. Hops only occur during non-Born-Oppenheimer dynamics.\nPlease check bo_dynamics_flag in header.'
        sys.exit()
    
    ## COLLECTION TIME ##
    if TYPEQ == 0: ## ALL TIME-STEPS
        TCOLL = (TSMAX - 1)*DT
    if TYPEQ == 1 or TYPEQ == 2: ## USER-DEFINED TIME OR TIME-STEPS AT HOPS ONLY
        TCOLL = input('Collect data up to what time in femtoseconds: ')
        if isinstance(TCOLL, int) == False and isinstance(TCOLL, float) == False:
            print 'Time must be integer or float.'
            sys.exit()
        if TCOLL < 0:
            print 'Time must be integer or float greater than zero.'
            sys.exit()
        TCOLL = np.float(TCOLL)
        NSTEPS = 1
        while NSTEPS*DT <= TCOLL:
            NSTEPS += 1
        TCOLL = (NSTEPS - 1)*DT
        if TCOLL > (TSMAX - 1)*DT:
            TCOLL = (TSMAX - 1)*DT

    ## NUMBER OF CLASSICAL STEPS ##
    TSCOL = 0
    while TSCOL*DT*ODATA <= TCOLL:
        TSCOL += 1

    ## DATA TYPE ##
    if BOFLAG == 0: ## NON-ADIABATIC
        DTYPEQ = input('Collect PESs [1], NACTs [2], or BOTH [3]: ')
    if BOFLAG == 1: ## ADIABATIC
        DTYPEQ = 1
    if DTYPEQ not in [1,2,3]:
        print 'Answer must be 1, 2, or 3.'
        sys.exit()

    ## LINE NUMBER ARRAY ##
    if VERB == 3:
        if TSTEPQ == 0 and NQSTEP != 1:
            LINENUMS = np.arange(NQSTEP, TSCOL*NQSTEP - (NQSTEP - 1), NQSTEP)
        else:
            LINENUMS = np.arange(0, TSCOL)
    else:
        LINENUMS = np.arange(0, TSCOL)

    ## TIME ARRAY ##
    if DTYPEQ == 1:
        TIMES = np.around(np.linspace(TINITH, TCOLL, TSCOL), decimals = 3)
    if DTYPEQ == 2 or DTYPEQ == 3:
        TIMES = np.around(np.linspace(TINITH + DT*ODATA, TCOLL, TSCOL - 1), decimals = 3)

    ## INDICES TO CUT EXTRANEOUS NACT DATA ##
    INDICES = np.array([])
    INDEX = NSTATES
    for TERM in np.split(np.arange(NSTATES*NSTATES),NSTATES):
        INDICES = np.append(INDICES,TERM[-INDEX::])
        INDEX -= 1
    INDICES = np.int_(np.insert(INDICES + 1, 0, 0, 0))

    ## CHECK FOR NEXMD FOLDERS ##
    NEXMDS = glob.glob('%s/NEXMD*/' % (NEXMDIR))
    NEXMDS.sort()
    if len(NEXMDS) == 0:
        print 'There are no NEXMD folders in %s.' % (NEXMDIR)
        sys.exit()

    ### ADIABATIC ###
    if BOFLAG == 1:
        ## ALL TIME-STEPS ##
        if TYPEQ == 0:
            ## GENERATE OUTPUT/ERROR FILES ##
            PESALL = open('%s/pes_raw_ensemble.out' % (CWD),'w')
            ERROR = open('%s/pesnact.err' % (CWD),'w')
            ## BEGIN LOOPING OVER DIRECTORIES ##
            TTRAJ = 0
            CTRAJ = 0
            ETRAJ = 0
            ERRFLAG = 0
            for NEXMD in NEXMDS:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %dirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(DIRLIST1,int) == True:
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
                    ## DETERMINE IF PES FILES EXIST AND OPEN THEM ##
                    if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                        print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                    PES = PES.readlines()
                    LINES = LINENUMS[0:len(PES):1]
                    ## COLLECT DATA ##
                    TFLAG = 0
                    INDEX = 0
                    for LINE in LINES:
                        PESS = np.float_(PES[LINE].split())
                        TIME = np.around(PESS[0], decimals = 3)
                        if TIME != TIMES[INDEX]:
                            print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                            TFLAG = 1
                            ERRFLAG = 1
                            break
                        print >> PESALL, '%s%04d' % (NEXMD,DIR), '%d' % (STATEINIT), '%d' % (STATEINIT), ' '.join(str('%.10f') % (x) for x in PESS)
                        INDEX += 1
                    if TFLAG == 0:
                        CTRAJ += 1
                        if TSTEPS == TSMAX:
                            ETRAJ += 1
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    TTRAJ += 1
            ## SUMMARY OF RESULTS ##
            if CTRAJ == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(TSMAX)), TCOLL)
                os.remove('%s/pes_raw_ensemble.out' % (CWD))
            else:
                print 'Total Trajectories:', '%04d' % (TTRAJ)
                print 'Completed Trajectories:', '%04d' % (CTRAJ)
                print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            if ERRFLAG == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(TSMAX)),TCOLL)
            else:
                os.remove('%s/pesnact.err' % (CWD))

        ## USER-DEFINED TIME-STEPS ##
        if TYPEQ == 1:
            ## GENERATE OUTPUT/ERROR FILES ##
            PESALL = open('%s/pes_raw_ensemble.out' % (CWD),'w')
            ERROR = open('%s/pesnact.err' % (CWD),'w')
            ## BEGIN LOOPING OVER DIRECTORIES ##
            TTRAJ = 0
            CTRAJ = 0
            ETRAJ = 0
            ERRFLAG = 0
            for NEXMD in NEXMDS:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %dirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(DIRLIST1,int) == True:
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
                    ## DETERMINE IF PES FILES EXIST AND OPEN THEM ##
                    if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                        print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                    PES = PES.readlines()
                    LINES = LINENUMS[0:len(PES):1]
                    ## COMPARE COMPLETED TIME-STEPS TO COLLECTION TIME-STEPS AND COLLECT DATA ##
                    if TSTEPS >= TSCOL:
                        TFLAG = 0
                        INDEX = 0
                        for LINE in LINES:
                            PESS = np.float_(PES[LINE].split())
                            TIME = np.around(PESS[0], decimals = 3)
                            if TIME != TIMES[INDEX]:
                                print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                                TFLAG = 1
                                ERRFLAG = 1
                                break
                            print >> PESALL, '%s%04d' % (NEXMD,DIR), '%d' % (STATEINIT), '%d' % (STATEINIT), ' '.join(str('%.10f') % (x) for x in PESS)
                            INDEX += 1
                        if TFLAG == 0:
                            CTRAJ += 1
                            if TSTEPS == TSMAX:
                                ETRAJ += 1
                    else:
                        print >> ERROR, '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                        ERRFLAG = 1
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    TTRAJ += 1
            ## SUMMARY OF RESULTS ##
            if CTRAJ == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(TSMAX)), TCOLL)
                os.remove('%s/pes_raw_ensemble.out' % (CWD))
            else:
                print 'Total Trajectories:', '%04d' % (TTRAJ)
                print 'Completed Trajectories:', '%04d' % (CTRAJ)
                print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            if ERRFLAG == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(TSMAX)),TCOLL)
            else:
                os.remove('%s/pesnact.err' % (CWD))

    ### NON-ADIABATIC ###
    if BOFLAG == 0:
        ## ALL TIME-STEPS ##
        if TYPEQ == 0:
            ## GENERATE OUTPUT/ERROR FILES ##
            if DTYPEQ == 1:
                PESALL = open('%s/pes_raw_ensemble.out' % (CWD),'w')
            if DTYPEQ == 2:
                NACTALL = open('%s/nact_raw_ensemble.out' % (CWD),'w')
            if DTYPEQ == 3:
                PESALL = open('%s/pes_raw_ensemble.out' % (CWD),'w')
                NACTALL = open('%s/nact_raw_ensemble.out' % (CWD),'w')
            ERROR = open('%s/pesnact.err' % (CWD),'w')
            ## BEGIN LOOPING OVER DIRECTORIES ##
            TTRAJ = 0
            CTRAJ = 0
            ETRAJ = 0
            ERRFLAG = 0
            for NEXMD in NEXMDS:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %sdirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(DIRLIST1,int) == True:
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
                    ## DETERMINE IF PES/NACT FILES EXIST AND OPEN THEM ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,DIR)):
                        print >> ERROR, '%s%04d/coeff-n.out' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    HOPS = open('%s/%04d/coeff-n.out' % (NEXMD,DIR),'r')
                    HOPS = HOPS.readlines()
                    HSTEPS = len(HOPS)
                    if DTYPEQ == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                        PES = PES.readlines()
                        LINES = LINENUMS[0:len(PES):1]
                    if DTYPEQ == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/nact.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        NACT = open('%s/%04d/nact.out' % (NEXMD,DIR),'r')
                        NACT = NACT.readlines()
                        LINES = LINENUMS[1:len(NACT) + 1:1]
                    if DTYPEQ == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                        PES = PES.readlines()
                        LINES = LINENUMS[1:len(PES):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/nact.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        NACT = open('%s/%04d/nact.out' % (NEXMD,DIR),'r')
                        NACT = NACT.readlines()
                    ## COLLECT DATA ##
                    TFLAG = 0
                    INDEX = 0
                    CSTATE = np.int(HOPS[0].split()[0])
                    for LINE in LINES:
                        if LINE <= HSTEPS - 1:
                            NSTATE = np.int(HOPS[LINE].split()[0])
                        if DTYPEQ in [1,3]:
                            PESS = np.float_(PES[LINE].split())
                            TIME = np.around(PESS[0], decimals = 3)
                            if TIME != TIMES[INDEX]:
                                print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                                TFLAG = 1
                                ERRFLAG = 1
                                break
                        if DTYPEQ in [2,3]:
                            if LINE <= HSTEPS - 1:
                                NACTS = np.float_(NACT[LINE - 1].split())[INDICES]
                                TIME = np.around(NACTS[0], decimals = 3)
                                if TIME != TIMES[INDEX]:
                                    print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                                    TFLAG = 1
                                    ERRFLAG = 1
                                    break
                        if DTYPEQ == 1:
                            print >> PESALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                            CSTATE = NSTATE
                            INDEX += 1
                            continue
                        if DTYPEQ == 2:
                            if LINE <= HSTEPS - 1:
                                print >> NACTALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                                CSTATE = NSTATE
                                INDEX += 1
                                continue
                        if DTYPEQ == 3:
                            print >> PESALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                            if LINE <= HSTEPS - 1:
                                print >> NACTALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                        CSTATE = NSTATE
                        INDEX += 1
                    if TFLAG == 0:
                        CTRAJ += 1
                        if TSTEPS == TSMAX:
                            ETRAJ += 1
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    TTRAJ += 1
            ## SUMMARY OF RESULTS ##
            if CTRAJ == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(TSMAX)), TCOLL)
                if DTYPEQ == 1:
                    os.remove('%s/pes_raw_ensemble.out' % (CWD))
                if DTYPEQ == 2:
                    os.remove('%s/nact_raw_ensemble.out' % (CWD))
                if DTYPEQ == 3:
                    os.remove('%s/pes_raw_ensemble.out' % (CWD))
                    os.remove('%s/nact_raw_ensemble.out' % (CWD))
            else:
                print 'Total Trajectories:', '%04d' % (TTRAJ)
                print 'Completed Trajectories:', '%04d' % (CTRAJ)
                print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            if ERRFLAG == 1:
                print 'One of more trajectories have experienced an error, check pesnact.err.' % (len(str(TSMAX)), TCOLL)
            else:
                os.remove('%s/pesnact.err' % (CWD))

        ## USER-DEFINED TIME-STEPS ##
        if TYPEQ == 1:
            ## GENERATE OUTPUT/ERROR FILES ##
            if DTYPEQ == 1:
                PESALL = open('%s/pes_raw_ensemble.out' % (CWD),'w')
            if DTYPEQ == 2:
                NACTALL = open('%s/nact_raw_ensemble.out' % (CWD),'w')
            if DTYPEQ == 3:
                PESALL = open('%s/pes_raw_ensemble.out' % (CWD),'w')
                NACTALL = open('%s/nact_raw_ensemble.out' % (CWD),'w')
            ERROR = open('%s/pesnact.err' % (CWD),'w')
            ## BEGIN LOOPING OVER DIRECTORIES ##
            TTRAJ = 0
            CTRAJ = 0
            ETRAJ = 0
            ERRFLAG = 0
            for NEXMD in NEXMDS:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %sdirlist1 does not exist.' % (NEXMD)
                    sys.exit()
                DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(DIRLIST1,int) == True:
                    DIRLIST1 = np.array([DIRLIST1])
                for DIR in DIRLIST1:
                    ## DETERMINE NUMBER OF COMPLETED TIME-STEPS ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,DIR)):
                        print 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,DIR)
                        sys.exit()
                    DATA = open('%s/%04d/energy-ev.out' % (NEXMD,DIR),'r')
                    DATA = DATA.readlines()
                    TSTEPS = len(DATA) - 1
                    ## DETERMINE IF PES/NACT FILES EXIST ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,DIR)):
                        print >> ERROR, '%s%04d/coeff-n.out' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    HOPS = open('%s/%04d/coeff-n.out' % (NEXMD,DIR),'r')
                    HOPS = HOPS.readlines()
                    HTSTEPS = len(HOPS)
                    if DTYPEQ == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                        PES = PES.readlines()
                        LINES = LINENUMS[0:len(PES):1]
                    if DTYPEQ == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/nact.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        NACT = open('%s/%04d/nact.out' % (NEXMD,DIR),'r')
                        NACT = NACT.readlines()
                        LINES = LINENUMS[1:len(NACT) + 1:1]
                    if DTYPE == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                        PES = PES.readlines()
                        LINES = LINENUMS[1:len(PES):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/nact.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        NACT = open('%s/%04d/nact.out' % (NEXMD,DIR),'r')
                        NACT = NACT.readlines()
                    ## COMPARE COMPLETED TIME-STEPS TO COLLECTION TIME-STEPS AND COLLECT DATA ##
                    if TSTEPS >= TSCOL:
                        TFLAG = 0
                        INDEX = 0
                        CSTATE = np.int(HOPS[0].split()[0])
                        for LINE in LINES:
                            if LINE <= HTSTEPS - 1:
                                NSTATE = np.int(HOPS[LINE].split()[0])
                            if DTYPEQ in [1,3]:
                                PESS = np.float_(PES[LINE].split())
                                TIME = np.around(PESS[0], decimals = 3)
                                if TIME != TIMES[INDEX]:
                                    print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                                    TFLAG = 1
                                    ERRFLAG = 1
                                    break
                            if DTYPEQ in [2,3]:
                                if LINE <= HTSTEPS - 1:
                                    NACTS = np.float_(NACT[LINE-1].split())[INDICES]
                                    TIME = np.around(NACTS[0], decimals = 3)
                                    if TIME != TIMES[INDEX]:
                                        print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                                        TFLAG = 1
                                        ERRFLAG = 1
                                        break
                            if DTYPEQ == 1:
                                print >> PESALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                                CSTATE = NSTATE
                                INDEX += 1
                                continue
                            if DTYPEQ == 2:
                                if LINE <= HTSTEPS - 1:
                                    print >> NACTALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                                    CSTATE = NSTATE
                                    INDEX += 1
                                    continue
                            if DTYPEQ == 3:
                                print >> PESALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                                if LINE <= HTSTEPS - 1:
                                    print >> NACTALL, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                            CSTATE = NSTATE
                            INDEX += 1
                        if TFLAG == 0:
                            CTRAJ += 1
                            if TSTEPS == TSMAX:
                                ETRAJ += 1
                    else:
                        print >> ERROR, '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                        ERRFLAG = 1
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    TTRAJ += 1
            ## SUMMARY OF RESULTS ##
            if CTRAJ == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(TSMAX)), TCOLL)
                if DTYPEQ == 1:
                    os.remove('%s/pes_raw_ensemble.out' % (CWD))
                if DTYPEQ == 2:
                    os.remove('%s/nact_raw_ensemble.out' % (CWD))
                if DTYPEQ == 3:
                    os.remove('%s/pes_raw_ensemble.out' % (CWD))
                    os.remove('%s/nact_raw_ensemble.out' % (CWD))
            else:
                print 'Total Trajectories:', '%04d' % (TTRAJ)
                print 'Completed Trajectories:', '%04d' % (CTRAJ)
                print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            if ERRFLAG == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(TSMAX)),TCOLL)
            else:
                os.remove('%s/pesnact.err' % (CWD))

        ## TIME-STEPS AT HOPS ONLY ##
        if TYPEQ == 2:
            ## GENERATE OUTPUT/ERROR FILES ##
            if DTYPEQ == 1:
                PESHOP = open('%s/pes_hop_ensemble.out' % (CWD),'w')
            if DTYPEQ == 2:
                NACTHOP = open('%s/nact_hop_ensemble.out' % (CWD),'w')
            if DTYPEQ == 3:
                PESHOP = open('%s/pes_hop_ensemble.out' % (CWD),'w')
                NACTHOP = open('%s/nact_hop_ensemble.out' % (CWD),'w')
            ERROR = open('%s/pesnact.err' % (CWD),'w')
            ## BEGIN LOOPING OVER DIRECTORIES ##
            TTRAJ = 0
            CTRAJ = 0
            ETRAJ = 0
            ERRFLAG = 0
            for NEXMD in NEXMDS:
                if not os.path.exists('%s/dirlist1' % (NEXMD)):
                    print 'Path %sdirlist1 does not exist' % (NEXMD)
                    sys.exit()
                DIRLIST1 = np.int_(np.genfromtxt('%s/dirlist1' % (NEXMD)))
                if isinstance(DIRLIST1,int) == True:
                    DIRLIST1 = np.array([DIRLIST1])
                for DIR in DIRLIST1:
                    ## DETERMINE NUMBER OF COMPLETED TIME-STEPS ##
                    if not os.path.exists('%s/%04d/energy-ev.out' % (NEXMD,DIR)):
                        print 'Path %s%04d/energy-ev.out does not exist.' % (NEXMD,DIR)
                        sys.exit()
                    DATA = open('%s/%04d/energy-ev.out' % (NEXMD,DIR),'r')
                    DATA = DATA.readlines()
                    TSTEPS = len(DATA) - 1
                    ## DETERMINE IF PES/NACT FILES EXIST ##
                    if not os.path.exists('%s/%04d/coeff-n.out' % (NEXMD,DIR)):
                        print >> ERROR, '%s%04d/coeff-n.out' % (NEXMD,DIR), 'does not exist'
                        ERRFLAG = 1
                        TTRAJ += 1
                        continue
                    HOPS = open('%s/%04d/coeff-n.out' % (NEXMD,DIR),'r')
                    HOPS = HOPS.readlines()
                    HSTEPS = len(HOPS)
                    if DTYPEQ == 1:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                        PES = PES.readlines()
                        LINES = LINENUMS[0:len(PES):1]
                    if DTYPEQ == 2:
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/nact.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        NACT = open('%s/%04d/nact.out' % (NEXMD,DIR),'r')
                        NACT = NACT.readlines()
                        LINES = LINENUMS[1:len(NACT) + 1:1]
                    if DTYPEQ == 3:
                        if not os.path.exists('%s/%04d/pes.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/pes.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        PES = open('%s/%04d/pes.out' % (NEXMD,DIR),'r')
                        PES = PES.readlines()
                        LINES = LINENUMS[1:len(PES):1]
                        if not os.path.exists('%s/%04d/nact.out' % (NEXMD,DIR)):
                            print >> ERROR, '%s%04d/nact.out' % (NEXMD,DIR), 'does not exist'
                            ERRFLAG = 1
                            TTRAJ += 1
                            continue
                        NACT = open('%s/%04d/nact.out' % (NEXMD,DIR),'r')
                        NACT = NACT.readlines()
                    ## COMPARE COMPLETED TIME-STEPS TO COLLECTION TIME-STEPS AND COLLECT DATA ##
                    if TSTEPS >= TSCOL:
                        TFLAG = 0
                        CSTATE = np.int(HOPS[0].split()[0])
                        INDEX = 0
                        for LINE in LINES:
                            if LINE <= HSTEPS - 1:
                                NSTATE = np.int(HOPS[LINE].split()[0])
                            if DTYPEQ in [1,3]:
                                PESS = np.float_(PES[LINE].split())
                                TIME = np.around(PESS[0], decimals = 3)
                                if TIME != TIMES[INDEX]:
                                    print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                                    TFLAG = 1
                                    ERRFLAG = 1
                                    break
                            if DTYPEQ in [2,3]:
                                if LINE <= HSTEPS - 1:
                                    NACTS = np.float_(NACT[LINE - 1].split())[INDICES]
                                    TIME = np.around(NACTS[0], decimals = 3)
                                    if TIME != TIMES[INDEX]:
                                        print >> ERROR, 'There is an inconsistency in time-step in %s%04d at %.3f fs' % (NEXMD,DIR,TIMES[INDEX])
                                        TFLAG = 1
                                        ERRFLAG = 1
                                        break
                            if NSTATE != CSTATE:
                                if DTYPEQ == 1:
                                    print >> PESHOP, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                                    CSTATE = NSTATE
                                    INDEX += 1
                                    continue
                                if DTYPEQ == 2:
                                    if LINE <= HSTEPS - 1:
                                        print >> NACTHOP, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                                        CSTATE = NSTATE
                                        INDEX += 1
                                        continue
                                if DTYPEQ == 3:
                                    print >> PESHOP, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in PESS)
                                    if LINE <= HSTEPS - 1:
                                        print >> NACTHOP, '%s%04d' % (NEXMD,DIR), '%d' % (CSTATE), '%d' % (NSTATE), ' '.join(str('%.10f') % (x) for x in NACTS)
                            CSTATE = NSTATE
                            INDEX += 1
                        if TFLAG == 0:
                            CTRAJ += 1
                            if TSTEPS == TSMAX:
                                ETRAJ += 1
                    else:
                        print >> ERROR, '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                        ERRFLAG = 1
                    print '%s%04d' % (NEXMD,DIR), '%0*.2f' % (len(str((TSMAX))) + 2, (TSTEPS - 1)*DT)
                    TTRAJ += 1
            ## SUMMARY OF RESULTS ##
            if CTRAJ == 0:
                print 'No trajectories completed within %0*.2f fs.' % (len(str(TSMAX)), TCOLL)
                if DTYPEQ == 1:
                    os.remove('%s/pes_hop_ensemble.out' % (CWD))
                if DTYPEQ == 2:
                    os.remove('%s/nact_hop_ensemble.out' % (CWD))
                if DTYPEQ == 3:
                    os.remove('%s/pes_hop_ensemble.out' % (CWD))
                    os.remove('%s/nact_hop_ensemble.out' % (CWD))
            else:
                print 'Total Trajectories:', '%04d' % (TTRAJ)
                print 'Completed Trajectories:', '%04d' % (CTRAJ)
                print 'Excellent Trajectories:', '%04d' % (ETRAJ)
            if ERRFLAG == 1:
                print 'One or more trajectories did not finish within %0*.2f femtoseconds, check pesnact.err.' % (len(str(TSMAX)),TCOLL)
            else:
                os.remove('%s/pesnact.err' % (CWD))
