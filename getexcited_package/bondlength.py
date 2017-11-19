#/usr/bin/python

'''

This function calculates a user-specified bond length in time.

The bond length is calculated between two atoms at every time-step
and may be tracked along a single trajectory or an ensemble of
trajectories during adiabatic or non-adiabatic dynamics.  The
user must supply the line numbers of two atoms, where the first atom
between '&coord' and '&endcoord' in 'input.ceon' is always labeled with
line number = 0.

Type of calculation:

[1] Single Trajectory
< collection time (fs)
> time (fs), bond length (Angstroms)

[2] Ensemble of Trajectories

[2a] Mean
< collection time (fs)
> time (fs), bond length (Angstroms), standard deviation (Angstroms)

[2b] All time-steps
> trajectory directory, bond length (Angstroms) at all time-steps and
trajectories

Output Files:
- bl_[type].out, where [type] = single, mean_ensemble, raw_ensemble

Error Files:
- bl_[type].err, where [type] = single, mean_ensemble, raw_ensemble

'''

import numpy as np
import os
import sys
import glob
import fileinput
import math

cwd = os.getcwd()

def BONDLENGTH():

    print 'Calculating a bond length as a function of time.'

    ## Type of calculation and directory check ##
    dynq = input('Calculate bond length along one trajectory or an ensemble of trajectories?\nAnswer one [1] or ensemble [0]: ')
    if dynq not in [1,0]:
        print 'Answer must be 1 or 0.'
        sys.exit()
    if dynq == 0: ## ensemble
        nexmdir = raw_input('Ensemble directory [e.g. nexmd]: ')
        if not os.path.exists(nexmdir):
            print 'Path %s does not exist.' % (nexmdir)
            sys.exit()
        ## Check if nexmd folders exist ##
        nexmds = glob.glob('%s/nexmd*/' % (nexmdir))
        nexmds.sort()
        if len(nexmds) == 0:
            print 'There are no nexmd folders in %s.' % (nexmdir)
            sys.exit()
        ## determine mean or all ##
        typeq = input('Output mean bla in time or output bla at all time-steps and trajectories?\nAnswer mean [0] or all [1]: ')
        if typeq not in [0,1]:
            print 'Answer must be 0 or 1.'
            sys.exit()
    if dynq == 1: ## single trajectory
        typeq = 0
        nexmdir = raw_input('Single trajectory directory: ')
        if not os.path.exists(nexmdir):
            print 'Path %s does not exist.' % (nexmdir)
            sys.exit()

    ## Information from header ##
    if dynq == 0: ## ensemble
        if not os.path.exists('%s/header' % (nexmdir)):
            print 'Path %s/header does not exist.' % (nexmdir)
            sys.exit()
        header = open('%s/header' % (nexmdir),'r')
        header = header.readlines()
    if dynq == 1: ## single trajectory
        if not os.path.exists('%s/input.ceon' % (nexmdir)):
            print 'Path %s/input.ceon does not exist.' % (nexmdir)
            sys.exit()
        header = open('%s/input.ceon' % (nexmdir),'r')
        header = header.readlines()
    for line in header:
        if 'time_init' in line:
            tinith = np.float(line.split()[0][len('time_init='):-1])
        if 'time_step' in line:
            dt = np.float(line.split()[0][len('time_step='):-1])
        if 'n_class_steps' in line:
            tsmax = np.int(line.split()[0][len('n_class_steps='):-1]) + 1
        if 'out_data_steps' in line:
            odata = np.int(line.split()[0][len('out_data_steps='):-1])
        if 'out_coords_steps' in line:
            cdata = np.int(line.split()[0][len('out_coords_steps='):-1])
        if 'natoms' in line:
            natoms = np.int(line.split()[0][len('natoms='):-1])

    ## Collection time ##
    if typeq == 0: ## mean bond length
        if dynq == 0: ## ensemble
            tcoll = input('Calculate bond length up to what time in femtoseconds?\nNote that averaged results will only include trajectories that are complete up to this time: ')
        if dynq == 1: ## single trajectory
            tcoll = input('Calculate bond length up to what time in femtoseconds? ')
        if isinstance(tcoll, int) == false and isinstance(tcoll, float) == false:
            print 'Time must be integer or float.'
            sys.exit()
        if tcoll < 0:
            print 'Time must be integer or float greater than zero.'
            sys.exit()
        tcoll = np.float(tcoll)
        if tcoll > (tsmax - 1)*dt:
            tcoll = (tsmax - 1)*dt
    if typeq == 1: ## all bla
        tcoll = (tsmax - 1)*dt

    ## Number of classical steps ##
    tscol = 0
    while tscol*dt*odata <= tcoll:
        tscol += 1

    ## Number of time-steps for coordinates ##
    ccoll = 0
    num = 0
    while ccoll <= tcoll:
        ccoll += dt*odata*cdata
        num += 1

    ## Collection time ##
    times = np.linspace(tinith, ccoll - dt*odata*cdata, num)

    ## Two unique atoms defined by user ##
    lines = input('Input the line numbers labeling the coordinates of the two atoms.\nInput an array of the form [[atom1, atom2], [atom3, atom4], .. ]: ')
    for line in lines:
        if isinstance(line, list) == false:
            print 'Subarray must be of the form [atom1, atom2], where atom# = line number of atom#.'
            sys.exit()
        if len(line) != 2:
            print 'Subarray must contain two elements labeling the line numbers of two atoms.'
            sys.exit()
        index = 0
        for i in line:
            if isinstance(i, int) == false:
                print 'Element number %d of subarray must be integer.\nuser inputted [%s, %s], which is not allowed.' % (index + 1, line[0], line[1])
                sys.exit()
            if i < 0:
                print 'Element number %d of subarray must be a positive integer.\nuser inputted [%s, %s], which is not allowed.' % (index + 1, line[0], line[1])
                sys.exit()
            if i > natoms - 1:
                print 'Element number %d of subarray must be less than the max number of atoms (-1).\nuser inputted [%s, %s], which is not allowed.' % (index + 1, line[0], line[1])
                sys.exit()
            index += 1
        if len(np.unique(line)) != 2:
            print 'All elements of subarray must be unique.\nuser inputted [%s, %s], which is not allowed.' % (line[0], line[1])
            sys.exit()
    nbonds = len(lines)

    ## Calculate bond length along a single trajectory ##
    if dynq == 1: ## single trajectory
        print 'Collecting bond length along single trajectory.  Please wait ...'
        ## genrate output file ##
        output = open('%s/bl_single.out' % (cwd),'w')
        etraj = 0
        ## Determine completed number of time-steps ##
        if not os.path.exists('%s/energy-ev.out' % (nexmdir)):
            print 'Path %s/energy-ev.out does not exist.' % (nexmdir)
            sys.exit()
        data = open('%s/energy-ev.out' % (nexmdir),'r')
        data = data.readlines()
        tsteps = len(data) - 1
        ## Generate array with indices of the coordinate blocks along a trajectory ##
        if tsteps >= tscol:
            if not os.path.exists('%s/coords.xyz' % (nexmdir)):
                print 'path %s/coords.xyz does not exist.' % (nexmdir)
                sys.exit()
            data = open('%s/coords.xyz' % (nexmdir),'r')
            data = data.readlines()
            lenc = len(data)
            ncoords = 0
            cindex = 0
            tflag1 = 0
            tflag2 = 0
            array = np.array([])
            for line in data:
                if 'time' in line:
                    if ncoords == 0:
                        tinit = np.float(line.split()[-1])
                        if tinit != tinith:
                            tflag1 = 1
                            break
                    else:
                        time = np.around(np.float(line.split()[-1]), decimals = 3)
                        if time > tcoll:
                            tflag3 = 1
                            break
                        if time != times[ncoords]:
                            tflag2 = 1
                            break
                    ncoords += 1
                    array = np.append(array,cindex)
                cindex += 1
            if tflag1 == 1:
                print 'Initial time in %s/coords.xyz does not match time_init in %s/input.ceon.' % (nexmdir,nexmdir)
                sys.exit()
            if tflag2 == 1:
                print 'There is an inconsistency in time-step in %s/coords.xyz.' % (nexmdir)
                sys.exit()
            ## Append lines for last coordinate set ##
            if tflag3 == 1:
                array = np.append(array,cindex)
            else:
                array = np.append(array, lenc + 1)
            array = np.int_(array)
            ## Checks to ensure bond length calculation ##
            if ncoords == 0:
                print 'No coordinates were found in %s/coords.xyz' % (nexmdir)
                sys.exit()
            if ncoords == 1:
                print 'Only initial coordinates, at %.2f fs, were found in %s/coords.xyz.' % (tinit,nexmdir)
                sys.exit()
            ## Calculate bond length along a single trajectory ##
            sbondlen = np.zeros((ncoords,nbonds))
            for ncoord in np.arange(ncoords):
                coords = data[array[ncoord]+1:array[ncoord+1]-1:1]
                index = 0
                for line in lines:
                    vec0 = np.float_(coords[line[0]].split()[1:])
                    vec1 = np.float_(coords[line[1]].split()[1:])
                    a = np.subtract(vec1, vec0)
                    sbondlen[ncoord,index] = np.linalg.norm(a)
                    index += 1
            print '%s' % (nexmdir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
            ctraj = 1
            if tsteps == tsmax:
                etraj = 1
        else:
            print '%s' % (nexmdir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
        ttraj = 1
        ## Summary of results ##
        if ctraj == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(tsmax)),tcoll)
        else:
            print 'Total trajectories:', '%04d' % (ttraj)
            print 'Completed trajectories:', '%04d' % (ctraj)
            print 'Excellent trajectories:', '%04d' % (etraj)
            print >> output, 'Total trajectories: ', '%04d' % (ttraj)
            print >> output, 'Completed trajectories: ', '%04d' % (ctraj)
            print >> output, 'Excellent trajectories: ', '%04d' % (etraj)
            for ncoord in np.arange(ncoords):
                print >> output, '%0*.2f' % (len(str((tsmax))) + 2,dt*odata*cdata*ncoord), ' '.join('%08.3f' % (bond) for bond in sbondlen[ncoord])

    ## Calculate bond length along an ensemble of trajectories ##
    if dynq == 0 and typeq == 0: ## mean from ensemble
        print 'Collecting mean bond length from ensemble.  Please wait ...'
        ## Determine total number of trajectories in ensemble ##
        with open('%s/totdirlist' % (nexmdir),'w') as data:
            for nexmd in nexmds:
                if not os.path.exists('%s/dirlist1' % (nexmd)):
                    print 'Path %nexmdirlist1 does not exist.' % (nexmd)
                    sys.exit()
                input = fileinput.input('%s/dirlist1' % (nexmd))
                data.writelines(input)
        dirlist1 = np.int_(np.genfromtxt('%s/totdirlist' % (nexmdir)))
        if isinstance(dirlist1,int) == true:
            dirlist1 = np.array([dirlist1])
        os.remove('%s/totdirlist' % (nexmdir))
        ## Generate output and error files ##
        output = open('%s/bl_mean_ensemble.out' % (cwd),'w')
        error = open('%s/bl_mean_ensemble.err' % (cwd),'w')
        ## Generate bond length array for final results ##
        fbondlen = np.zeros((len(times), nbonds))
        ebondlen = np.zeros((len(times), len(dirlist1), nbonds))
        ttraj = 0
        ctraj = 0
        etraj = 0
        errflag = 0
        for nexmd in nexmds:
            if not os.path.exists('%s/dirlist1' % (nexmd)):
                print 'Path %dirlist1 does not exist.' % (nexmd)
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (nexmd)))
            if isinstance(dirlist1, int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Determine completed number of time-steps ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (nexmd,dir)):
                    print >> error, '%s%04d/energy-ev.out' % (nexmd,dir), 'does not exist'
                    errflag = 1
                    ttraj += 1
                    continue
                data = open('%s/%04d/energy-ev.out' % (nexmd,dir),'r')
                data = data.readlines()
                tsteps = len(data) - 1
                ## Generate array with indices of the coordinate blocks along trajectory ##
                if tsteps >= tscol:
                    if not os.path.exists('%s/%04d/coords.xyz' % (nexmd,dir)):
                        print >> error, '%s%04d/coords.xyz' % (nexmd,dir), 'does not exist'
                        errflag = 1
                        ttraj += 1
                        continue
                    data = open('%s/%04d/coords.xyz' % (nexmd,dir),'r')
                    data = data.readlines()
                    lenc = len(data)
                    ncoords = 0
                    cindex = 0
                    tflag1 = 0
                    tflag2 = 0
                    tflag3 = 0
                    array = np.array([])
                    for line in data:
                        if 'time' in line:
                            if ncoords == 0:
                                tinit = np.float(line.split()[-1])
                                if tinit != tinith:
                                    tflag1 = 1
                                    continue
                            else:
                                time = np.around(np.float(line.split()[-1]), decimals = 3)
                                if time > tcoll:
                                    tflag3 = 1
                                    continue
                                if time != times[ncoords]:
                                    tflag2 = 1
                                    continue
                            ncoords += 1
                            array = np.append(array,cindex)
                        cindex += 1
                    if tflag1 == 1:
                        print >> error, 'Initial time in %s%04d/coords.xyz does not match time_init in %s/header.' % (nexmd,dir,nexmdir)
                        errflag = 1
                        ttraj += 1
                        continue
                    if tflag2 == 1:
                        print >> error, 'There is an inconsistency in time-step in %s%04d/coords.xyz.' % (nexmd,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Append lines for last coordinate set ##
                    if tflag3 == 1:
                        array = np.append(array,cindex)
                    else:
                        array = np.append(array, lenc + 1)
                    array = np.int_(array)
                    ## Checks to ensure bond length calculation ##
                    if ncoords == 0:
                        print >> error, 'No coordinates were found in %s%04d/coords.xyz' % (nexmd,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    if ncoords == 1:
                        print >> error, 'Only initial coordinates, at %.2f fs, were found in %s%04d/coords.xyz.' % (tinit,nexmd,dir)
                        errflag = 1
                        ttraj += 1
                        continue
                    ## Calculate bond length along a single trajectory ##
                    sbondlen = np.zeros((ncoords, nbonds))
                    for ncoord in np.arange(ncoords):
                        coords = data[array[ncoord]+1:array[ncoord+1]-1:1]
                        index = 0
                        for line in lines:
                            vec0 = np.float_(coords[line[0]].split()[1:])
                            vec1 = np.float_(coords[line[1]].split()[1:])
                            a = np.subtract(vec1, vec0)
                            sbondlen[ncoord,index] = np.linalg.norm(a)
                            ebondlen[ncoord,ctraj,index] = sbondlen[ncoord,index]
                            index += 1
                    fbondlen += sbondlen
                    print '%s%04d' % (nexmd,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    ctraj += 1
                    if tsteps == tsmax:
                        etraj += 1
                else:
                    print '%s%04d' % (nexmd,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    print >> error, '%s%04d' % (nexmd,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    errflag = 1
                ttraj += 1
        ## Summary of results ##
        if ctraj == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(tsmax)),tcoll)
        else:
            ## Mean and standard deviation for bond length ##
            ebondlen = np.delete(ebondlen, np.arange(ctraj, ttraj), axis = 1)
            ebondlen = np.std(ebondlen, axis = 1)
            fbondlen = fbondlen/ctraj
            print 'Total trajectories:', '%04d' % (ttraj)
            print 'Completed trajectories:', '%04d' % (ctraj)
            print 'Excellent trajectories:', '%04d' % (etraj)
            print >> output, 'Total trajectories: ', '%04d' % (ttraj)
            print >> output, 'Completed trajectories: ', '%04d' % (ctraj)
            print >> output, 'Excellent trajectories: ', '%04d' % (etraj)
            for ncoord in np.arange(ncoords):
                print >> output, '%0*.2f' % (len(str((tsmax))) + 2,dt*odata*cdata*ncoord), ' '.join('%08.3f' % (bond) for bond in fbondlen[ncoord]), ' '.join('%07.3f' % (bond) for bond in ebondlen[ncoord])
        if errflag == 1:
            print 'One or more trajectories did not finish within %0*.2f femtoseconds, check bl_mean_ensemble.err.' % (len(str(tsmax)),tcoll)
        else:
            os.remove('%s/bl_mean_ensemble.err' % (cwd))

    ## Calculate bond length from ensemble of trajectories at all time-steps ##
    if dynq == 0 and typeq == 1: ## all from ensemble
        print 'Collecting all bond lengths from ensemble.  Please wait ...'
        ## Generate output and error files ##
        output = open('%s/bl_raw_ensemble.out' % (cwd),'w')
        error = open('%s/bl_raw_ensemble.err' % (cwd),'w')
        ttraj = 0
        etraj = 0
        errflag = 0
        for nexmd in nexmds:
            if not os.path.exists('%s/dirlist1' % (nexmd)):
                print 'Path %dirlist1 does not exist.' % (nexmd)
                sys.exit()
            dirlist1 = np.int_(np.genfromtxt('%s/dirlist1' % (nexmd)))
            if isinstance(dirlist1, int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                ## Determine number of time-steps completed ##
                if not os.path.exists('%s/%04d/energy-ev.out' % (nexmd,dir)):
                    print >> error, '%s%04d/energy-ev.out' % (nexmd,dir), 'does not exist'
                    errflag = 1
                    ttraj += 1
                    continue
                data = open('%s/%04d/energy-ev.out' % (nexmd,dir),'r')
                data = data.readlines()
                tsteps = len(data) - 1
                ## Generate array with indices of the coordinate blocks along trajectory ##
                if not os.path.exists('%s/%04d/coords.xyz' % (nexmd,dir)):
                    print >> error, '%s%04d/coords.xyz' % (nexmd,dir), 'does not exist'
                    errflag = 1
                    ttraj += 1
                    continue
                data = open('%s/%04d/coords.xyz' % (nexmd,dir),'r')
                data = data.readlines()
                lenc = len(data)
                ncoords = 0
                cindex = 0
                tflag1 = 0
                tflag2 = 0
                tflag3 = 0
                array = np.array([])
                for line in data:
                    if 'time' in line:
                        if ncoords == 0:
                            tinit = np.float(line.split()[-1])
                            if tinit != tinith:
                                tflag1 = 1
                                continue
                        else:
                            time = np.around(np.float(line.split()[-1]), decimals = 3)
                            if time > tcoll:
                                tflag3 = 1
                                continue
                            if time != times[ncoords]:
                                tflag2 = 1
                                continue
                        ncoords += 1
                        array = np.append(array,cindex)
                    cindex += 1
                if tflag1 == 1:
                    print >> error, 'Initial time in %s%04d/coords.xyz does not match time_init in %s/header.' % (nexmd,dir,nexmdir)
                    errflag = 1
                    ttraj += 1
                    continue
                if tflag2 == 1:
                    print >> error, 'There is an inconsistency in time-step in %s%04d/coords.xyz.' % (nexmd,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                ## Append lines for last coordinate set ##
                if tflag3 == 1:
                    array = np.append(array,cindex)
                else:
                    array = np.append(array, lenc + 1)
                array = np.int_(array)
                ## Checks to ensure bond length calculation ##
                if ncoords == 0:
                    print >> error, 'No coordinates were found in %s%04d/coords.xyz' % (nexmd,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                if ncoords == 1:
                    print >> error, 'Only initial coordinates, at %.2f fs, were found in %s%04d/coords.xyz.' % (tinit,nexmd,dir)
                    errflag = 1
                    ttraj += 1
                    continue
                ## Calculate bond length along a single trajectory ##
                for ncoord in np.arange(ncoords):
                    coords = data[array[ncoord]+1:array[ncoord+1]-1:1]
                    sbondlen = np.zeros(nbonds)
                    index = 0
                    for line in lines:
                        vec0 = np.float_(coords[line[0]].split()[1:])
                        vec1 = np.float_(coords[line[1]].split()[1:])
                        a = np.subtract(vec1, vec0)
                        sbondlen[index] = np.linalg.norm(a)
                        index += 1
                    print >> output, '%s%04d' % (nexmd,dir), '%0*.2f' % (len(str((tsmax))) + 2,dt*odata*cdata*ncoord), ' '.join('%08.3f' % (bond) for bond in sbondlen)
                print '%s%04d' % (nexmd,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                if tsteps == tsmax:
                    etraj += 1
                else:
                    print '%s%04d' % (nexmd,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    print >> error, '%s%04d' % (nexmd,dir), '%0*.2f' % (len(str((tsmax))) + 2, (tsteps - 1)*dt)
                    errflag = 1
                ttraj += 1
        ## Summary of results ##
        if ttraj == 0:
            print 'No trajectories completed within %0*.2f.' % (len(str(tsmax)),tcoll)
        else:
            print 'Total trajectories:', '%04d' % (ttraj)
            print 'Excellent trajectories:', '%04d' % (etraj)
        if errflag == 1:
            print 'One or more trajectories have experienced an error, check bl_raw_ensemble.err.'
        else:
            os.remove('%s/bl_raw_ensemble.err' % (cwd))
