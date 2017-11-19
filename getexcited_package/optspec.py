#/usr/bin/python

'''
 ___________________________________________________________________
|                                                                   |
| This function generates an optical spectrum from single-point     |
| calculations.                                                     |
|                                                                   |
| One of two files will be generated in the current working         |
| directory if an optical spectrum is requested, 'ceo_<type>.out'   |
| or 'ceo.err'.  Here, <type> is 'gauss' or 'lorentz'. These refer  |
| to Gaussian or Lorentzian lineshapes, respectively.  In           |
| 'ceo.out', the spectrum is generated with one of these            |
| lineshapes.  Both the lineshape and the width are defined by the  |
| user.  The outputted spectrum is a sum of all spectra determined  |
| from initial geometries, divided the number of geometries.  The   |
| 'collectceo.sh' script, located in the getexcited_package is      |
| called to grep excitation energies and oscillator strengths from  |
| single-point calculations.  In 'ceo.out', first column is energy  |
| in eV, followed by the spectrum associated to each state,         |
| followed by the total spectrum.  Spectra should be interpreted as |
| relative absorbance in arbitrary units.  In 'ceo.err',            |
| directories of single-point calculations that are incomplete are  |
| listed. The 'ceo.out' file cannot be generated unless all         |
| single-point calculations are complete.                           |
|___________________________________________________________________|

'''

import numpy as np
import os
import sys
import subprocess
import shlex
import glob

cwd = os.getcwd()

def OPTSPEC(pathtoceo):

    print 'Generating optical spectrum.'

    ## Directory names ##
    spdir = raw_input('Single-point calculations directory: ')
    if not os.path.exists(spdir):
        print 'Path %s does not exist.' % (spdir)
        sys.exit()
    
    ## Check excitation energies ##
    print 'Checking energies and oscillator strengths. Please wait ...'
    if not os.path.exists('%s/getexcited_package/collectceo.sh' % (pathtoceo)):
        print 'The script, collectceo.sh, must be in the getexcited_package.'
        sys.exit()
    if not os.path.exists('%s/header' % (spdir)):
        print 'Path %s/header does not exist.' % (spdir)
        sys.exit()
    header = open('%s/header'% (spdir),'r')
    for line in header:
        if 'n_exc_states_propagate' in line:
            nstates = np.int(line.split()[0][len('n_exc_states_propagate='):-1])
            break
    NEXMDs = glob.glob('%s/NEXMD*/' % (spdir))
    NEXMDs.sort()
    if len(NEXMDs) == 0:
        print 'There are no NEXMD folders in %s.' % (spdir)
        sys.exit()
    error = open('%s/ceo.err' % (cwd),'w')
    ceoflag = 0
    for NEXMD in NEXMDs:
        if not os.path.exists('%s/%s/dirlist1' % (cwd,NEXMD)):
            print 'Path %s/%sdirlist1 does not exist.' % (cwd,NEXMD)
            sys.exit()
        dirlist1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (cwd,NEXMD)))
        if isinstance(dirlist1,int) == true:
            dirlist1 = np.array([dirlist1])
        for dir in dirlist1:
            if not os.path.exists('%s/%s/%04d' % (cwd,NEXMD,dir)):
                print >> error, '%s%04d' % (NEXMD,dir), 'does not exist'
                ceoflag = 1
                continue
            os.chdir('%s/%s/%04d' % (cwd,NEXMD,dir))
            if not os.path.exists('%s/%s/%04d/md.out' % (cwd,NEXMD,dir)):
                print >> error, '%s%04d/md.out' % (NEXMD,dir), 'does not exist'
                ceoflag = 1
                continue
            subprocess.call(shlex.split('sh %s/getexcited_package/collectceo.sh %d' % (pathtoceo,nstates+1)))
            if not os.path.exists('%s/%s/%04d/ceo.out' % (cwd,NEXMD,dir)):
                print >> error, '%s/%04d/ceo.out' % (NEXMD,dir), 'does not exist'
                ceoflag = 1
                continue
            with open('%s/%s/%04d/ceo.out' % (cwd,NEXMD,dir),'r') as data:
                if len(data.readlines()) != nstates:
                    print >> error, '%s%04d/ceo.out' % (NEXMD,dir), 'is incomplete'
                    ceoflag = 1
    if ceoflag == 1:
        print 'One or more single-point calculations did not finish, check ceo.err.'
        sys.exit()
    else:
        os.remove('%s/ceo.err' % (cwd))

    ## Generate optical spectrum ##
    stype = input('Spectral lineshape? Answer Gaussian [0] or Lorentzian [1]: ')
    if stype not in [0,1]:
        print 'Answer must be 0 or 1.'
        sys.exit()
    if stype == 0:
        specb = input('Spectral broadening (i.e. Gaussian standard deviation) in eV [e.g. 0.15]: ')
    else:
        specb = input('Spectral broadening (i.e. Lorentzian fwhm) in eV [e.g. 0.36]: ')
    if isinstance(specb, int) == false and isinstance(specb, float) == false:
        print 'Spectral broadening must be integer or float.'
        sys.exit()
    if specb < 0:
        print 'Spectral broadening must be integer or float greater than zero.'
        sys.exit()
    npoints = 100000
    dpoints = np.zeros((nstates+2,npoints))
    dpoints[0] = np.linspace(0.0,10.0,npoints)
    traj = 0
    if stype == 0:
        for NEXMD in NEXMDs:
            dirlist1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (cwd,NEXMD)))
            if isinstance(dirlist1,int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                data = np.genfromtxt('%s/%s/%04d/ceo.out' % (cwd,NEXMD,dir))
                for state, energy, osx, osy, osz, oscstren in data:
                    dpoints[state] = dpoints[state] + oscstren*np.exp(-(energy-dpoints[0])**(2.0)/(2.0*specb**(2.0)))/np.sqrt(2.0*np.pi*specb**(2.0))
                print '%s%04d' % (NEXMD,dir)
                traj += 1
    else:
        for NEXMD in NEXMDs:
            dirlist1 = np.int_(np.genfromtxt('%s/%s/dirlist1' % (cwd,NEXMD)))
            if isinstance(dirlist1,int) == true:
                dirlist1 = np.array([dirlist1])
            for dir in dirlist1:
                data = np.genfromtxt('%s/%s/%04d/ceo.out' % (cwd,NEXMD,dir))
                for state, energy, osx, osy, osz, oscstren in data:
                    dpoints[state] = dpoints[state] + oscstren/(1.0+((energy-dpoints[0])**2.0)/(specb/2.0)**(2.0))/(specb*np.pi/2.0)
                print '%s%04d' % (NEXMD,dir)
                traj += 1
    dpoints[0] = dpoints[0]*traj
    dpoints[-1] = np.sum(dpoints[1:-1:1], 0)
    np.savetxt('%s/ceo_%s.out' % (cwd,'gauss' if stype == 0 else 'lorentz'), np.transpose(dpoints/traj), fmt='%10.5e')
