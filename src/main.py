import numpy as np
from subprocess import call
import os
from pathlib import Path

from .tracing import trace
from .musictools import musicgadget3
from .read_config import getproperties


def createIC(haloid, parentlabel, zoominlevel=11, halolistfolder='./halotracing', skiptracing=True,
             outputdir='default'):
    """
    Script to trace a halo in the parent simulation and creating the zoom-in ICs for it.
    You can set additional parameters for the single routines are in the code.

     haloid         : ID of the halo for which you want to creat the zoom-in IC.
     parentlabel    : Label of the parent simulation, to be used in getproperties().
     zoominlevel    : Maximum resolution level of the zoom-in simulation.
     halolistfolder : Folder containing the halolist.txt file and the halo tracing files.
     skiptracing    : If True, the tracing is skipped when there is already a position list for the specified halo.
    """
    if outputdir == "default":
        outputdir = str(Path(__file__).parent) +"/../"+ parentlabel
        if not os.path.exists(outputdir):
            os.mkdir(outputdir)
    poslistdir, snapfilebase, ICfile, cosmology, boxsize, zstart, seedsset, parentres, lowestres = getproperties(
        parentlabel)
    if os.path.isfile(poslistdir + "/poslist" + str(haloid) + ".txt") and skiptracing:
        print("Position list file for this halo already present; skipping halo tracing.")
    else:
        trace(haloid, parentlabel=parentlabel, halolocfile=halolistfolder + '/halo' + str(haloid) + '.txt',
              outputdir=outputdir)
        call(["cp", "poslist" + str(haloid) + ".txt", poslistdir])

    musicgadget3(haloid, parentlabel=parentlabel, highestres=zoominlevel, initialpad=6, regionmode='ellipsoid',
                 outputdir=outputdir)


def multipleICcreator(parentlabel, zoominlevel=11, halolistfolder='./halotracing', startingatnumber=0,
                      skiptracing=True, outputdir='default'):
    """
    Creates ICs for all the halos listed in a halo list file.

    Parameters:
     parentlabel      : Label of the parent simulation, to be used in getproperties().
     zoominlevel      : Maximum resolution level of this series of zoom-in ICs.
     halolistfolder     : File with the list of halos (as generated by rockstarselect()).
     startingatnumber : Ignores the halos in the list with id smaller than this number.
     skiptracing      : If True, the tracing is skipped when there is already a position list for the specified halo.
    """
    call(['mkdir', 'ZoomICsLevel' + str(zoominlevel)])
    haloids = np.loadtxt(halolistfolder + '/halolist.txt')[:, 0]
    for i in haloids:
        if i >= startingatnumber:
            print('\n \n Constructing IC number: ' + str(i))
            try:
                createIC(int(i), parentlabel, zoominlevel, halolistfolder, skiptracing, outputdir)
                call(['mv', 'IC' + str(int(i)) + '.gdt', 'ZoomICsLevel' + str(zoominlevel)])
            except MemoryError:
                print("Could not create IC for this halo. Possibly a memory error in the tracing.")