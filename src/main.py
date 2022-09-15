import numpy as np
from subprocess import call
import os
from tracing import trace
from musictools import musicgadget3
import yaml
from pathlib import Path


def createIC(haloid, parentlabel, zoominlevel=11, halolistfolder='./halotracing', skiptracing=True, outputdir='default'):
    """
    Script to trace a halo in the parent simulation and creating the zoom-in ICs for it.
    You can set additional parameters for the single routines are in the code.

     haloid         : ID of the halo for which you want to creat the zoom-in IC.
     parentlabel    : Label of the parent simulation, to be used in getproperties().
     zoominlevel    : Maximum resolution level of the zoom-in simulation.
     halolistfolder : Folder containing the halolist.txt file and the halo tracing files.
     skiptracing    : If True, the tracing is skipped when there is already a position list for the specified halo.
    """
    poslistdir, snapfilebase, ICfile, cosmology, boxsize, zstart, seedsset, parentres, lowestres = getproperties(
        parentlabel)
    if os.path.isfile(poslistdir + "/poslist" + str(haloid) + ".txt") and skiptracing:
        print("Position list file for this halo already present; skipping halo tracing.")
    else:
        trace(haloid, parentlabel=parentlabel, halolocfile=halolistfolder + '/halo' + str(haloid) + '.txt')
        call(["cp", "poslist" + str(haloid) + ".txt", poslistdir])
    if outputdir == "default":
        outputdir = Path(__file__).parent / "../"+parentlabel
        os.mkdir(outputdir)

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


def getproperties(parentlabel):
    """
    Function that reads the properties of the parent simulations in use from config.yaml.
    For a given parent label, it returns:

     poslistdir   : Directory of the position list files.
     snapfilebase : Name/path of the snapshot files (without the number). If you have multiple subsnapshots, put it
                         in the form ["snapshot_dir","/snapshot_name"] (both without the number).
     ICfile       : Name/path of the initial condition file of the parent simulation.
     cosmology    : Choice of cosmological parameters. Can be "planck" or "wmap".
     boxsize      : Size of the cosmological box, in comoving Mpc/h.
     zstart       : Initial redshift of the IC.
     seedsset     : Set of seeds. Can be "matteo" or "ludwig", but other options can easily be added in the code.
     parentres    : Resolution level of the parent simulation.
     lowestres    : Lowest resolution level for the zoom-in simulations based on this parent.

    Feel free to add options for more parent sims.
    """
    config_path = Path(__file__).parent / "../config.yaml"
    with open(config_path, "r") as ymlfile:
        cfg = yaml.load(ymlfile)
    for section in cfg["parent_setups"]:
        if parentlabel == section.key:
            return section["poslistdir"], section["snapfilebase"], section["ICfile"], section["cosmology"], \
                   section["boxsize"], section["zstart"], section["seedsset"], section["parentres"], section["lowestres"]
    raise IOError("parentlabel " + str(parentlabel) + " does not match any label in config.yaml")
