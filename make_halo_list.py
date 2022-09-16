import os
from src.rockstartools import rockstarselect
from src.read_config import getproperties

parentlabel = "parent3"
treefile = '/ptmp/mpa/mfrigo/musicsims/parent3/Rockstar/trees/tree_0_0_0.dat'
workdir, snapfilebase, parent_ic, cosmology, boxsize, zstart, seedsset, parentres, lowestres = getproperties(
    parentlabel)
if not os.path.exists(workdir):
    os.mkdir(workdir)
halolistfile = workdir + "/halotracing/halolist.txt"
if not os.path.exists(halolistfile):
    rockstarselect(treefile, mvirmin=10e11, mvirmax=10e12, mindist=10.)
else:
    print("Halo list file already exists in "+halolistfile+"; delete it if you want to recreate it")