import os
from src.main import createIC, multipleICcreator
from src.rockstartools import rockstarselect
from src.read_config import getproperties

parentlabel = "parent3"
treefile = '/ptmp/mpa/mfrigo/musicsims/parent3/Rockstar/trees/tree_0_0_0.dat'
workdir, snapfilebase, parent_ic, cosmology, boxsize, zstart, seedsset, parentres, lowestres = getproperties(
    parentlabel)
if not os.path.exists(workdir):
    os.mkdir(workdir)
if not os.path.exists(workdir + "/halotracing/halolist.txt"):
    rockstarselect(treefile, mvirmin=10e11, mvirmax=10e12, mindist=10.)
createIC(40, parentlabel, halolistfolder=workdir + '/halotracing', zoominlevel=12)
# multipleICcreator(parentlabel, halolistfolder=workdir+'/halotracing', zoominlevel=12, startingatnumber=40)
