import os
from src.main import createIC, multipleICcreator
from src.read_config import getproperties
from src.diagnostics import zoomindiag

parentlabel = "parent3"
haloid = 40
zoom_highest_res = 12
workdir, snapfilebase, parent_ic, cosmology, boxsize, zstart, seedsset, parentres, lowestres = getproperties(
    parentlabel)
if not os.path.exists(workdir):
    os.mkdir(workdir)
if not os.path.exists(workdir + "/halotracing/halolist.txt"):
    raise IOError("No halolist file found - run make_halo_list.py first")
createIC(haloid, parentlabel, halolistfolder=workdir + '/halotracing', zoominlevel=zoom_highest_res)
# multipleICcreator(parentlabel, halolistfolder=workdir+'/halotracing', zoominlevel=12, startingatnumber=40)