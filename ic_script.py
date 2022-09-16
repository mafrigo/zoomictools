from src.main import createIC, multipleICcreator
from src.rockstartools import rockstarselect
from src.read_config import getproperties

parentlabel="parent3"
#rockstarselect('/ptmp/mpa/mfrigo/musicsims/parent3/Rockstar/trees/tree_0_0_0.dat', mvirmin=10e11, mvirmax=10e12, mindist=10.)
poslistdir, snapfilebase, ICfile, cosmology, boxsize, zstart, seedsset, parentres, lowestres = getproperties(parentlabel)
createIC(40, parentlabel, halolistfolder=poslistdir+'/halotracing', zoominlevel=12)
#multipleICcreator(parentlabel, halolistfolder='./parent3', zoominlevel=12, startingatnumber=40)

