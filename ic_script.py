from src.main import createIC, multipleICcreator
from src.rockstartools import rockstarselect

#rockstarselect('/ptmp/mpa/mfrigo/musicsims/parent3/Rockstar/trees/tree_0_0_0.dat', mvirmin=10e11, mvirmax=10e12, mindist=10.)
createIC(40, 'parent3', halolistfolder='/ptmp/mpa/mfrigo/musicsims/parent3/halotracing', zoominlevel=12)
#multipleICcreator('parent3', halolistfolder='./parent3', zoominlevel=12, startingatnumber=40)

