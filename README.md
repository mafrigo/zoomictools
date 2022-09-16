# Required tools:
- Music (https://www-n.oca.eu/ohahn/MUSIC/)
- Gadget2/3
- Rockstar (https://bitbucket.org/gfcstanford/rockstar)
- Python >3.4
- Python modules in requirements.txt (notably pygad: https://bitbucket.org/broett/pygad/)

# Workflow:
1) Create IC for dark matter only parent simulation with MUSIC.
2) Run parent simulation with Gadget.
3) Run rockstar on parent simulation to create a halo tree file.
4) Add configuration for your IC setup in config.yaml. music_exec should link to your music executable.
In particular you should use the same cosmology, zstart, boxsize and seeds (in the resolutions levels where parent and zoom overlap) as in the parent.
5) Run rockstarselect() with the tree file to select halos you are interested in, 
based on min/max virial mass and distance to other halos.
6) Create zoom initial condition with createIC(). Use id of a galaxy from the rockstarselect() output (they are different from the rockstar ids)
, your resolution level of choice, and the label of your new config.yaml setup. You can also use multipleICcreator() to create many ICs at the same time, 
starting from a given galaxy ID and going to higher numbers.
7) Run the zoom simulation with Gadget, including gas, star formation, etc.
8) Check whether virial mass of zoomed galaxy matches value in parent, and whether there are intruder low resolution particles (with diagnostics.py)
Note: ic_script.py does steps 5 and 6