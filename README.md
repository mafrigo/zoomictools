# Required tools:
- Music (https://www-n.oca.eu/ohahn/MUSIC/)
- Gadget2/3
- Rockstar (https://bitbucket.org/gfcstanford/rockstar)
- Python modules in requirements.txt

# Workflow:
1) Create IC for dark matter only parent simulation with MUSIC.
2) Run parent simulation with Gadget.
3) Run rockstar on parent simulation to create a halo tree file.
4) Run rockstarselect() with the tree file to select halos you are interested in, 
based on min/max virial mass and distance to other halos.
5) Create zoom initial condition with createIC(). Use id of a galaxy from the rockstarselect() output (they are different from the rockstar ids).
Or use multipleICcreator() to create many ICs at the same time, starting from a given galaxy ID and going to higher numbers.
6) Run the zoom simulation with Gadget, including gas, star formation, etc.