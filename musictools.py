import numpy as np
from subprocess import call
from snap_format_adapter import shift_parts


def musicgadget3(haloid, parentlabel, highestres=10, initialpad=8, regionmode='ellipsoid', parentres=None,
                 lowestres=None, poslistdir=None):
    """
    Convenience function that creates MUSIC parameter file for a given halo, runs MUSIC and then merges the
    particle groups 3,4,5 in order to obtain an initial condition file to be used with Gadget3.

    REQUIRES Music executable to be in the current directory or in the path.

    Parameters:
     haloid      : Id of the halo. The position list file must be in the form poslist[haloid].txt .
     parentlabel : Properties of the parent simulation, taken from getproperties().
     highestres  : Maximum resolution level of the IC file.
     initialpad  : Size of the intermediate resolution levels. If too large, will use the maximum padding allowed.
     regionmode  : How the refinement region is defined. Can be "ellipsoid", "convex_hull" or "box".
     parentres   : Resolution level of the parent simulation. If None, picks the one from getproperties.
     lowestres   : Minimum resolution level of the IC file. If None, picks the one from getproperties.
     poslistdir  : Directory where the position list files are stored. If None, picks the one from getproperties.
    """
    import os.path
    poslistdir0, snapfilebase, ICfile, cosmology, boxsize, zstart, seedsset, parentres0, lowestres0 = getproperties(
        parentlabel)
    if poslistdir == None:
        poslistdir = poslistdir0
    if parentres == None:
        parentres = parentres0
    if lowestres == None:
        lowestres = lowestres0
    pad = initialpad
    foundmaxpad = False
    call("rm ICraw" + str(haloid) + "_" + str(highestres) + ".gdt*", shell=True)
    numfiles = 1
    if highestres >= 13:
        numfiles = 16
    while foundmaxpad == False:
        writemusicparam(haloid, highestres, parentres, lowestres, padding=pad, filename="param_zoominC.inp",
                        poslistdir=poslistdir, regionmode=regionmode, seedsset=seedsset, boxsize=boxsize,
                        cosmology=cosmology, zstart=zstart, numfiles=numfiles,
                        outname='ICraw' + str(haloid) + '_' + str(highestres) + '.gdt')
        call(["MUSIC", "param_zoominC.inp"])
        if os.path.isfile("ICraw" + str(haloid) + "_" + str(highestres) + ".gdt") or os.path.isfile(
                "ICraw" + str(haloid) + "_" + str(highestres) + ".gdt.0"):
            foundmaxpad = True
        else:
            pad = pad - 1
            if pad == 0:
                print("MUSIC fails to create the ICs with every padding value (perhaps the halo is too big for this boxsize?)")
                break
    shift_parts(infname='ICraw' + str(haloid) + '_' + str(highestres) + '.gdt',
                outfname="IC" + str(haloid) + "_" + str(highestres) + ".gdt")
    print("Zoom-in IC created: IC" + str(haloid) + "_" + str(highestres) + ".gdt")
    print("Resolution: " + str(highestres))
    print("Padding: " + str(pad))
#  call('rm ICraw'+str(haloid)+'_'+str(highestres)+'.gdt*', shell=True)


def writemusicparam(haloid, highestres=11, parentres=9, lowestres=4, regionmode='ellipsoid',
                    filename="param_zoominC.inp", padding=8, poslistdir='/ptmp/mpa/mfrigo/musicsims/parent3',
                    baryons='yes', seedsset='matteo', boxsize=72, cosmology='planck', zstart=43.029476,
                    numfiles=1, outname='./ICraw.gdt'):
    """
    Writes a music parameter file.

    Parameters:
     haloid     : ID of the halo to resimulate, as in the name of the poslist file (poslist<haloid>.txt).
     highestres : Highest resolution level (refinement region).
     parentres  : Resolution level of the parent simulation.
     lowestres  : Lowest resolution level.
     regionmode : How to define the refinement region. Can be 'ellipsoid', 'convex_hull', or 'box'.
     filename   : Name of the parameter file to generate.
     padding    : Padding (Number of cells for each intermediate resolution level).
     poslistdir : Directory of the position list files.
     baryons    : If to include baryons in the IC. Can be "yes" or "no".
     seedsset   : Set of seeds. Can be "matteo" or "ludwig", but other options can easily be added in the code.
     boxsize    : Size of the cosmological box, in comoving Mpc/h.
     cosmology  : Choice of cosmological parameters. Can be "planck" or "wmap".
     zstart     : Initial redshift of the IC.
     numfiles   : The output Gadget ICs will be subdivided in this number of subfiles.
     outname    : Name of the output file(s).

    Additional parameters are hard-coded into this function: initial redshift, cosmological parameters,...
    """
    if seedsset == 'matteo':
        seeds = [12345, 23336, 34567, 12146, 23111, 24343, 45454, 87654, 98765,
                 34987]  # seeds for levels 7,8,9,10,11,12,13,14,15,16
    if seedsset == 'ludwig':
        seeds = [12345, 23456, 34567, 90341, 56789, 67890, 45454, 87654, 98765,
                 34987]  # (LUDWIG sim) seeds for levels 7,8,9,10,11,12; I added the others
    if cosmology == 'planck':
        cospam = [0.307, 0.693, 0.0483, 67.77, 0.8288, 0.9611]  # Omega_m, Omega_L, Omega_b, H0, sigma_8, nspec
    if cosmology == 'wmap':
        cospam = [0.26, 0.74, 0.044, 72., 0.84, 0.93]  # Omega_m, Omega_L, Omega_b, H0, sigma_8, nspec

    print("Creating MUSIC parameter file")

    param = []
    param.append(r"[setup]" + "\n")
    param.append(r"boxlength		= " + str(boxsize) + "\n")
    param.append(r"zstart			= " + str(zstart) + "\n")
    param.append(r"levelmin		= " + str(lowestres) + "\n")
    param.append(r"levelmin_TF		= " + str(parentres) + "\n")
    param.append(r"levelmax		= " + str(highestres) + "\n")
    param.append(r"padding			= " + str(padding) + "\n")
    param.append(r"overlap			= 4  ##4" + "\n")
    if regionmode == 'box':
        poslist = np.loadtxt(str(poslistdir) + "/poslist" + str(haloid) + ".txt")
        boxcenter = [np.mean(poslist[:][0]), np.mean(poslist[:][1]), np.mean(poslist[:][2])]
        boxextent = [np.max(abs(poslist[:][0] - boxcenter[0])), np.max(abs(poslist[:][1] - boxcenter[1])),
                     np.max(abs(poslist[:][2] - boxcenter[2]))]
        param.append(
            r"ref_center		= " + str(boxcenter[0]) + ", " + str(boxcenter[1]) + ", " + str(boxcenter[2]) + "\n")
        param.append(
            r"ref_extent		= " + str(boxextent[0]) + ", " + str(boxextent[1]) + ", " + str(boxextent[2]) + "\n")
    param.append(r"align_top		= no" + "\n")
    param.append(r"baryons			= " + str(baryons) + "\n")
    param.append(r"use_2LPT		= yes" + "\n")
    param.append(r"use_LLA			= no" + "\n")
    param.append(r"periodic_TF		= yes" + "\n")
    if regionmode == 'ellipsoid':
        param.append(r"region			= ellipsoid" + "\n")
        param.append(r"region_point_file    = " + str(poslistdir) + "/poslist" + str(haloid) + ".txt" + "\n")
    if regionmode == 'convex_hull':
        param.append(r"region			= convex_hull" + "\n")
        param.append(r"region_point_file	= " + str(poslistdir) + "/poslist" + str(haloid) + ".txt" + "\n")
    param.append(r"" + "\n")
    param.append(r"[cosmology]" + "\n")
    param.append(r"Omega_m			= " + str(cospam[0]) + "\n")
    param.append(r"Omega_L			= " + str(cospam[1]) + "\n")
    param.append(r"Omega_b			= " + str(cospam[2]) + "\n")
    param.append(r"H0			= " + str(cospam[3]) + "\n")
    param.append(r"sigma_8			= " + str(cospam[4]) + "\n")
    param.append(r"nspec			= " + str(cospam[5]) + "\n")
    param.append(r"transfer		= eisenstein" + "\n")
    param.append(r"" + "\n")
    param.append(r"[random]" + "\n")
    for i in np.arange(17):
        if i > 6 and i <= highestres:
            if i >= parentres:
                param.append(r"seed[" + str(i) + "]     = " + str(seeds[i - 7]) + "\n")
            else:
                param.append(r"#seed[" + str(i) + "]     = " + str(seeds[i - 7]) + "\n")

    param.append(r"" + "\n")
    param.append(r"[output]" + "\n")
    param.append(r"##Gadget-2 (type=1: high-res particles, type=5: rest)" + "\n")
    param.append(r"format			= gadget2" + "\n")
    param.append(r"filename		= " + outname + "\n")
    param.append(r"gadget_lunit		= kpc" + "\n")
    param.append(r"gadget_num_files	= " + str(numfiles) + "\n")
    param.append(r"gadget_spreadcoarse	= yes" + "\n")
    param.append(r"#gadget_coarsetype	= 2" + "\n")
    param.append(r"" + "\n")
    param.append(r"[poisson]" + "\n")
    param.append(r"fft_fine		= yes" + "\n")
    param.append(r"accuracy		= 1e-5" + "\n")
    param.append(r"pre_smooth		= 3" + "\n")
    param.append(r"post_smooth		= 3" + "\n")
    param.append(r"smoother		= gs" + "\n")
    param.append(r"laplace_order		= 6" + "\n")
    param.append(r"grad_order		= 6" + "\n")
    of = open(filename, "w")
    for l in param:
        of.write(l)
    of.close()
