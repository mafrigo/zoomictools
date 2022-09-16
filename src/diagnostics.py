import numpy as np
import matplotlib.pyplot as plt
import pygad


def boxview(file, npixel=400):
    """
    Creates an image of the whole simulation box, color-coded by density.
    Requires pygadplus.py (ask Matteo Frigo).

    Parameters:
     file   : Gadget file of which you want to take a box picture.
     npixel : npixel x npixel is the number of pixels of the generated image.
    """
    s = pygad.Snap(file)
    try:
        mpcboxsize = 0.001 * s.boxsize.value
    except AttributeError:
        mpcboxsize = 0.001 * s.boxsize.item()
    print("Box size: " + str(mpcboxsize) + "cMpc h**-1")
    center = [1000 * mpcboxsize / 2., 1000 * mpcboxsize / 2., 1000 * mpcboxsize / 2.]
    s['pos'] -= center
    for i in np.arange(1):
        try:
            pygad.plotting.image(s, qty='mass', av=None, extent=str(mpcboxsize) + ' cMpc h_0**-1',
                                                Npx=npixel, showcbar=True)  # pygadplus
            plt.savefig("boxview" + str(i))
        except IOError:
            print("No particles of type: " + str(i))
    plt.show()
    plt.clf()


def zoomindiag(file, shrinkon='stars'):
    """
    Calculates mass statistics for zoom-in simulations: halo mass, stellar mass, gas mass,
    number of low resolution "intruder" particles, percentage of high resolution particles
    that end up in the halo.

    Parameters:
     file     : Name/path of the Gadget file.
     shrinkon : Type of particles used to determine the center of the halo. Can be "stars", "dm", "gas", "bh".
    """
    pygad.gadget.families['gr1'] = [1]
    pygad.gadget.families['gr2'] = [2]
    pygad.gadget.families['gr3'] = [3]

    def printstuff(s0):
        print("Dark halo mass: " + str("%e" % sum(s0.dm["mass"])))
        print("Luminous mass: " + str("%e" % (sum(s0.gas["mass"]) + sum(s0.stars["mass"]) + sum(s0.bh["mass"]))))
        print("Gas mass: " + str("%e" % sum(s0.gas["mass"])))
        print("Stellar mass: " + str("%e" % sum(s0.stars["mass"])))
        print("BH mass: " + str("%e" % sum(s0.bh["mass"])))
        print("Percentage of mass from low res particles: " + str(100 * (sum(s0.lowres["mass"])) / sum(s0["mass"])))
        print("Number of low res particles: " + str(len(s0.lowres["ID"])))
        print()

    try:
        s, halo, asf = pygad.prepare_zoom(file, shrink_on=shrinkon)
        Rvir, Mvir = pygad.analysis.virial_info(s)
        sb = s[pygad.BallMask(Rvir)]
        sc = s[pygad.BallMask(0.1 * Rvir)]
        Reff = halfmassradius(file)
        sd = s[pygad.BallMask(Reff)]
        print("\n \n")
        #    print "TOTAL: "
        #    printstuff(s)
        print("HALO: (Rvir=" + str(Rvir) + ")")
        printstuff(sb)
        print("GALAXY (within 0.1*Rvir=" + str(0.1 * Rvir) + "): ")
        printstuff(sc)
        print("INNER GALAXY (within r1/2=" + str(Reff) + " kpc): ")
        printstuff(sd)

        haloHighresDMMass = sum(sb.gr1["mass"])
        totalHighresDMMass = sum(s.gr1["mass"])
        print("Percentage of high res DM particles in the virial radius of the halo: ")
        print(str(100 * haloHighresDMMass / totalHighresDMMass) + " % ")

    except ValueError:

        print("No particles of type <" + shrinkon + "> available, or no overdensity to zoom into. Trying with dark matter. \n \n")
        s = pygad.Snap(file)
        s["mass"].convert_to("Msol")
        s1 = s.gr1
        center = pygad.analysis.center_of_mass(s1)
        s["pos"] = s["pos"] - center
        Rvir, Mvir = pygad.analysis.virial_info(s[pygad.BallMask('1 Mpc', fullsph=True)])
        sb = s[pygad.BallMask(Rvir, fullsph=True)]
        sc = s[pygad.BallMask(0.1 * Rvir, fullsph=True)]
        print()
        print()
        #    print "TOTAL: "
        #    printstuff(s)
        print("HALO: (Rvir=" + str(Rvir) + ")")
        printstuff(sb)
        print("GALAXY (within 0.1*Rvir=" + str(0.1 * Rvir) + "): ")
        printstuff(sc)
        haloHighresDMMass = sum(sb.gr1["mass"])
        totalHighresDMMass = sum(s.gr1["mass"])
        print("Percentage of high res DM particles in the virial radius of the halo: %.f " %(100 * haloHighresDMMass / totalHighresDMMass)) + "%"


def particlemasses(file):
    """
    Gives the minimum, average and maximum particle mass for each particle group in the specified file.
    """
    pygad.gadget.families['gr1'] = [1]
    pygad.gadget.families['gr2'] = [2]
    pygad.gadget.families['gr3'] = [3]
    s = pygad.Snap(file)
    print("\nParticles in each group: " + str(s.parts))
    s["mass"].convert_to("Msol")
    print("Particle masses:      (min)              (mean)               (max)")
    if s.parts[1] > 0:
        print("DM 1              " + str(np.min(s.gr1["mass"])) + " " + str(np.mean(s.gr1["mass"])) + " " + str(
            np.max(s.gr1["mass"])))
    else:
        print("No dm1")
    if s.parts[2] > 0:
        print("DM 2              " + str(np.min(s.gr2["mass"])) + " " + str(np.mean(s.gr2["mass"])) + " " + str(
            np.max(s.gr2["mass"])))
    else:
        print("No dm2")
    if s.parts[3] > 0:
        print("DM 3              " + str(np.min(s.gr3["mass"])) + " " + str(np.mean(s.gr3["mass"])) + " " + str(
            np.max(s.gr3["mass"])))
    else:
        print("No dm3")
    if s.parts[0] > 0:
        print("Gas              " + str(np.min(s.gas["mass"])) + " " + str(np.mean(s.gas["mass"])) + " " + str(
            np.max(s.gas["mass"])))
    else:
        print("No gas")
    if s.parts[4] > 0:
        print("Stars            " + str(np.min(s.stars["mass"])) + " " + str(np.mean(s.stars["mass"])) + " " + str(
            np.max(s.stars["mass"])))
    else:
        print("No stars")
    if s.parts[5] > 0:
        print("BH               " + str(np.min(s.bh["mass"])) + " " + str(np.mean(s.bh["mass"])) + " " + str(
            np.max(s.bh["mass"])))
    else:
        print("No black holes")


def halfmassradius(file, onlystars=True):
    """
    Calculates half mass radius of the galaxy in a zoom-in simulation.

    Parameters:
     file      : Name/path of the Gadget file.
     onlystars : If True, only stars are considered in computing the half-mass radius
                   (so it's more like a 3D half-light radius).
    """
    s, halo, asf = pygad.prepare_zoom(file)
    Rvir, Mvir = pygad.analysis.virial_info(s)
    sc = s[pygad.BallMask(0.1 * Rvir)]
    if onlystars:
        sc = sc.stars
    stellarm = sum(sc["mass"])
    sortindex = np.argsort(sc["r"])
    radii = np.array(sc["r"])
    radii = radii[sortindex]
    masses = np.array(sc["mass"])
    masses = masses[sortindex]
    cumumass = np.cumsum(masses)
    Reff1 = (radii[cumumass < 0.5 * stellarm])[-1]
    Reff2 = (radii[cumumass > 0.5 * stellarm])[0]
    m1 = (cumumass[cumumass < 0.5 * stellarm])[-1] - 0.5 * stellarm
    m2 = (cumumass[cumumass > 0.5 * stellarm])[0] - 0.5 * stellarm
    Reff = Reff1 + (abs(m1) / (m2 - m1)) * (Reff2 - Reff1)
    print("Half mass radius: " + str(Reff))
    return Reff