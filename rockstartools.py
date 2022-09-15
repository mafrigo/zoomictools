import numpy as np
import os
from subprocess import call


def readRockstarHaloTree(file, haloid):
    """
    Reads the Rockstar merger tree file and finds the properties of the selected halo (haloid) across the simulation. Used in trace().
    """
    data = np.loadtxt(file, skiprows=45)
    print('Rockstar tree data loaded')
    scale = data[:, 0]
    id = data[:, 1]
    desc_id = data[:, 3]
    mvir = data[:, 10]
    rvir = data[:, 11]
    x = data[:, 17]
    y = data[:, 18]
    z = data[:, 19]
    return _buildhaloloclist(haloid, scale, id, desc_id, mvir, rvir, x, y, z)


def _buildhaloloclist(haloid, scale, id, desc_id, mvir, rvir, x, y, z):
    mainrvir = []
    mainmvir = []
    mainx = []
    mainy = []
    mainz = []
    scalelist = list(reversed(np.sort(list(set(scale)))))
    prev_id = -1  # for final remnant
    for sc in scalelist:
        if prev_id == -1:
            most_massive_prog = haloid
            if desc_id[id == haloid] == -1.:
                print("Good choice.")
            else:
                print("Well, halo " + str(haloid) + " is not a final remnant, you should reconsider your choice.")
        else:
            progs = id[desc_id == prev_id]
            if len(progs) != 0:
                most_massive_prog = id[np.logical_and(desc_id == prev_id, mvir == np.max(mvir[desc_id == prev_id]))]
            else:
                continue
        mainrvir.append(rvir[id == most_massive_prog][0])
        mainmvir.append(mvir[id == most_massive_prog][0])
        mainx.append(x[id == most_massive_prog][0] * 1000.)  # in kpc
        mainy.append(y[id == most_massive_prog][0] * 1000.)
        mainz.append(z[id == most_massive_prog][0] * 1000.)
        prev_id = most_massive_prog
    return scalelist[0:len(mainmvir)], mainmvir, mainrvir, [mainx, mainy, mainz]


def rockstarselect(mvirmin, mvirmax, mindist, file='tree_0_0_0.dat', type='tree', savelist=True,
                   tracehalos=True, tracefolder='./halotracing/', listfile='halolist.txt'):
    """
    Reads the Rockstar merger tree or halo list file and finds halos with mvirmin < halo mass < mvirmax
    and whose closest halo is at least mindist away. Prints the selected halos on screen.

    Parameters:
     mvirmin     : Minimum mass of the halo.
     mvirmax     : Maximum mass of the halo.
     mindist     : Minimum distance to the next significant (M>10**11 Msol) halo.
     file        : (Optional) Path of the Rockstar merger tree / halo list file.
     type        : Type of Rockstar file. Can be 'tree' (merger trees) or 'list' (halo list).
     savelist    : If true, saves the list of selected halos in tracefolder/listfile.
     tracehalos  : If true, traces each selected halo across the simulation and saves its
                   info in the "tracefolder" folder. Only possible if type='tree'.
     tracefolder : Folder where the halo list and the halo tracing files are saved.
     listfile    : Name of the halo list file, saved in the "tracefolder" folder.
    """
    print("Opening Rockstar file...")
    if type == 'tree':
        try:
            data = np.loadtxt(file, skiprows=45)
        except ValueError:
            print("Trying with genfromtxt...")
            data = np.genfromtxt(file, skip_header=45)
        print("Merger tree file loaded.")
        desc_id = data[:, 3]
        scale = data[:, 0][desc_id == -1.]
        id = data[:, 1][desc_id == -1.]
        mvir = data[:, 10][desc_id == -1.]
        rvir = data[:, 11][desc_id == -1.]
        x = data[:, 17][desc_id == -1.]
        y = data[:, 18][desc_id == -1.]
        z = data[:, 19][desc_id == -1.]
    if type == 'list':
        data = np.loadtxt(file)
        print("Halo list file loaded.")
        id = data[:, 0]
        mvir = data[:, 2]
        rvir = data[:, 4]
        x = data[:, 8]
        y = data[:, 9]
        z = data[:, 10]
    if savelist:
        print("Opening output file: " + tracefolder + listfile)
        if os.path.exists(tracefolder + listfile):
            lasthalo = np.loadtxt(tracefolder + listfile)[:, 0][-1]
            f = open(tracefolder + listfile, 'a')
        else:
            call(["mkdir", tracefolder])
            f = open(tracefolder + listfile, 'a')
            f.write(
                "#NewID   RockstarID    M200(10^12Msol)     R200(kpc)   Closest halo ID     Closest halo distance(Mpc) \n")
            lasthalo = -1
    goodids = id[np.logical_and(mvir >= mvirmin, mvir <= mvirmax)]
    print("NewID   RockstarID    M200(10^12Msol)     R200(kpc)   Closest halo ID     Closest halo distance(Mpc)")
    newid = lasthalo + 1
    for i in goodids:
        if mindist != 0.:
            closeid, closedist = _closesthalo(i, id, x, y, z, mvir)
        else:
            closeid = 0.
            closedist = 10000.
        if closedist >= mindist:
            print(newid, i, mvir[id == i][0] / (10 ** 12), rvir[id == i][0], closeid, closedist)
            if savelist:
                f.write(str(newid) + " " + str(i) + " " + str(mvir[id == i][0] / (10 ** 12)) + " " + str(
                    rvir[id == i][0]) + " " + str(closeid[0]) + " " + str(closedist) + '\n')
                if tracehalos:
                    treemask = (data[:, 29] == i)  # mask for selecting the data of this tree only
                    scaleT, mvirT, rvirT, centerT = _buildhaloloclist(i, scale=data[:, 0][treemask],
                                                                      id=data[:, 1][treemask],
                                                                      desc_id=data[:, 3][treemask],
                                                                      mvir=data[:, 10][treemask],
                                                                      rvir=data[:, 11][treemask],
                                                                      x=data[:, 17][treemask], y=data[:, 18][treemask],
                                                                      z=data[:, 19][treemask])
                    outtable = np.transpose(np.array(
                        [scaleT, mvirT, rvirT, np.array(centerT)[0, :], np.array(centerT)[1, :],
                         np.array(centerT)[2, :]]))
                    np.savetxt(tracefolder + 'halo' + str(int(newid)) + '.txt', outtable)
            newid = newid + 1
    if savelist:
        f.close()


def _closesthalo(i, id, x, y, z, mvir, minhalomasstobeconsidered=10 ** 11):
    """
    Finds the closest halo to the halo with id=i and its distance.

    Parameters:
     i     : id of the main halo
     id    : list of halo ids
     x,y,z : list of halo positions
     mvir  : list of halo masses
     minhalomasstobeconsidered : Minimum halo mass for which we check the distance.
    """
    closedist = 10 ** 30
    closeid = 0.
    r = np.sqrt((x - x[id == i]) ** 2 + (y - y[id == i]) ** 2 + (z - z[id == i]) ** 2)
    try:
        closedist = np.min(r[np.logical_and(r != 0, mvir >= minhalomasstobeconsidered)])
        closeid = id[r == closedist]
    except ValueError:
        print("No halo with mass larger than " + str(minhalomasstobeconsidered) + " !")
    return closeid, closedist