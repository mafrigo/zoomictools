import numpy as np
import os
from subprocess import call
import pygad

from .rockstartools import readRockstarHaloTree
from .read_config import getproperties

def trace(haloid, parentlabel, eps=10., eps2=300., tracefactor=2., preemptivecut=True, nmax=94, npp=0,
          haloloctype='zit', halolocfile=None, snaps_to_be_traced="all", forcegc=False, outputdir='.'):
    """
    Traces a halo in an uniform dark matter only simulation from a given snapshot to the initial conditions,
               and saves the positions of the particles which end up forming the halo.

    Parameters:
      haloid           : ID of the halo you want to trace. To be selected with rockstarselect().
      parentlabel      : Label of the parent simulation, to be used in getproperties().
      eps              : Adjustment to the size of the sphere inside which particle are traced.
      eps2             : Size of the pseudoparticle cloud applied at the end of the tracing to each selected particle.
      tracefactor      : Halo particles are traced within tracefactor*r200 (+eps)
      preemptivecut    : If True, cuts the snapshot at 50*R200 to shorten the process. Used only if haloloctype='list'.
      nmax             : Snapshot at which the tracing starts.
      npp              : Number of pseudopartiles generated at the end of the tracing for each original particle. If 0,
                          no pseudoparticle expansion is performed.
      haloloctype      : How the halo is located; can be 'list', 'tree', or 'zit'.
                         'list': Gets the position of the halo at redshift 0 from the Rockstar halolist file, then tracks
                                  the halo back in time; more computationally expensive.
                         'tree': Gets the position of the halo at every snapshot from the Rockstar merger tree file.
                         'zit' : Gets the position of the halo at every snapshot from the rockstarselect() output (with
                                  tracehalos = True). Computationally the fastest option.
      halolocfile      : Name/path of the halo location file, of the kind required for the selected haloloctype. If None,
                          a guess is made according to haloloctype (e.g. 'halos_0.0.ascii' if haloloctype='list').
      snaps_to_be_traced : List of snapshots in which the halo is traced. Must be in growing order (ex: [30,60,94]). If "all",
                          all available snapshots are used.
      forcegc          : If True, forces the python garbage collection after each analyzed snapshot. Useful only when
                          memory problems show up (huge snapshot files).
    """
    poslistdir0, filebase, ICfile, cosmology0, boxsize0, zstart0, seedsset0, parentres0, lowestres0 = getproperties(
        parentlabel)
    nmin = 0  # indices of the snapshots to trace (from nmax to nmin)
    if snaps_to_be_traced == 'all':
        snaps_to_be_traced = range(nmin, nmax + 1)  # all snaps

    if haloloctype == 'tree' or haloloctype == 'zit':
        print('Tracing halo number ' + str(haloid))
        if haloloctype == 'tree':
            if halolocfile is None:
                halolocfile = './Rockstar/trees/tree_0_0_0.txt'
            scaleRS, m200RS, r200RS, halocenterRS = readRockstarHaloTree(halolocfile, haloid)
        if haloloctype == 'zit':
            if halolocfile is None:
                halolocfile = './halotracing/halo' + str(haloid) + '.txt'
            scaleRS, m200RS, r200RS, xRS, yRS, zRS = np.transpose(np.loadtxt(halolocfile))
            halocenterRS = [xRS, yRS, zRS]
        nmin = nmax - len(r200RS) + 1
        snaps_to_be_traced = [i for i in snaps_to_be_traced if i >= nmin]
    if haloloctype == 'list':
        if halolocfile is None:
            halolocfile = 'halos_0.0.ascii'
        data = np.loadtxt(halolocfile)  # Opening rockstar halolist to get halo properties
        m200 = data[haloid, 2] / 10 ** 12
        r200 = data[haloid, 4]
        xinit = data[haloid, 8] * 1000.  # in kpc
        yinit = data[haloid, 9] * 1000.
        zinit = data[haloid, 10] * 1000.
        halocenter = [xinit, yinit, zinit]
        print('Tracing halo number ' + str(haloid) + ', with initial center:')
        print(xinit, yinit, zinit)
        print('initial virial radius: ' + str(r200))
        print('and virial mass: ' + str(m200) + ' *10^12 solar masses')

    idlist = list([])
    snaps_to_be_traced = list(reversed(snaps_to_be_traced))  # reversing snapshot list
    for i in snaps_to_be_traced:  # main loop over snapshots
        print('\n Tracing halo in snapshot ' + str(i))
        if np.shape(filebase) == ():
            filename = filebase + '0%02i' % i
        else:
            filename = ""
            for k in np.arange(len(filebase)):
                filename = filename + filebase[k] + '0%02i' % i
        print("Loaded snapshot at " + filename)
        if haloloctype == 'tree' or haloloctype == 'zit':
            j = nmax - i
            halocenter = np.array(halocenterRS)[:, j]
            r200 = r200RS[j]
            m200 = m200RS[j]
            if os.path.exists(filename + '.0') and os.path.exists(filename + '.1'):
                i = 0
                nextfileexists = True
                idlist0 = []
                while nextfileexists:
                    print("Opening subfile number " + str(i))
                    s = _partialSnap(filename + '.' + str(i))
                    if np.min(np.sqrt(np.sum((s["pos"] - halocenter) ** 2, axis=1))) < (1.5 * tracefactor * r200 + eps):
                        ball = s[pygad.BallMask(R=str(tracefactor * r200 + eps) + ' ckpc h_0**-1',
                                                center=halocenter)]  # Cutting snapshot at 2 virial radii from halo center
                        idlist0 = idlist0 + list(ball["ID"])  # getting the halo IDs in this snapshot
                        print("Updated number of halo particles: " + str(len(idlist0)))
                    del s
                    i = i + 1
                    nextfileexists = os.path.exists(filename + '.' + str(i))
            else:
                s = pygad.Snap(filename) # Loading snapshot
                ball = s[pygad.BallMask(R=str(tracefactor * r200 + eps) + ' ckpc h_0**-1',
                                        center=halocenter)]  #Cutting snapshot at 2 virial radii from halo center
                idlist0 = list(ball["ID"])  # getting the halo IDs in this snapshot
                del s
        else:
            s = pygad.Snap(filename)
            if r200 < 50.:
                r200 = 50.  # just in case r200 goes to 0
            if preemptivecut:
                s = s[pygad.BallMask(R=str(np.max([50. * r200, 20000])) + ' ckpc h_0**-1',
                                     center=halocenter)]  # cutting snapshot at 20 Mpc for faster computation
            if len(idlist) > 0:  # Calculating new halo center
                halocenter = pygad.analysis.center_of_mass(s[pygad.IDMask(idlist)])  # oldball)
                halocenter = pygad.analysis.shrinking_sphere(s, halocenter, R=tracefactor * r200, periodic=True,
                                                             shrink_factor=0.93, stop_N=13)
            r200, m200 = pygad.analysis.virial_info(s, center=halocenter)  # Calculating virial radius
            try:
                r200 = r200.value()
            except AttributeError:
                r200 = r200.item()
            ball = s[pygad.BallMask(R=str(tracefactor * r200 + eps) + ' ckpc h_0**-1',
                                    center=halocenter)]  # ' kpc',center=halocenter,fullsph=True)]   #Cutting snapshot at 2 virial radii from halo center
            idlist0 = list(ball["ID"])  # getting the halo IDs in this snapshot
            del s
        print('Ball center: ' + str(halocenter))
        print('r200: ' + str(r200) + ', m200: ' + str(m200))
        print(str(len(idlist0)) + ' particles within ' + str(tracefactor) + '*R200')
        idlist = list(set(idlist + idlist0))
        print("Total number of IDs in refinement region: " + str(len(idlist)))
        if forcegc:
            pygad.gc_full_collect()

    if os.path.exists(ICfile):
        print("Opening IC file:")
        call(["rm", "idlist" + str(haloid) + ".txt"])
        call(["rm", "poslist" + str(haloid) + ".txt"])
        s = pygad.Snap(ICfile)
        try:
            boxsizempc = 0.001 * s.boxsize.value
        except AttributeError:
            boxsizempc = 0.001 * s.boxsize.item()
        highres = s[pygad.IDMask(idlist)]
        np.savetxt('idlist' + str(haloid) + '.txt', highres["ID"])  # list of the IDs
        if eps2 is not 0:  # pseudoparticle expansion to increase the size of the refinement region
            lenpseudosnap = len(highres["pos"][:, 0]) * npp
            xg = np.ravel(highres["pos"][:, 0].repeat(npp) + np.clip(
                np.random.normal(np.zeros(lenpseudosnap), eps2, lenpseudosnap), 0., 1.5 * eps2))
            yg = np.ravel(highres["pos"][:, 1].repeat(npp) + np.clip(
                np.random.normal(np.zeros(lenpseudosnap), eps2, lenpseudosnap), 0., 1.5 * eps2))
            zg = np.ravel(highres["pos"][:, 2].repeat(npp) + np.clip(
                np.random.normal(np.zeros(lenpseudosnap), eps2, lenpseudosnap), 0., 1.5 * eps2))
            poslist = np.transpose(np.array([xg, yg, zg]) / (boxsizempc * 1000.))
        else:
            poslist = highres["pos"] / (boxsizempc * 1000.)
    else:
        if os.path.exists(ICfile + '.0'):
            print("Opening IC subfiles separately:")
            call(["rm", "idlist" + str(haloid) + ".txt"])
            call(["rm", "poslist" + str(haloid) + ".txt"])
            np.savetxt('idlist' + str(haloid) + '.txt', idlist)  # list of the IDs
            poslist = np.empty((0, 3))
            i = 0
            nextfileexists = True
            while nextfileexists:
                s = _partialSnap(ICfile + '.' + str(i))
                print(ICfile + '.' + str(i))
                try:
                    boxsizempc = 0.001 * s.boxsize.value
                except AttributeError:
                    boxsizempc = 0.001 * s.boxsize.item()
                highres = s[pygad.IDMask(idlist)]
                if npp > 0:  # pseudoparticle expansion to increase the size of the refinement region
                    lenpseudosnap = len(highres["pos"][:, 0]) * npp
                    xg = np.ravel(highres["pos"][:, 0].repeat(npp) + np.clip(
                        np.random.normal(np.zeros(lenpseudosnap), eps2, lenpseudosnap), 0., 1.5 * eps2))
                    yg = np.ravel(highres["pos"][:, 1].repeat(npp) + np.clip(
                        np.random.normal(np.zeros(lenpseudosnap), eps2, lenpseudosnap), 0., 1.5 * eps2))
                    zg = np.ravel(highres["pos"][:, 2].repeat(npp) + np.clip(
                        np.random.normal(np.zeros(lenpseudosnap), eps2, lenpseudosnap), 0., 1.5 * eps2))
                    poslist0 = list(np.transpose(np.array([xg, yg, zg]) / (boxsizempc * 1000.)))
                else:
                    poslist0 = highres["pos"] / (boxsizempc * 1000.)
                if i == 0:
                    poslist = poslist0
                else:
                    poslist = poslist + poslist0
                i = i + 1
                nextfileexists = os.path.exists(ICfile + '.' + str(i))
        else:
            raise IOError('IC file not found')

    print('Printing final list ')
    np.savetxt('poslist' + str(haloid) + '.txt', poslist)  # list of the positions normalized to the box size
    print('Halo number: ' + str(haloid))
    print('Number of unique IDs: ' + str(len(idlist)))
    print('Number of poslist entries: ' + str(len(poslist)))


def _partialSnap(filename, physical=False, load_double_prec=False, cosmological=None,
                 gad_units=None, unclear_blocks=None):
    '''
    Create a snapshot from a single subfile, rather than trying to open
     all subfiles at once like the default _Snap function in pygad.

    Args:
        filename (str):         The path to the snapshot. If it is distributed
                                over several files, you shall omit the trailing
                                (of inbetween in case of an HDF5 file) '.0'.
        physical (bool):        Whether to convert to physical units on loading.
        load_double_prec (bool):Force to load all blocks in double precision.
                                Equivalent with setting the snapshots attribute.
        cosmological (bool):    Explicitly tell if the simulation is a
                                cosmological one.
        gad_units (dict):       Alternative base units (LENGTH, VELOCITY, MASS)
                                for this snapshot. The default base units units
                                are updated, meaning one can also just change one
                                them.
        unclear_blocks (str):   What to do the blocks for which the block info is
                                unclear (cannot be infered). Possible modes are:
                                * exception:    raise an IOError
                                * warning:      print a warning to the stderr
                                * ignore:       guess what
                                If it is None, the value from the `gadget.cfg` is
                                taken.

    Raises:
        IOError:            If the snapshot does not exist.
        RuntimeError:       If the information could not be infered or if the
                            given dtype of the given family is unknown.
    '''
    import pygad.snapshot as snapshot
    from pygad.snapshot.sim_arr import SimArr
    import pygad.gadget as gadget
    import pygad.physics as physics
    import pygad.environment as environment

    #    import pygad.utils as utils
    filename = os.path.expandvars(filename)
    filename = os.path.expanduser(filename)
    # handle different filenames, e.g. cases where the snapshot is distributed
    # over severale files
    base, suffix = os.path.splitext(filename)
    if suffix != '.hdf5':
        base, suffix = filename, ''
    if not os.path.exists(filename):
        filename = base + '.0' + suffix
        if not os.path.exists(filename):
            raise IOError('Snapshot "%s%s" does not exist!' % (base, suffix))

    s = snapshot.snapshot._Snap(gad_units=gad_units, physical=physical, cosmological=cosmological)
    s._filename = os.path.abspath(base + suffix)
    s._descriptor = os.path.basename(base) + suffix
    # need first (maybe only) file for basic information
    greader = gadget.FileReader(filename, unclear_blocks=unclear_blocks)
    s._file_handlers = [greader]

    s._block_avail = {
        block.name: block.ptypes
        for block in greader.infos()
        if block.dtype is not None}

    s._N_part = map(int, greader.header['N_part'])  # _all'])
    s._time = greader.header['time']
    s._redshift = greader.header['redshift']
    s._boxsize = SimArr(greader.header['boxsize'],
                        units=s._gad_units['LENGTH'],
                        snap=s)
    s._cosmology = physics.FLRWCosmo(  # Note: Omega_b at default!
        h_0=greader.header['h_0'],
        Omega_Lambda=greader.header['Omega_Lambda'],
        Omega_m=greader.header['Omega_m'])
    s._properties = {k: v for k, v in greader.header.iteritems()
                     if k not in ['N_part', 'mass', 'time', 'redshift', 'N_part_all',
                                  'N_files', 'h_0', 'Omega_Lambda', 'Omega_m', 'boxsize',
                                  'unused']}
    if s._cosmological is None:
        s._cosmological = abs(s.scale_factor - s.time) < 1e-6
    s._load_double_prec = bool(load_double_prec)

    if physical:
        s._boxsize.convert_to(s._boxsize.units.free_of_factors(['a', 'h_0']),
                              subs=s)

    # Process block names: make standard names lower case (except ID) and replace
    # spaces with underscores for HDF5 names. Also strip names
    s._load_name = {}
    for name, block in s._block_avail.items():
        if s._file_handlers[0]._format == 3 \
                and '%-4s' % name not in gadget.std_name_to_HDF5:
            new_name = name.strip()
        else:
            new_name = name.strip().lower()

        # some renaming
        if new_name in ['id', 'z']:
            new_name = new_name.upper()
        elif new_name == 'age':
            new_name = 'form_time'

        s._load_name[new_name] = name
        s._block_avail[new_name] = s._block_avail[name]
        if name != new_name:
            #           if environment.verbose >= environment.VERBOSE_TALKY \
            #                and new_name.lower() != name.strip().lower():
            #            print 'renamed block "%s" to %s' % (name, new_name)
            del s._block_avail[name]  # blocks should not appear twice
    # now the mass block is named 'mass' for all cases (HDF5 or other)
    s._block_avail['mass'] = [n > 0 for n in s._N_part]

    s.fill_derived_rules()

    s._descriptor = '"' + s._descriptor + '"'

    s._N_files = greader.header['N_files']

    return s