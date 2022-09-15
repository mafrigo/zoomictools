import numpy as np
import os


def shift_parts(infname='ICraw.gdt', outfname='IC.gdt', toGroup=[3, 4, 5], icform=1):
    """
    Merges particle groups in a Gadget snapshot. Useful to make the MUSIC output file compatible with Gadget3.

    Parameters:
     infname  : Name/path of the input file.
     outfname : Name/path of the output file.
     toGroup  : List of groups to be merged. All particles will be shifted to the group with the lower number in the list.
     icform   : Format of the Gadget snapshot (1, 2 or 3).
    """
    if not (np.sum(np.diff(toGroup) == 1) == (len(toGroup) - 1)):
        raise ValueError("The groups to be merged are not continous and monotonically increasing!")

    if np.sum(np.array(toGroup) > 5) + np.sum(np.array(toGroup) < 0):
        raise ValueError("The groups to be merged are not in the interval 0-5!")

    if not os.path.exists(infname):
        if not os.path.exists(infname + '.0'):
            raise IOError("File " + infname + " does not exist")
        else:
            header = _read_header(infname + '.0', icform=icform)
            Nfiles = header["num_files"]
            print("Number of subfiles: " + str(Nfiles))
            infnamelist = [infname + '.' + str(i) for i in np.arange(Nfiles)]
            outfnamelist = [outfname + '.' + str(i) for i in np.arange(Nfiles)]
    else:
        Nfiles = 1
        infnamelist = [infname]
        outfnamelist = [outfname]

    for j in np.arange(Nfiles):
        infname = infnamelist[j]
        outfname = outfnamelist[j]

        with open(outfname, "w") as outf, open(infname, "r") as inf:

            print("Merging groups %s in file \"%s\", writing to \"%s\"" % (str(toGroup), infname, outfname))
            header = _read_header(infname, icform=icform)
            print(header)
            _print_header_simple(header, infname)

            inf.seek(0, 2)
            inflen = inf.tell()
            inf.seek(0, 0)

            # read header block from infile to modify
            newHeader = _read_block(inf, headerType, icform=icform)

            floatType = np.dtype("f8") if newHeader['flag_doubleprecision'] else np.dtype("f4")
            idType = np.dtype("u4")

            # check whether we already need/have a mass block for a given type
            withMassBlock = np.array([newHeader["npart%i" % i] and not newHeader["mass%i" % i] for i in range(6)],
                                     dtype="bool")
            print(withMassBlock)
            oldNums = np.array([newHeader["npart%i" % i] for i in range(6)]).flatten()
            print(oldNums)

            masses = []
            for i in toGroup[1:]:
                # change header
                newHeader["npart%i" % toGroup[0]] += newHeader["npart%i" % i]
                newHeader["npart%i" % i] = 0
                newHeader["npartTotal%i" % toGroup[0]] += newHeader["npartTotal%i" % i]
                newHeader["npartTotal%i" % i] = 0
                newHeader['npartTotalHighWord%d' % toGroup[0]] += newHeader['npartTotalHighWord%d' % i]
                newHeader['npartTotalHighWord%d' % i] = 0

            for i in toGroup:
                masses.append(float(newHeader["mass%i" % i]))
                newHeader["mass%i" % i] = 0

            _print_header_simple(newHeader, outfname)

            allMassesEqual = not np.sum(np.abs(np.diff(masses)))

            # write changed header
            _write_block(outf, newHeader, icform=icform)

            # read/write positions
            _write_block(outf, _read_block(inf, floatType), icform=icform)
            print(inf.tell())
            # read/write velocities
            _write_block(outf, _read_block(inf, floatType), icform=icform)
            print(inf.tell())
            # read/write IDs
            _write_block(outf, _read_block(inf, idType), icform=icform)
            print(inf.tell())
            # read masses if there are any
            if np.sum(withMassBlock):
                oldMasses = _read_block(inf, floatType, icform=icform)
            print(inf.tell())
            # create new sub-blocks if needed
            if not allMassesEqual:
                if np.sum(withMassBlock):
                    oldMassSubBlocks = np.split(oldMasses, np.cumsum(oldNums[withMassBlock])[:-1])
                newMassSubBlocks = [np.array([m], floatType).repeat(oldNums[i]) for m, i in zip(masses, toGroup)]
                newMasses = []
                o = 0
                n = 0
                for i in range(6):
                    if withMassBlock[i]:
                        newMasses.append(oldMassSubBlocks[o])
                        o += 1
                    else:
                        if i in toGroup:
                            newMasses.append(newMassSubBlocks[n])
                    if i in toGroup:
                        n += 1

                newMasses = np.hstack(newMasses)
                print(newMasses)
            else:
                print("All masses are equal, not gonna write new mass sub-blocks.")
                newMasses = oldMasses

            # if there is a reason to write a mass block, do so
            if np.sum(withMassBlock) or not allMassesEqual:
                print(newMasses)
                if len(newMasses) > 0 or newHeader["npart3"] != 0:
                    _write_block(outf, newMasses, icform=icform)
                else:
                    print("There might be a problem")

            # read/write all other blocks until finished
            if newHeader["npart0"] != 0:
                while not inf.tell() >= inflen:
                    print("Copied additional block.")
                    additionalblocks = _read_block(inf, np.dtype("S1"), icform=icform)
                    print(additionalblocks)
                    if len(additionalblocks) > 0:
                        _write_block(outf, additionalblocks, icform=icform)
    print("IC file created!")


headerType = np.dtype(
    [('npart%i' % i, 'i4') for i in range(6)] +
    [('mass%i' % i, 'f8') for i in range(6)] +
    [('time', 'f8')] +
    [('redshift', 'f8')] +
    [('flag_sfr', 'i4')] +
    [('flag_feedback', 'i4')] +
    [('npartTotal%i' % i, 'i4') for i in range(6)] +
    [('flag_cooling', 'i4')] +
    [('num_files', 'i4')] +
    [('BoxSize', 'f8')] +
    [('Omega0', 'f8')] +
    [('OmegaLambda', 'f8')] +
    [('HubbleParam', 'f8')] +
    [('flag_stellarage', "i4")] +
    [('flag_metals', "i4")] +
    [('npartTotalHighWord%i' % i, 'u4') for i in range(6)] +
    [('flag_entropy_instead_u', "i4")] +
    [('flag_doubleprecision', "i4")] +
    [('flag_lpt_ics', "i4")] +
    [('lpt_scalingfactor', "f4")] +
    [('flag_tracer_field', "i4")] +
    [('composition_vector_length', "i4")] +
    [('fill', 'S40')]
)


def _read_header(infname, icform=1):
    with open(infname, "r") as f:
        f.read(4)
        if (icform == 2):  # Read identifying string for format 2
            f.read(4)
        header = np.fromfile(f, headerType, 1)
        return header


def _read_block(file, dtype, icform=1):
    numBytes = np.fromfile(file, "i4", 1)
    if icform == 2:
        blockname = file.read(4)
    print(np.shape(numBytes), np.shape(dtype.itemsize))
    data = np.fromfile(file, dtype, numBytes / dtype.itemsize)
    numBytes2 = np.fromfile(file, "i4", 1)
    if numBytes != numBytes2:
        if (icform == 2):
            raise ValueError("Byte numbers at beginning and end of block \"%s\" do not match!" % blockname)
        else:
            raise ValueError("Byte numbers at beginning and end of block do not match!")
    return data


def _write_block(file, data, icform=1, blockname=""):
    numBytes = np.array([data.dtype.itemsize * len(data)], "i4")

    numBytes.tofile(file)
    if (icform == 2):
        if len(blockname) != 4:
            raise ValueError('The blockname has to contain exactly 4 characters!')
        f.write(blockname)
    data.tofile(file)
    numBytes.tofile(file)


def _write_unif_block(file, number, value, icform=1, blockname=""):
    numBytes = np.array([value.dtype.itemsize * number], "i4")
    data = value.repeat(number)

    numBytes.tofile(file)

    if (icform == 2):
        if len(blockname) != 4:
            raise ValueError('The blockname has to contain exactly 4 characters!')
        f.write(blockname)

    data.tofile(file)
    numBytes.tofile(file)


def _print_header_simple(header, fname):
    import string
    print("File \"{:s}\" is part of a set of {:d} files\n".format(
        fname.split('/')[-1], int(header['num_files'])))

    print("============= Particle Contents =============\n")
    print("\n".join(["Type {:d} : {:d}\t of {:d}".format(i, int(header['npart%d' % i]),
                                                         int(header['npartTotal%d' % i] + (
                                                                 header['npartTotalHighWord%d' % i].astype(
                                                                     "i8") << 31)))
                     for i in range(6)]))
    print("\n")
    print("\n".join(["Mass {:d} : {:e}".format(i, float(header['mass%d' % i]))
                     for i in range(6)]))
    print("\n")

    print("=============   Box Properties  =============\n")
    print("\n".join(["{:s} : \t\t {:f}".format(string.ljust(name, 16), float(header[name]))
                     for name in ["time", "redshift", "BoxSize", "Omega0", "OmegaLambda", "HubbleParam"]]))
    print("\n")


def subfiletofile(file, outfname=None, icform=1):
    if outfname is None:
        outfname = "out" + file
    with open(outfname, "w") as outf, open(file, "r") as inf:
        inf.seek(0, 2)
        inflen = inf.tell()
        inf.seek(0, 0)
        newHeader = _read_block(inf, headerType, icform=icform)  # read header block from infile to modify
        newHeader["num_files"] = 1
        for i in np.arange(6):
            newHeader["npartTotal%i" % i] = newHeader["npart%i" % i]
        _print_header_simple(newHeader, outfname)
        _write_block(outf, newHeader, icform=icform)  # write changed header
        while not inf.tell() >= inflen:
            print("Copied additional block.")
            additionalblocks = _read_block(inf, np.dtype("S1"), icform=icform)
            _write_block(outf, additionalblocks, icform=icform)
