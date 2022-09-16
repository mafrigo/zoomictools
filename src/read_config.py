import yaml
from pathlib import Path


def getproperties(parentlabel):
    """
    Function that reads the properties of the parent simulations in use from config.yaml.
    For a given parent label, it returns:

     workdir      : Directory of the position list files and where the output ICs are saved.
     parent_snaps : Name/path of the snapshot files (without the number). If you have multiple subsnapshots, put it
                         in the form ["snapshot_dir","/snapshot_name"] (both without the number).
     parent_ic    : Name/path of the initial condition file of the parent simulation.
     cosmology    : Choice of cosmological parameters. Can be "planck" or "wmap".
     boxsize      : Size of the cosmological box, in comoving Mpc/h.
     zstart       : Initial redshift of the IC.
     seedsset     : Set of seeds. Can be "matteo" or "ludwig", but other options can easily be added in the code.
     parentres    : Resolution level of the parent simulation.
     lowestres    : Lowest resolution level for the zoom-in simulations based on this parent.

    Feel free to add options for more parent sims.
    """
    config_path = Path(__file__).parent / "../config.yaml"
    with open(config_path, "r") as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)
    for section in cfg["parent_setups"]:
        if parentlabel == section:
            config_dict = cfg["parent_setups"][section]
            return config_dict["workdir"], config_dict["parent_snaps"], config_dict["parent_ic"], \
                   config_dict["cosmology"], config_dict["boxsize"], config_dict["zstart"], config_dict["seedsset"], \
                   config_dict["parentres"], config_dict["lowestres"]
    raise IOError("parentlabel " + str(parentlabel) + " does not match any label in config.yaml")
