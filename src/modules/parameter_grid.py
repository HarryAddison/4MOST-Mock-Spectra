from sklearn.model_selection import ParameterGrid
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.table import QTable


def create_sn_ia_param_grid(config):
    '''
    Form a grid of parameters to form SN spectra from.

    Input:
        > "params" = input parameters

    Output:
        > "param_grid" = gird of parameters
    '''

    # Create the grids for the SALT templates
    args = {
            "template": config["salt_templates"],
            "z": config["z"],
            "t0": config["t0"],
            "amp_passband": config["amp_passband"],
            "x1": config["x1"],
            "c": config["c"],
            "hostebv": [0],
            "zenith": config["zenith"] * u.deg, 
            "seeing": config["seeing"] * u.arcsec,
            "moon_brightness": config["moon_brightness"],
            "t_exp": config["t_exp"] * u.min,
            "ra": config["ra"] * u.deg,
            "dec": config["dec"] * u.deg,
            "amplitude": np.linspace(config["amp_max"], config["amp_min"], config["amp_num"]),
            "observer_phase": np.linspace(config["phase_min"], config["phase_max"], config["phase_num"]) * u.day,
           }
    param_grid = list(ParameterGrid(args))
    
    param_grid = convert_param_grid_to_table(param_grid)

    return param_grid


def convert_param_grid_to_table(param_grid):

    # convert parameter grid from dictionary to dataframe/table format.
    param_grid = pd.DataFrame.from_dict(param_grid, orient="columns")

    # Correcting format of columns that are actually multiple columns.
    template = np.array(list(param_grid["template"])).T
    passband = np.array(list(param_grid["amp_passband"])).T
    # Adding data to the table.
    param_grid = QTable(data=[list(template[0]),
                              list(template[1]),
                              list(template[2]),
                              list(param_grid["z"]),
                              list(param_grid["t0"]),
                              list(passband[0]),
                              list(passband[1]),
                              list(param_grid["x1"]),
                              list(param_grid["c"]),
                              list(param_grid["hostebv"]),
                              list(param_grid["zenith"]) * u.deg,
                              list(param_grid["seeing"]) * u.arcsec,
                              list(param_grid["moon_brightness"]),
                              list(param_grid["t_exp"]) * u.min,
                              list(param_grid["ra"]) * u.deg,
                              list(param_grid["dec"]) * u.deg,
                              list(param_grid["amplitude"]),
                              list(param_grid["observer_phase"])],
                    names=["template_name",
                           "template_version",
                           "sn_type",
                           "z",
                           "t0",
                           "amp_passband",
                           "passband_system",
                           "x1",
                           "c",
                           "hostebv",
                           "zenith",
                           "seeing",
                           "moon_brightness",
                           "t_exp",
                           "ra",
                           "dec",
                           "amplitude",
                           "observer_phase"])

    return param_grid
