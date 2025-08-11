from modules.input_output import load_config
from modules.parameter_grid import create_sn_ia_param_grid
from modules.extinction import get_mwebv
from modules.template_spectrum import template_spectrum
from modules.l1_spectrum import l1_spectrum
import numpy as np
import astropy.units as u


if __name__ == "__main__":

    config = load_config()

    param_grid = create_sn_ia_param_grid(config)

    param_grid["id"] = np.arange(1, 1+ len(param_grid), 1)
    param_grid["rest_phase"] = param_grid["observer_phase"] / (1 + param_grid["z"]) * u.day

    param_grid["mwebv"] = 0.0 * u.mag
    for i, params in enumerate(param_grid):
        param_grid["mwebv"][i] = get_mwebv(params["ra"], params["dec"]) * u.mag

    param_grid["snr"] = -1.0
    param_grid["snr_4500_8000"] = -1.0
    param_grid["snr_per_15A"] = -1.0
    param_grid["snr_4500_8000_per_15A"] = -1.0
    for i, params in enumerate(param_grid):

        # Set the save paths for figures and data#
        template_path = f"{config['template_spec_dir']}/template_{params['id']}.fits"
        template_spec = template_spectrum(params, template_path)

        l1_path = f"{config['l1_spec_dir']}/l1_spectrum_{params['id']}.fits"
        spec, snr_values = l1_spectrum(template_spec, params, l1_path)

        # add the SNR to the parameter grid
        param_grid["snr"][i] = snr_values[0]
        param_grid["snr_4500_8000"][i] = snr_values[1]
        param_grid["snr_per_15A"][i] = snr_values[2]
        param_grid["snr_4500_8000_per_15A"][i] = snr_values[3]

    # Save the parameters grid of the spectra to produced.
    param_grid.write(config["param_grid_path"], overwrite=True)

    print("Full parameter grid saved.")
