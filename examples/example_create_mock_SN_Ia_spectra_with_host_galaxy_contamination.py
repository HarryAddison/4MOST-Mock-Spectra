from qmost_mock_spectra import *
import numpy as np
import astropy.units as u
from astropy.table import QTable
from scipy import interpolate
import os
from qmostetc import Spectrum


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
    param_grid["gal_type"] = "None"
    param_grid["gal_mag_sn_mag_frac"] = -1.0
    for i, params in enumerate(param_grid):

        # SN spectrum
        sn_template_path = f"{config['sn_template_spec_dir']}/template_{params['id']}.fits"
        sn_template_spec = template_spectrum(params, sn_template_path)
        sn_template_spec["wave"] = sn_template_spec["wave"].to(u.nm)

        sn_spec = Spectrum(sn_template_spec["wave"], sn_template_spec["flux"].value * u.erg / (u.cm ** 2 * u.s * u.angstrom))
        sn_mag = sn_spec.get_mag(u.ABmag, "LSST_LSST.r")

        # Galaxy parameters and spectrum
        gal_type = assign_host(config["host_type"])

        gal_template_path = f"{config['gal_template_spec_dir']}/{gal_type}_template.fits"
        gal_template_spec = QTable.read(gal_template_path)
        gal_template_spec["WAVELENGTH"] *= u.angstrom
        gal_template_spec["FLUX"] *= u.erg / (u.cm ** 2 * u.s * u.angstrom)

        gal_template_spec["WAVELENGTH"] = gal_template_spec["WAVELENGTH"].to(u.nm)
        gal_template_spec["WAVELENGTH"] *= (1 + params["z"])  # template is in rest frame so redshift it.
        
        gal_spec = Spectrum(gal_template_spec["WAVELENGTH"], gal_template_spec["FLUX"].value * u.erg / (u.cm ** 2 * u.s * u.angstrom))
        gal_mag = gal_spec.get_mag(u.ABmag, "LSST_LSST.r")

        gal_mag_required = sn_mag - 2.5 * np.log10(config["sn_gal_flux_ratio"]) * u.ABmag  # Calculate magnitude if its x times brighter than SN mag
        
        flux_ratio = 10**((gal_mag_required - gal_mag).value/ -2.5)  # Scale the spectrum to the desired magnitude
        gal_template_spec["FLUX"] *= flux_ratio

        # Interpolate the galaxy spectrum to the SN spectrum wavelength values
        spl = interpolate.splrep(gal_template_spec["WAVELENGTH"].value, gal_template_spec["FLUX"].value)
        gal_spec_flux_interpolated = interpolate.splev(sn_template_spec["wave"].value, spl) * gal_template_spec["FLUX"].unit

        # Contaminate the SN spectrum with the galaxy flux (linearly)
        sn_template_spec["flux"] += gal_spec_flux_interpolated

        # Create L1 spectrum
        l1_path = f"{config['l1_spec_dir']}/l1_spectrum_{params['id']}.fits"
        spec, snr_values = l1_spectrum(sn_template_spec, params, l1_path, "point")

        # add the SNR and galaxy parameters to the parameter grid
        param_grid["snr"][i] = snr_values[0]
        param_grid["snr_4500_8000"][i] = snr_values[1]
        param_grid["snr_per_15A"][i] = snr_values[2]
        param_grid["snr_4500_8000_per_15A"][i] = snr_values[3]
        param_grid["gal_type"][i] = gal_type
        param_grid["gal_mag_sn_mag_frac"][i] = config["sn_gal_flux_ratio"]

    # Save the parameters grid of the spectra to produced.
    param_grid.write(config["param_grid_path"], overwrite=True)

    print("Full parameter grid saved.")
