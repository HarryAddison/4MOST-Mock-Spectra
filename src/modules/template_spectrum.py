import sncosmo as sn
import numpy as np
import astropy.units as u
from astropy.table import QTable


def template_spectrum(params, save_path):

    model = _sn_model(params)
    spectrum = _spectrum(model, params, save_path)

    return spectrum


def _sn_model(params):

    '''
    Produce the model of the SN.

    Model can be used to produce the SN spectra at a given
    time, or it can be used to produce synthetic photometry.

    Model includes effects of Galactic and extragalactic extinction.
    '''

    # Setting the model parameters
    source = sn.get_source(params["template_name"], params["template_version"])

    # Setting up the model and inclusion of dust (host+Galactic).
    dust = sn.F99Dust()

    model = sn.Model(source=source,
                        effects=[dust, dust],
                        effect_names=["host", "mw"],
                        effect_frames=["rest", "obs"])

    # Set the parameters of the template spectrum
    for param in model.param_names:
        if param != "amplitude" and param != "x0":
            model[param] = params[param]

    # Set the absoulte amplitude of the model spectrum
    # Must be done after other params as SALT x1 and c alter it.
    model.set_source_peakabsmag(params["amplitude"],
                                params["amp_passband"],
                                params["passband_system"])

    return model


def _spectrum(model, params, save_path):
    '''
    Use the SN model (self.model) to produce fluxes over for the
    given wavelength range.

    Fluxes produced include effects of Galactic and extragalactic
    extinction.
    '''

    # Create wavelengths of spectrum. Resoultion of 4 per angstrom
    # in the wavelength range of the model.
    res = 4 * int(model.maxwave() - model.minwave())
    wl = (np.linspace(model.minwave(), model.maxwave(), res, endpoint=True)) * u.angstrom

    # Produce fluxes at the given wavelengths.
    flux = (model.flux(params["observer_phase"], wl) * (u.erg / u.s / (u.cm**2) / u.angstrom))

    # Putting spectrum data into a table.
    spec = QTable([wl, flux], names=("wave", "flux"))

    # Check to make sure all fluxes are > 0. If not then shift the flux
    # of the entire spectrum.
    min_flux = min(flux)
    if min_flux <= 0:
            spec["flux"] -= min_flux  # subtract the -ve flux to make it positive
            spec["flux"] += 10**-200 * spec["flux"].unit  # Add a very small number to make zero values non zero

    # Writing the spectrum data to file
    spec.write(save_path, format="fits", overwrite=True)

    return spec
