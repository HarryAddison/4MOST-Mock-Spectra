from qmostetc import Spectrum, QMostObservatory, L1DXU, Atmosphere
import numpy as np
import astropy.units as u
from modules.snr import calc_snr


def l1_spectrum(template_spec, params, save_path):
    spec = _l1_spectrum(template_spec, params, save_path)

    snr, snr_4500_8000 = calc_snr(spec, n_bins=None)
    snr_per_15A, snr_4500_8000_per_15A = calc_snr(spec)

    return spec, [snr, snr_4500_8000, snr_per_15A, snr_4500_8000_per_15A]


def _dummy_atmosphere():
    wavelength = np.arange(300, 950, 10) * u.nm
    zenith_transmission = np.ones(len(wavelength))
    flux_unit = u.ph / (u.s * u.m**2 * u.arcsec**2 * u.nm) 
    dark_sky_flux = np.zeros((len(wavelength), 2)) * flux_unit
    moon_flux = np.zeros((len(wavelength), 2)) * flux_unit
    airmasses = np.array([1., 100.])
    def moon_sun_coeff(sep, wl): return np.zeros(len(wl))
    return Atmosphere(wavelength, zenith_transmission, dark_sky_flux,
                      moon_flux, airmasses, moon_sun_coeff)


def _l1_spectrum(template_spec, params, save_path):

    '''
    From the given template spectrum produced the spectrum as if
    it were observed by 4MOST with the given parameters
    (i.e sky position, moon brightness, etc.).
    '''

    # put template spectrum into 4MOST ETC Spectrum object
    template_spec = Spectrum(template_spec["wave"], template_spec["flux"])

    # Remove the atomosphere
    observatory = QMostObservatory('lrs')
    observatory.atmosphere = _dummy_atmosphere()
    for arm_observatory in observatory.values():
        arm_observatory.atmosphere = observatory.atmosphere

    obs = observatory(params["zenith"], params["seeing"], params["moon_brightness"])
    obs.set_target(template_spec, 'point')
    tbl = obs.expose(params["t_exp"])
    dxu = L1DXU(observatory, tbl, params["t_exp"], with_noise=True)

    l1_spec = dxu.joined_spectrum()

    l1_spec.rename_column("WAVE", "wave")
    l1_spec.rename_column("FLUX", "flux")
    l1_spec.rename_column("ERR_FLUX", "flux_err")
    l1_spec.write(save_path, format="fits", overwrite=True)

    return l1_spec

