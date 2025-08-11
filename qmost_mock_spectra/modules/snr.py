import numpy as np

def convolveTo15A(x,nBins=60):
    return np.convolve(x, np.ones(nBins), 'valid') / nBins


def convolveNoiseTo15A(x,nBins=60):
    for m in range(len(x)-(nBins-1)):
        yield np.sqrt(sum((np.ones(nBins) * x[m:m+nBins])**2)) / nBins


def calc_snr(spec, n_bins=60):
    '''
    Calculate the mean SNR of the given spectrum data.

    As qmostetc spectra are binned as 0.25A, and we want to bin the
    spectrum per 15A, then 60 bins are needed, which is the default
    value.
    '''

    # convert data to dimensionless quantities
    wave = spec["wave"].value
    signal = spec["flux"].value
    noise = spec["flux_err"].value

    if n_bins != None:
        # rebin the data
        # Convert fluxes to flux * wave before binning then convert back.
        wave15 = convolveTo15A(wave, n_bins)
        signal = convolveTo15A(signal * wave, n_bins) / wave15
        noise = np.array(list(convolveNoiseTo15A(noise * wave, n_bins))) / wave15
        wave = wave15

    # Mean SNR of entire spectrum.
    mean_snr_spec = np.mean(signal/noise)
    
    # Mean SNR in range 4500A to 8000A
    snrMask = (wave>=4500) & (wave<=8000)
    mean_snr_4500_8000 = np.mean(signal[snrMask]/noise[snrMask])

    return mean_snr_spec, mean_snr_4500_8000
