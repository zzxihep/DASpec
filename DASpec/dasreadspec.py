import os
import numpy as np
from PyAstronomy import pyasl
from scipy.interpolate import interp1d
# import matplotlib.pyplot as plt
from PyAstronomy.pyaC.pyaErrors.pyaValErrs import PyAValError
from astropy.io import fits
from specutils import Spectrum1D


def get_airmass(fn):
    try:
        air = fits.getval(fn, 'AIRMASS')
        if air > 3:
            air=3.0
        return air
    except Exception:
        return 1.0


def get_snr(flux):
    signal = np.median(flux)
    n = len(flux)
    noise = 0.6052697 * np.median(np.abs(2*flux[2:n-2] - flux[0:n-4] - flux[4:n]))
    snr = signal / noise
    return snr


# def get_err(wave, flux, redshift, airmass):
#     snr = get_snr(flux)
#     rat = np.sqrt(np.median(flux)) / snr
#     err = np.sqrt(flux) * rat
#     return err


def get_err(wave, flux, redshift, airmass):
    print('redshift = ' + str(redshift))
    curdir = os.path.dirname(os.path.realpath(__file__))
    extfile = 'kpnoextinct.dat'
    dataext = np.loadtxt(curdir+os.sep+extfile)
    wext = dataext[:, 0]
    magext = dataext[:, 1] * airmass
    fscale = 10**(magext/2.5)
    # funext = interp1d(wext, fscale, kind='cubic')
    scale = np.interp(wave, wext, fscale)
    # scale = funext(wave)
    # print scale[0], scale[-1]
    counts = wave*(flux/scale)
    cerr = 1 / np.sqrt(counts)
    # plt.plot(wave, cerr)
    # plt.show()
    rerr = flux * cerr
    conl = (4740-15)*(1+redshift)
    conr = (4740+15)*(1+redshift)
    arg = np.where((wave>conl)&(wave<conr))
    # print '-' * 44
    # partw = wave[arg]
    partf = flux[arg]
    partre = rerr[arg]
    reale = np.std(partf)
    rat = reale / np.median(partre)
    # print rat
    err = rerr * rat
    return err


def is_sdss_spec(fn):
    try:
        key = fits.getval(fn, 'TELESCOP').lower()
        if 'sdss' in key:
            return True
        else:
            return False
    except KeyError:
        return False


def read_fits_spec(fn, redshift=0):
    try:
        if is_sdss_spec(fn):
            spec = Spectrum1D.read(fn)
            wave = spec.wavelength.value.astype('float64')
            flux = spec.flux.value.astype('float64')
            ivar = spec.uncertainty.array
            arg = np.where(ivar == 0)
            arg2 = np.where(ivar != 0)
            res2 = ivar[arg2]
            ivar2 = np.min(res2) * 0.01
            ivar[arg] = ivar2
            err = 1/ivar
            err = err.astype('float64')
            return wave, flux, err
        else:
            wave, flux = pyasl.read1dFitsSpec(fn)
        # arg = np.where((wave<8000) & (wave>3500))
        # wave = wave[arg]
        # flux = flux[arg]
        # rerr = np.sqrt(wave*flux)
        # ind = int(len(wave)/2)
        # con = 0.01*flux[ind]/rerr[ind]
        # err = con*rerr
        # redshift = 0.158339
        # redshift = 0.0
        # redshift = 0.022189
        # print '-'* 44
        airmass = get_airmass(fn)
        # err = np.ones(flux.shape) * np.mean(flux) * 0.01
        # print '-'* 44
        err = get_err(wave, flux, redshift, airmass)
        return wave, flux, err
    except Exception:
        hl = fits.open(fn)
        naxis1 = hl[0].header['NAXIS1']
        ref = hl[0].header['CRVAL1']
        step = hl[0].header['CD1_1']
        naxis = hl[0].header['NAXIS']
        if naxis == 2:
            flux = hl[0].data[0]
            err = hl[0].data[-1]
        elif naxis == 3:
            flux = hl[0].data[0, -1]
            err = hl[0].data[-1, -1]
        wave = ref+np.arange(naxis1)*step
        return wave, flux, err


def readspec(fn, redshift=0):
    filetype = os.path.splitext(fn)[1]
    if filetype == '.fits':
        return read_fits_spec(fn)
    else:
        data = np.loadtxt(fn)
        wave = data[:, 0]
        flux = data[:, 1]
        # import matplotlib.pyplot as plt
        # plt.plot(wave, flux)
        # plt.show()
        if data.shape[1] == 3:
            err = data[:, 2]
        else:
            err = get_err(wave, flux, redshift, 1)
        return wave, flux, err
