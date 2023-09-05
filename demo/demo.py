import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


def read_spec(fitname):
    hdul = fits.open(fitname)
    data = hdul[1].data
    wave = 10**data['loglam'].astype(float)
    flux = data['flux'].astype(float)
    ivar = data['ivar'].astype(float)
    arg = ivar > 0
    err = np.sqrt(1.0 / ivar[arg])
    wave = wave[arg]
    flux = flux[arg]
    return wave, flux, err


def get_redshift(fitname):
    hdul = fits.open(fitname)
    data = hdul[2].data
    return float(data['Z'][0])


def main():
    fitname = 'spec0627-52144-0143.fits'
    wave, flux, err = read_spec(fitname)
    redshift = get_redshift(fitname)
    wave_rest = wave / (1 + redshift)
    plt.plot(wave_rest, flux)
    plt.show()


if __name__ == '__main__':
    main()