# DASpec
## Decomposition of AGN Spectrum (DASpec)

This package is forked from [https://github.com/PuDu-Astro/DASpec](https://github.com/PuDu-Astro/DASpec). I have made some changes to make it work on python3 and added some features. Because I cannot install the pybwidget package in python 3 enviroment, I have delete the code of GUI interface and only keep the core code of fitting. If you want to use the GUI interface, please use the original version.

Users could use this package to fit the AGN spectra. I will write a demo to show how to use this package later.

---

A multi-component spectral fitting code for active galactic nuclei.

based on Python 3

Version: 0.0.1

### Dependency (please install them first):
1. GSL
2. swig
3. cmpfit 
(https://www.physics.wisc.edu/~craigm/idl/cmpfit.html). `cmpfit` has been included in this package. You do not need to care about it now.

For Mac users, you can install them with homebrew:
```shell
brew install gsl
brew install swig
```
For Ubuntu or Debian users, you can install them with apt-get:
```shell
sudo apt-get install libgsl-dev
sudo apt-get install swig
```
For arch or manjaro users, you can install them with pacman:
```shell
sudo pacman -S gsl
sudo pacman -S swig
```

### Install:
You can install the package by using pip:
```bash
pip install git+https://github.com/zzxihep/DASpec.git
```
Some system may raise an error `error: externally-managed-environment`, you can try to add `--break-system-packages` to the command above.
```bash
pip install git+https://github.com/zzxihep/DASpec.git --break-system-packages
```
or clone the repository and install it:
```bash
git clone https://github.com/zzxihep/DASpec.git
cd DASpec
python setup.py install --user
```

<!-- ### Install:
1. run "python setup.py build_ext --inplace"
2. Add the path to your $PYTHONPATH
3. Add "DASpec_TEMPLATE_PATH" in .bashrc as the directory of template files, e.g., export DASpec_TEMPLATE_PATH="..."

### Usage:
1. run "DASpec_GUI.py" or "DASpec_GUI_4k.py" (for 4k screen) directly
2. "DASpec_GUI.py -s spectrum_file.txt" (spectrum_file.txt with three columns: wavelength, flux, err)
3. "DASpec_GUI.py -b list.txt" (list.txt is a list of many spectrum files) -->

### Interface:

<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/Interface.png" width="720" height="438">

### Tutorial:

1. Set fitting windows: select the checkboxes and click "Update"

<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/fit_wins.png" width="309" height="325">

2. Set components: add or delete the components, and then select the checkboxes (each component has three lines: input parameters, lower limits and upper limits)

<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/Components.png" width="684" height="583">

Here, I used 

(a) a power law to model the AGN continuum with the parameters of flux and power law index (keyword "5100.0" means the first parameter is the flux at 5100A)

(b) an Fe II template convolved by a Gaussian to model the Fe II emission, with the parameters of flux, width (km/s), and shift (km/s) (keyword "fetemplate_" is the name of the template file which has two columns: wavelength and flux, keywords "4434.0" and "4684.0" mean that the first parameter is integrated flux from 4434.0 and 4684.0)

(c) a double-Gaussians for the broad component of Hbeta emission line with the parameters of flux, width of the first Gaussian, shift of the first Gaussian, width of the second Gaussian, shift of the second Gaussian, and the ratio of the first to the total line flux (keyword "4861.0" means the line center is located at 4861.0A)

(d) a Gaussian for each of the narrow emission lines with the parameters of flux, width, and shift (Hbeta, [OIII]4959,5007). 

3. If you want to tie the profiles of the emission lines: for example, tie component 7 to have the same profile as component 5, and tie the profile of component 6 to 5 with flux ratio 0.3333, then select the checkboxes and click "Update"

<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/fit_tie.png" width="571" height="325">

The narrow Hbeta and [OIII]4959 are contrainted to have the same profile as [OIII]5007. And the flux ratio of [OIII]5007/4959 is fixed to be 3.

4. Begin fitting: Click lmfit (levenberg-Marquardt method) or mix_fit (similar to Basinhopping in Scipy)

<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/fit_buttons.png" width="309" height="145">

5. Check the fitting result

<img src="https://github.com/PuDu-Astro/images_for_doc/blob/master/Spectrum_window.png" width="564" height="513">

6. You will have an .out file in your working directory

7. If you want to begin from the last fitting: DASpec_GUI.py -s spectrum_file.txt -m spectrum_file.txt.out

### Acknowledgement:
I will appreciate if you can cite DASpec in your paper: e.g., \software{\cb DASpec \url{https://github.com/PuDu-Astro/DASpec}}
