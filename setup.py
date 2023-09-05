import os
from setuptools import setup


curdir = os.getcwd()
sdir = os.path.dirname(os.path.abspath(__file__))
mpdir = os.path.join(sdir, 'mpfit/cmpfit-1.3a')
os.chdir(mpdir)
os.system('make')
os.chdir(curdir)
DAS_dir = os.path.join(sdir, 'DASpec')
os.chdir(DAS_dir)
os.system('make')

if os.path.exists('_swigDASpec.so') == False:
    raise RuntimeError('_swigDASpec.so not found, please check the makefile')
if os.path.exists('_carray.so') == False:
    raise RuntimeError('_carray.so not found, please check the makefile')

os.chdir(curdir)
setup()
