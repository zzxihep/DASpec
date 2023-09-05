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
os.chdir(curdir)
setup()
