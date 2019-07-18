<<<<<<< HEAD
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
 
setuptools.setup(
    name='eproj',
    packages=setuptools.find_packages(), # Mismo nombre que en la estructura de carpetas de arriba
    version='0.1',
    license='LGPL v3', # La licencia que tenga tu paqeute
    description='A random test lib',
    author='Esther',
    author_email='esvemira@hotmail.com',
    url='https://github.com/sthhher/eproj', # Usa la URL del repositorio de GitHub
    keywords='test example develop', # Palabras que definan tu paquete
    classifiers=['Programming Language :: Python'],
=======
from distutils.core import setup
setup(
    name = 'mypackage',
    packages = ['mypackage'], # this must be the same as the name above
    version = '0.1',
    description = 'Pass the file of github to pipy',
    author = 'Esther Vendrell',
    author_email = 'esvemira@hotmail.com',
    url = 'https://github.com/sthhher/eproj.git', # use the URL to the github repo
    download_url = 'https://github.com/sthhher/eproj.git/tarball/0.1',
    keywords = ['testing', 'logging', 'example'],
    classifiers = [],
>>>>>>> 385e2b8d45c70a2f7876ab852418c3344cc9b2c0
)
