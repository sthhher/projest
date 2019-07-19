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
    url='https://github.com/sthhher/projest', # Usa la URL del repositorio de GitHub
    keywords='test example develop', # Palabras que definan tu paquete
    classifiers=['Programming Language :: Python'],
)
