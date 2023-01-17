from setuptools import setup, find_packages

setup(name='spontaneous_gamma_circuit',
      version='0.1',
      description='AR(2) model to linear E-I circuit',
      url='https://github.com/matteosaponati/spontaneous-gamma-circuit.git',
      author='Matteo Saponati',
      author_email='matteosaponati@gmail.com',
      packages=find_packages(),
      install_requires=[
        'numpy'
      ],
      zip_safe=False,
      )