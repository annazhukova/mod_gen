__author__ = 'anna'

from setuptools import setup, find_packages
from sys import version

if version < '2.2.3':
    from distutils.dist import DistributionMetadata

    DistributionMetadata.classifiers = None
    DistributionMetadata.download_url = None


setup(name='SBMLModelGeneralization',
      description='SBML model generalization utilities.',
      long_description=open('README.txt').read(),
      author='Anna Zhukova',
      author_email='zhutchok@gmail.com',
      url='http://metamogen.gforge.inria.fr',
      version='1.0',
      packages=find_packages(exclude=['sbml_generalization/runner/', 'sbml_generalization/runner/path2models']),
      include_package_data=True,
      platform=['MacOS', 'Linux', 'Windows'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      download_url='http://metamogen.gforge.inria.fr',
      requires=['python-libsbml-experimental', 'mod_sbml']
)