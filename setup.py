import os
from distutils.core import setup

setup(
    name='sbml_generalization',
    description='Knowledge-based generalization for metabolic models in SBML format.',
    long_description=open('README.md').read(),
    author='Anna Zhukova',
    author_email='zhutchok@gmail.com',
    url='https://github.com/annazhukova/mod_gen',
    version='0.1.1',
    packages=['sbml_generalization'],
    platform=['MacOS', 'Linux', 'Windows'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    package_data={'sbml_generalization': [os.path.join('generalization', '*.py'),
                                          os.path.join('merge', '*.py'),
                                          os.path.join('runner', '*.py'),
                                          os.path.join('sbml', '*.py'),
                                          os.path.join('..', 'README.md')]},
    include_package_data=True,
    download_url='https://github.com/annazhukova/mod_gen/archive/0.1.1.zip',
    install_requires=['python-libsbml-experimental', 'mod_sbml']
)
