from setuptools import setup
import sys
import glob
import pkgutil


setup(
    name='CombiVEP',
    version='0.1.2',
    author='Jessada Thutkawkorapin',
    author_email='jessada.thutkawkorapin@gmail.com',
    packages=['combivep',
              'combivep.engine',
              'combivep.refdb',
              'combivep.preproc',
              ],
    scripts=['bin/CombiVEP_reference_updater',
             'bin/CombiVEP_predictor',
             'bin/CombiVEP_trainer',
             ],
    package=['CombiVEP'],
    package_data={'': ['data/CBV/*.cbv']
                  },
    data_files=[('combivep/data/CBV', ['combivep/data/CBV/training.cbv',
                                       'combivep/data/CBV/test.cbv',
                                       ]),
                ],
    url='http://pypi.python.org/pypi/CombiVEP/',
    license='LICENSE.txt',
    description='Combined Variant Effect Predictors',
    long_description=open('README.txt').read(),
    install_requires=["pysam >= 0.7"],
)
