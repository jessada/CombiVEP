from setuptools import setup
import sys
import glob
import pkgutil

setup(
	name='CombiVEP',
	version='0.1.1',
	author='Jessada Thutkawkorapin',
	author_email='tester@test.com',
	packages=['combivep', 'combivep.engine', 'combivep.engine.test', 'combivep.refdb', 'combivep.refdb.test', 'combivep.preproc', 'combivep.preproc.test', 'combivep.devtools'],
	scripts=['bin/CombiVEP_reference_updater',
             'bin/CombiVEP_predictor',
             'bin/CombiVEP_trainer',
             'bin/dev_report_preproc',
             'bin/dev_demo_training',
             'bin/dev_demo_predicting',
             'bin/dev_report_performance'],
#	scripts=['bin/CombiVEP_reference_updater', 'bin/CombiVEP_predictor', 'bin/CombiVEP_trainer', 'bin/CombiVEP_demo'],
    package=['CombiVEP'],
#    package_dir={'':'combivep'},
    package_data={'': ['data/CBV/*.cbv'],
                  '': ['data/CBV/*.scores'],
                  },
    data_files=[('CombiVEP/CBV', ['combivep/data/CBV/training.cbv', 'combivep/data/CBV/test.cbv', 'combivep/data/CBV/training.cbv.scores', 'combivep/data/CBV/test.cbv.scores']),
                ],
	url='http://pypi.python.org/pypi/combivep/',
	license='LICENSE.txt',
	description='Combined Variant Effect Predictors',
	long_description=open('README.txt').read(),
	install_requires=[
	    "pysam >= 0.7",
	    ],
	)


