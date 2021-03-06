=============
Pre-requisite
=============

* Python 2.7
* pip (http://pypi.python.org/pypi/pip)
* numpy (http://pypi.python.org/pypi/numpy)
* matplotlib (http://matplotlib.org)

============
Installation
============

To install this module, run the following commands:

    sudo python setup.py install

===========================
Download reference database
===========================

After the installation is complete, type

    CombiVEP_reference_updater

This application will automatically check with UCSC and LJB database, and see
if it is required to download the new one. The original database size is
around 1GB each. The total operation time for each database should be
around 30-60 mins.

========
Training
========

After having reference database installed, the CombiVEP model can be trained
using

    CombiVEP_trainer <training_data_file>

<training_data_file> must be in CBV format: CHROM, POS, REF, ALT,
ACTUAL_DELETERIOUS_EFFECT. Each field is separated by a tab. SNP Position(POS)
is 1-based index. The VariBench training data file in CBV format can be found at

    combivep/data/CBV/training.cbv

==========
Prediction
==========

To use the trained model to predict the effect, you can do it using

    CombiVEP_predictor <input_file> [-F FORMAT]

The input file can be either in VCF or above CBV format. Default is in
VCF format. So if you want to use input file in VCF format, simply type

    CombiVEP_predictor <vcf_file>

If you to do the prediction using file in CBV format, you can 

    CombiVEP_predictor <cbv_file> -F CBV

The VariBench test data file in CBV format can be found at

    combivep/data/CBV/test.cbv



