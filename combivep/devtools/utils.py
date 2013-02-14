import os
import sys
import numpy as np
import combivep.devtools.settings as devtools_settings
import combivep.settings as combivep_settings
from combivep.preproc.dataset import DataSetManager
from combivep.cfg import Configure
from combivep.template import CombiVEPBase
from combivep.engine.wrapper import Trainer
from combivep.engine.wrapper import Predictor
from collections import namedtuple
from combivep.devtools.settings import PRECISION_MEASURES

PrecisionPerformance = namedtuple("precision_performance", PRECISION_MEASURES)

class ScoresReader(CombiVEPBase):
    """

    to read parsed SNPS file
    The format are CHROM, POS, REF, ALT, EFFECT.
    All fields are tab separated.

    """


    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, scores_file):
        self.scores_file_name = scores_file

    def fetch_array_snps(self):
        scores_file = open(self.scores_file_name)
        for line in scores_file:
            if line[0] == '#':
                continue
            yield line.rstrip('\n').split('\t')

    def fetch_hash_snps(self):
        for rec in self.fetch_array_snps():
            snp_info = {devtools_settings.KEY_SCORES_CHROM : rec[devtools_settings.SCORES_0_INDEX_CHROM],
                        devtools_settings.KEY_SCORES_POS   : rec[devtools_settings.SCORES_0_INDEX_POS],
                        devtools_settings.KEY_SCORES_REF   : rec[devtools_settings.SCORES_0_INDEX_REF],
                        devtools_settings.KEY_SCORES_ALT   : rec[devtools_settings.SCORES_0_INDEX_ALT],
                        }
            prediction = {devtools_settings.KEY_SCORES_TARGETS : rec[devtools_settings.SCORES_0_INDEX_TARGETS]}
            scores     = {combivep_settings.KEY_PHYLOP_SCORE : rec[devtools_settings.SCORES_0_INDEX_PHYLOP_SCORE],
                          combivep_settings.KEY_SIFT_SCORE   : rec[devtools_settings.SCORES_0_INDEX_SIFT_SCORE],
                          combivep_settings.KEY_PP2_SCORE    : rec[devtools_settings.SCORES_0_INDEX_PP2_SCORE],
                          combivep_settings.KEY_LRT_SCORE    : rec[devtools_settings.SCORES_0_INDEX_LRT_SCORE],
                          combivep_settings.KEY_MT_SCORE     : rec[devtools_settings.SCORES_0_INDEX_MT_SCORE],
                          combivep_settings.KEY_GERP_SCORE   : rec[devtools_settings.SCORES_0_INDEX_GERP_SCORE],
                          }
            yield {combivep_settings.KEY_SNP_INFO_SECTION   : snp_info,
                   combivep_settings.KEY_PREDICTION_SECTION : prediction,
                   combivep_settings.KEY_SCORES_SECTION     : scores,
                   }

class FastDataSetManager(DataSetManager):


    def __init__(self, config_file=combivep_settings.COMBIVEP_CONFIGURATION_FILE):
        DataSetManager.__init__(self, config_file=combivep_settings.COMBIVEP_CONFIGURATION_FILE)

    def load_data(self, file_name, file_type=devtools_settings.FILE_TYPE_SCORES):
        if file_type == devtools_settings.FILE_TYPE_SCORES:
            return self.__load_scores_data(file_name)
        else:
            return DataSetManager.load_data(self, file_name, file_type)

    def __load_scores_data(self, file_name):
        self.clear_data()
        scores_reader = ScoresReader()
        scores_reader.read(file_name)
        for rec in scores_reader.fetch_hash_snps():
            snp_data   = {combivep_settings.KEY_CHROM : rec[combivep_settings.KEY_SNP_INFO_SECTION][devtools_settings.KEY_SCORES_CHROM],
                          combivep_settings.KEY_POS   : rec[combivep_settings.KEY_SNP_INFO_SECTION][devtools_settings.KEY_SCORES_POS],
                          combivep_settings.KEY_REF   : rec[combivep_settings.KEY_SNP_INFO_SECTION][devtools_settings.KEY_SCORES_REF],
                          combivep_settings.KEY_ALT   : rec[combivep_settings.KEY_SNP_INFO_SECTION][devtools_settings.KEY_SCORES_ALT],
                          }
            prediction = {combivep_settings.KEY_TARGETS : rec[combivep_settings.KEY_PREDICTION_SECTION][devtools_settings.KEY_SCORES_TARGETS]}
            self.dataset.append({combivep_settings.KEY_SNP_INFO_SECTION : snp_data,
                                 combivep_settings.KEY_PREDICTION_SECTION : prediction,
                                 combivep_settings.KEY_SCORES_SECTION : rec[combivep_settings.KEY_SCORES_SECTION],
                                 })



def filter_cbv_data(cbv_file,
                    config_file=combivep_settings.COMBIVEP_CONFIGURATION_FILE):
    (dir_name, file_name) = os.path.split(cbv_file)
    print
    print "> > >  " + file_name
    dm = DataSetManager()
    dm.load_data(cbv_file, file_type=combivep_settings.FILE_TYPE_CBV)
    print "%-25s: %5d\n" % ("Original", len(dm.dataset))

    dm.validate_data()
    f_clean = open(cbv_file + '.clean', 'w')
    for item in dm.dataset:
        f_clean.write("%s\t%s\t%s\t%s\t%s\n" % (item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_CHROM],
                                                 item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_POS],
                                                 item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_REF],
                                                 item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_ALT],
                                                 item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS],
                                                 )
                        )
    f_clean.close()
    print "%-25s: %5d" % ("Clean pathogenic", len([item for item in dm.dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Clean neutral", len([item for item in dm.dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '0']))
    print "%-25s: %5d\n" % ("Total", len(dm.dataset))

    dm.calculate_scores()
    f_scores = open(cbv_file + '.scores', 'w')
    for item in dm.dataset:
        f_scores.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_CHROM],
                                                                         item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_POS],
                                                                         item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_REF],
                                                                         item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_ALT],
                                                                         item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS],
                                                                         item[combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_PHYLOP_SCORE],
                                                                         item[combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_SIFT_SCORE],
                                                                         item[combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_PP2_SCORE],
                                                                         item[combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_LRT_SCORE],
                                                                         item[combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_MT_SCORE],
                                                                         item[combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_GERP_SCORE],
                                                                         )
                        )
    f_scores.close()
    print "%-25s: %5d" % ("Scored pathogenic", len([item for item in dm.dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Scored neutral", len([item for item in dm.dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '0']))
    print "%-25s: %5d\n" % ("Total", len(dm.dataset))

    dm.set_shuffle_seed(combivep_settings.DEMO_SEED)
    dm.shuffle_data()
    dm.partition_data()

    #partition data
    training_dataset      = dm.get_training_data() 
    validation_dataset    = dm.get_validation_data() 

    print "%-25s: %5d" % ("Training pathogenic", len([item for item in training_dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Training neutral", len([item for item in training_dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '0']))
    print "%-25s: %5d\n" % ("Total", len(training_dataset))

    print "%-25s: %5d" % ("Validation pathogenic", len([item for item in validation_dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Validation neutral", len([item for item in validation_dataset if item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS] == '0']))
    print "%-25s: %5d\n" % ("Total", len(validation_dataset))

def calculate_roc(pathogenic_dataset, neutral_dataset, roc_range):
    false_positive_rates = np.zeros([len(roc_range), pathogenic_dataset.shape[1]])
    true_positive_rates  = np.zeros([len(roc_range), neutral_dataset.shape[1]])
    pathogenic_data_size = pathogenic_dataset.shape[0]
    neutral_data_size    = neutral_dataset.shape[0]
    for i in xrange(len(roc_range)):
        false_positive_rates[i, :] = np.matrix(np.sum(neutral_dataset > roc_range[i], axis=0).astype(np.float))/ neutral_data_size
        true_positive_rates[i, :]  = np.matrix(np.sum(pathogenic_dataset > roc_range[i], axis=0).astype(np.float))/ pathogenic_data_size
    return (false_positive_rates, true_positive_rates)

def print_preproc(col1, col2, col3):
    print "%-40s%18s%18s" % (col1, col2, col3)

def print_precision(col1, col2, col3):
    print "%-25s%9s%10s" % (col1, col2, col3)

def fast_training(training_data_file, 
                  params_out_file=combivep_settings.USER_PARAMETERS_FILE,
                  random_seed=combivep_settings.DEFAULT_SEED,
                  n_hidden_nodes=combivep_settings.DEFAULT_HIDDEN_NODES,
                  figure_dir=combivep_settings.DEFAULT_FIGURE_DIR,
                  iterations=combivep_settings.DEFAULT_ITERATIONS,
                  config_file=combivep_settings.COMBIVEP_CONFIGURATION_FILE,
                  ):
    """

    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
    CBV has 5 fields, CHROM, POS, REF, ALT, EFFECT (1=deleterious, 0=neutral). All are tab separated
    Required arguments
    - neutral_data_file : list of SNPs with no harmful effect, CBV format
    - pathognice_data_file : list of SNPs with deleterious effect, CBV format

    """
    #pre-processing dataset
    print >> sys.stderr, 'pre-processing dataset, this may take a while (around 750 SNPs/mins). . . . '
    dm = FastDataSetManager(config_file=config_file)
    dm.load_data(training_data_file, file_type=devtools_settings.FILE_TYPE_SCORES)
    dm.set_shuffle_seed(random_seed)
    dm.shuffle_data()
    dm.partition_data()

    #partition data
    training_dataset      = dm.get_training_data() 
    validation_dataset    = dm.get_validation_data() 

    #train !!!
    print >> sys.stderr, 'Training CombiVEP, please wait (around 500 SNPs/mins) . . . . '
    trainer = Trainer(training_dataset, validation_dataset, random_seed, n_hidden_nodes, figure_dir)
    trainer.train(iterations)
    if not os.path.exists(combivep_settings.USER_PARAMETERS_DIR):
        os.makedirs(combivep_settings.USER_PARAMETERS_DIR)
    trainer.export_best_parameters(params_out_file)

def fast_predict(SNPs_file,
                 params_file=combivep_settings.USER_PARAMETERS_FILE,
                 file_type=combivep_settings.FILE_TYPE_VCF,
                 output_file=None,
                 config_file=combivep_settings.COMBIVEP_CONFIGURATION_FILE,
                 ):
    """

    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
    CBV has 5 fields, CHROM, POS, REF, ALT, EFFECT (1=deleterious, 0=neutral). All are tab separated
    Required arguments
    - SNPs_file : list of SNPs to be predicted, can be either VCF or CBV (default is VCF)

    """
    #pre-processing test dataset
    print >> sys.stderr, 'pre-processing dataset, this may take a while (around 750 SNPs/mins). . . . '
    dm = FastDataSetManager(config_file=config_file)
    dm.load_data(SNPs_file, file_type=devtools_settings.FILE_TYPE_SCORES)

    #predict
    predictor = Predictor()
    predictor.import_parameters(params_file=params_file)
    out = (np.array(predictor.predict(dm.dataset)).reshape(-1,))

    #print output
    if output_file is not None:
        sys.stdout = open(output_file, 'w')
    print >> sys.stdout, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("CHROM", "POS", "REF", "ALT", "ACTUAL_DELETERIOUS_EFFECT", "PREDICTED_DELETERIOUS_PROBABILITY", "PHYLOP_SCORE", "SIFT_SCORE", "PP2_SCORE", "LRT_SCORT", "MT_SCORE", "GERP_SCORE")
    for i in xrange(len(dm.dataset)):
        print >> sys.stdout, "%s\t%s\t%s\t%s\t%s\t%6.4f\t%s\t%s\t%s\t%s\t%s\t%s" % (dm.dataset[i][combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_CHROM],
                                                                                    dm.dataset[i][combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_POS],
                                                                                    dm.dataset[i][combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_REF],
                                                                                    dm.dataset[i][combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_ALT],
                                                                                    dm.dataset[i][combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS],
                                                                                    out[i],
                                                                                    dm.dataset[i][combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_PHYLOP_SCORE],
                                                                                    dm.dataset[i][combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_SIFT_SCORE],
                                                                                    dm.dataset[i][combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_PP2_SCORE],
                                                                                    dm.dataset[i][combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_LRT_SCORE],
                                                                                    dm.dataset[i][combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_MT_SCORE],
                                                                                    dm.dataset[i][combivep_settings.KEY_SCORES_SECTION][combivep_settings.KEY_GERP_SCORE],
                                                                                    )
    sys.stdout = sys.__stdout__

def measure_precision(ratio,
                      positive_samples_scores,
                      negative_samples_scores,
                      ):
    true_positive  = len(positive_samples_scores[positive_samples_scores > ratio])
    false_negative = len(positive_samples_scores[positive_samples_scores <= ratio])
    true_negative  = len(negative_samples_scores[negative_samples_scores < ratio])
    false_positive = len(negative_samples_scores[negative_samples_scores >= ratio])

    accuracy         = float(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative)
    sensitivity      = float(true_positive)/(true_positive+false_negative)
    specificity      = float(true_negative)/(true_negative+false_positive)
    balance_accuracy = (sensitivity+specificity)/2

    return PrecisionPerformance(true_positive, false_negative, true_negative, false_positive, accuracy, sensitivity, specificity, balance_accuracy)
