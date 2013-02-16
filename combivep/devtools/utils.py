import os
import sys
import numpy as np
import combivep.devtools.settings as dev_const
import combivep.settings as cbv_const
from combivep.preproc.dataset import DataSetManager
from combivep.config import Configure
from combivep.template import CombiVEPBase
from combivep.engine.wrapper import Trainer
from combivep.engine.wrapper import Predictor
from collections import namedtuple
from combivep.devtools.settings import PRECISION_MEASURES

PrecisionPerformance = namedtuple("precision_performance", PRECISION_MEASURES)

class ScoresRecord(CombiVEPBase):
    """ to automatically parse data with scores """


    def __init__(self, array_snp):
        self.__array_snp = array_snp

    @property
    def chrom(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_CHROM]

    @property
    def pos(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_POS]

    @property
    def ref(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_REF]

    @property
    def alt(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_ALT]

    @property
    def targets(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_TARGETS]

    @property
    def phylop_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_PHYLOP_SCORE]

    @property
    def sift_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_SIFT_SCORE]

    @property
    def pp2_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_PP2_SCORE]

    @property
    def lrt_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_LRT_SCORE]

    @property
    def mt_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_MT_SCORE]

    @property
    def gerp_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_GERP_SCORE]

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
            yield ScoresRecord(line.rstrip('\n').split('\t'))

    def fetch_hash_snps(self):
        for rec in self.fetch_array_snps():
            snp_info = {dev_const.KEY_SCORES_CHROM : rec.chrom,
                        dev_const.KEY_SCORES_POS   : rec.pos,
                        dev_const.KEY_SCORES_REF   : rec.ref,
                        dev_const.KEY_SCORES_ALT   : rec.alt,
                        }
            prediction = {dev_const.KEY_SCORES_TARGETS : rec.targets}
            scores     = {cbv_const.KEY_PHYLOP_SCORE : rec.phylop_score,
                          cbv_const.KEY_SIFT_SCORE   : rec.sift_score,
                          cbv_const.KEY_PP2_SCORE    : rec.pp2_score,
                          cbv_const.KEY_LRT_SCORE    : rec.lrt_score,
                          cbv_const.KEY_MT_SCORE     : rec.mt_score,
                          cbv_const.KEY_GERP_SCORE   : rec.gerp_score,
                          }
            yield {cbv_const.KEY_SNP_INFO_SECTION   : snp_info,
                   cbv_const.KEY_PREDICTION_SECTION : prediction,
                   cbv_const.KEY_SCORES_SECTION     : scores,
                   }

class FastDataSetManager(DataSetManager):


    def __init__(self, cfg_file=cbv_const.CBV_CFG_FILE):
        DataSetManager.__init__(self, cfg_file=cbv_const.CBV_CFG_FILE)

    def load_data(self, file_name, file_type=dev_const.FILE_TYPE_SCORES):
        if file_type == dev_const.FILE_TYPE_SCORES:
            return self.__load_scores_data(file_name)
        else:
            return DataSetManager.load_data(self, file_name, file_type)

    def __load_scores_data(self, file_name):
        self.clear_data()
        scores_reader = ScoresReader()
        scores_reader.read(file_name)
        for rec in scores_reader.fetch_hash_snps():
            snp_info = rec[cbv_const.KEY_SNP_INFO_SECTION]
            snp_data = {cbv_const.KEY_CHROM : snp_info[dev_const.KEY_SCORES_CHROM],
                        cbv_const.KEY_POS   : snp_info[dev_const.KEY_SCORES_POS],
                        cbv_const.KEY_REF   : snp_info[dev_const.KEY_SCORES_REF],
                        cbv_const.KEY_ALT   : snp_info[dev_const.KEY_SCORES_ALT],
                        }
            prediction = {cbv_const.KEY_TARGETS : rec[cbv_const.KEY_PREDICTION_SECTION][dev_const.KEY_SCORES_TARGETS]}
            self.dataset.append({cbv_const.KEY_SNP_INFO_SECTION : snp_data,
                                 cbv_const.KEY_PREDICTION_SECTION : prediction,
                                 cbv_const.KEY_SCORES_SECTION : rec[cbv_const.KEY_SCORES_SECTION],
                                 })



def filter_cbv_data(cbv_file,
                    cfg_file=cbv_const.CBV_CFG_FILE):
    (dir_name, file_name) = os.path.split(cbv_file)
    print
    print "> > >  " + file_name
    dm = DataSetManager()
    dm.load_data(cbv_file, file_type=cbv_const.FILE_TYPE_CBV)
    print "%-25s: %5d\n" % ("Original", len(dm.dataset))

    dm.validate_data()
    f_clean = open(cbv_file + '.clean', 'w')
    for item in dm.dataset:
        f_clean.write("%s\t%s\t%s\t%s\t%s\n" % (item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CHROM],
                                                 item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_POS],
                                                 item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_REF],
                                                 item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_ALT],
                                                 item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS],
                                                 )
                        )
    f_clean.close()
    print "%-25s: %5d" % ("Clean pathogenic", len([item for item in dm.dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Clean neutral", len([item for item in dm.dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '0']))
    print "%-25s: %5d\n" % ("Total", len(dm.dataset))

    dm.calculate_scores()
    f_scores = open(cbv_file + '.scores', 'w')
    for item in dm.dataset:
        f_scores.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CHROM],
                                                                         item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_POS],
                                                                         item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_REF],
                                                                         item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_ALT],
                                                                         item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS],
                                                                         item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_PHYLOP_SCORE],
                                                                         item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_SIFT_SCORE],
                                                                         item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_PP2_SCORE],
                                                                         item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_LRT_SCORE],
                                                                         item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_MT_SCORE],
                                                                         item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_GERP_SCORE],
                                                                         )
                        )
    f_scores.close()
    print "%-25s: %5d" % ("Scored pathogenic", len([item for item in dm.dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Scored neutral", len([item for item in dm.dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '0']))
    print "%-25s: %5d\n" % ("Total", len(dm.dataset))

    dm.set_shuffle_seed(cbv_const.DEMO_SEED)
    dm.shuffle_data()
    dm.partition_data()

    #partition data
    training_dataset      = dm.get_training_data() 
    validation_dataset    = dm.get_validation_data() 

    print "%-25s: %5d" % ("Training pathogenic", len([item for item in training_dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Training neutral", len([item for item in training_dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '0']))
    print "%-25s: %5d\n" % ("Total", len(training_dataset))

    print "%-25s: %5d" % ("Validation pathogenic", len([item for item in validation_dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '1']))
    print "%-25s: %5d" % ("Validation neutral", len([item for item in validation_dataset if item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] == '0']))
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
                  params_out_file=cbv_const.USER_PARAMS_FILE,
                  random_seed=cbv_const.DEFAULT_SEED,
                  n_hidden_nodes=cbv_const.DEFAULT_HIDDEN_NODES,
                  figure_dir=cbv_const.DEFAULT_FIGURE_DIR,
                  iterations=cbv_const.DEFAULT_ITERATIONS,
                  cfg_file=cbv_const.CBV_CFG_FILE,
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
    dm = FastDataSetManager(cfg_file=cfg_file)
    dm.load_data(training_data_file, file_type=dev_const.FILE_TYPE_SCORES)
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
    if not os.path.exists(cbv_const.USER_PARAMS_DIR):
        os.makedirs(cbv_const.USER_PARAMS_DIR)
    trainer.export_best_parameters(params_out_file)

def fast_predict(SNPs_file,
                 params_file=cbv_const.USER_PARAMS_FILE,
                 file_type=cbv_const.FILE_TYPE_VCF,
                 output_file=None,
                 cfg_file=cbv_const.CBV_CFG_FILE,
                 ):
    """

    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
    CBV has 5 fields, CHROM, POS, REF, ALT, EFFECT (1=deleterious, 0=neutral). All are tab separated
    Required arguments
    - SNPs_file : list of SNPs to be predicted, can be either VCF or CBV (default is VCF)

    """
    #pre-processing test dataset
    print >> sys.stderr, 'pre-processing dataset, this may take a while (around 750 SNPs/mins). . . . '
    dm = FastDataSetManager(cfg_file=cfg_file)
    dm.load_data(SNPs_file, file_type=dev_const.FILE_TYPE_SCORES)

    #predict
    predictor = Predictor()
    predictor.import_parameters(params_file=params_file)
    out = (np.array(predictor.predict(dm.dataset)).reshape(-1,))

    #print output
    if output_file is not None:
        sys.stdout = open(output_file, 'w')
    tmp_rec = []
    tmp_rec.append("CHROM")
    tmp_rec.append("POS")
    tmp_rec.append("REF")
    tmp_rec.append("ALT")
    tmp_rec.append("ACTUAL_DELETERIOUS_EFFECT")
    tmp_rec.append("PREDICTED_DELETERIOUS_PROBABILITY")
    tmp_rec.append("PHYLOP_SCORE")
    tmp_rec.append("SIFT_SCORE")
    tmp_rec.append("PP2_SCORE")
    tmp_rec.append("LRT_SCORT")
    tmp_rec.append("MT_SCORE")
    tmp_rec.append("GERP_SCORE")
    print "#" + "\t".join(tmp_rec)
    for i in xrange(len(dm.dataset)):
        del tmp_rec[:]
        snp_info   = dm.dataset[i][cbv_const.KEY_SNP_INFO_SECTION]
        prediction = dm.dataset[i][cbv_const.KEY_PREDICTION_SECTION]
        scores     = dm.dataset[i][cbv_const.KEY_SCORES_SECTION]
        tmp_rec.append(snp_info[cbv_const.KEY_CHROM])
        tmp_rec.append(snp_info[cbv_const.KEY_POS])
        tmp_rec.append(snp_info[cbv_const.KEY_REF])
        tmp_rec.append(snp_info[cbv_const.KEY_ALT])
        tmp_rec.append(prediction[cbv_const.KEY_TARGETS])
        tmp_rec.append("%6.4f" % out[i])
        tmp_rec.append(scores[cbv_const.KEY_PHYLOP_SCORE])
        tmp_rec.append(scores[cbv_const.KEY_SIFT_SCORE])
        tmp_rec.append(scores[cbv_const.KEY_PP2_SCORE])
        tmp_rec.append(scores[cbv_const.KEY_LRT_SCORE])
        tmp_rec.append(scores[cbv_const.KEY_MT_SCORE])
        tmp_rec.append(scores[cbv_const.KEY_GERP_SCORE])
        print "\t".join(tmp_rec)
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
