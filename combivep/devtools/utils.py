import os
import sys
import numpy as np
from combivep.preproc.dataset import DataSetManager
import combivep.settings as combivep_settings
from combivep.cfg import Configure


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
        f_scores.write("%s\t%s\t%s\t%s\t%s\n" % (item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_CHROM],
                                                 item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_POS],
                                                 item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_REF],
                                                 item[combivep_settings.KEY_SNP_INFO_SECTION][combivep_settings.KEY_ALT],
                                                 item[combivep_settings.KEY_PREDICTION_SECTION][combivep_settings.KEY_TARGETS],
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






