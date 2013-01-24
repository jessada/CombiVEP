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
#def ucsc_join(cbv_file):
#    cfg = Configure()
#    cfg.load_config()
#    fh = open(cbv_file);
#    tabix_file = pysam.Tabixfile(cfg.config_values[combivep_settings.LATEST_UCSC_FILE_NAME])
#    for line in fh:
#        cbv_data = line.strip().split('\t')
#        for line in tabix_file.fetch('chr' + cbv_data[combivep_settings.CBV_0_INDEX_CHROM],
#                                     int(cbv_data[combivep_settings.CBV_0_INDEX_POS]) - 2, 
#                                     int(cbv_data[combivep_settings.CBV_0_INDEX_POS]) + 2,
#                                     ):
##            print cbv_data, line.strip()
#            break
#    fh.close()

#        if chromosome.startswith('chr'):
#            chrom = chromosome
#        else:
#            chrom = 'chr' + chromosome
#            yield line.rstrip('\n').split('\t')
#
#    def fetch_hash_snps(self, chromosome, start_pos, end_pos):
#        for rec in self.fetch_array_snps(chromosome, start_pos, end_pos):
#            if len(rec) != combivep_settings.UCSC_EXPECTED_LENGTH :
#                raise Exception("Invalid formatting is found in file '%s'>> Chrom : %s\tStart pos : %s\tEnd pos : %s" % (self.db_file_name, rec[combivep_settings.UCSC_0_INDEX_CHROM], rec[combivep_settings.UCSC_0_INDEX_START_POS], rec[combivep_settings.UCSC_0_INDEX_END_POS]))
#            snp_info = {combivep_settings.KEY_UCSC_CHROM     : rec[combivep_settings.UCSC_0_INDEX_CHROM], 
#                        combivep_settings.KEY_UCSC_START_POS : rec[combivep_settings.UCSC_0_INDEX_START_POS],
#                        combivep_settings.KEY_UCSC_END_POS   : rec[combivep_settings.UCSC_0_INDEX_END_POS],
#                        combivep_settings.KEY_UCSC_STRAND    : rec[combivep_settings.UCSC_0_INDEX_STRAND],
#                        combivep_settings.KEY_UCSC_REF       : rec[combivep_settings.UCSC_0_INDEX_REF],
#                        combivep_settings.KEY_UCSC_OBSERVED  : rec[combivep_settings.UCSC_0_INDEX_OBSERVED],
#                        }
#            yield {combivep_settings.KEY_SNP_INFO_SECTION : snp_info}






