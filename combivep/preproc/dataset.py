import math
import random
import numpy as np
from combivep.template import CombiVEPBase
import combivep.settings as cbv_const
from combivep.preproc.reader import VcfReader
from combivep.preproc.reader import CbvReader
from combivep.preproc.referer import Referer


class DataSet(list):


    def __init__(self, *args):
        list.__init__(self, *args)
        self.shuffle_seed     = None

    def clear(self):
        del self[:]

    def shuffle(self):
        random.seed(self.shuffle_seed)
        random.shuffle(self)

    def set_shuffle_seed(self, shuffle_seed):
        self.shuffle_seed = shuffle_seed

    def __add__(self, other):
        for item in other:
            super(DataSet, self).append(item)
        return self

    @property
    def feature_vectors(self):
        return self.__get_feature_vectors()

    def __get_feature_vectors(self):
        feature_vector_arrays = []
        for item in self:
            tmp_array = []
            tmp_array.append(float(item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_PHYLOP_SCORE]))
            tmp_array.append(float(item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_SIFT_SCORE]))
            tmp_array.append(float(item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_PP2_SCORE]))
            tmp_array.append(float(item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_LRT_SCORE]))
            tmp_array.append(float(item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_MT_SCORE]))
            tmp_array.append(float(item[cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_GERP_SCORE]))
            feature_vector_arrays.append(tmp_array)
        return np.matrix(feature_vector_arrays).T

    @property
    def targets(self):
        return self.__get_targets()

    def __get_targets(self):
        if len(self) == 0:
            return None
        if self[0][cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS] is None:
            return None
        target_array = []
        for item in self:
            target_array.append(item[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS])
        return np.array(target_array).astype(np.int)

    @property
    def n_features(self):
        return self.__get_n_features()

    def __get_n_features(self):
        for item in self:
            return len(item[cbv_const.KEY_SCORES_SECTION].keys())
            break

    @property
    def n_data(self):
        return self.__get_n_data()

    def __get_n_data(self):
        return len(self)


class DataSetManager(CombiVEPBase):


    def __init__(self, config_file=cbv_const.CBV_CONFIG_FILE):
        CombiVEPBase.__init__(self)

        self.referer = Referer()
        self.referer.config_file = config_file
        self.referer.load_config()
        self.dataset = DataSet()

    def clear_data(self):
        self.dataset.clear()

    def load_data(self, file_name, file_type=cbv_const.FILE_TYPE_VCF):
        if file_type == cbv_const.FILE_TYPE_VCF:
            return self.__load_vcf_data(file_name)
        if file_type == cbv_const.FILE_TYPE_CBV:
            return self.__load_cbv_data(file_name)

    def __load_vcf_data(self, file_name):
        self.clear_data()
        vcf_reader = VcfReader()
        vcf_reader.read(file_name)
        for rec in vcf_reader.fetch_hash_snps():
            snp_data = {cbv_const.KEY_CHROM : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_VCF_CHROM],
                        cbv_const.KEY_POS   : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_VCF_POS],
                        cbv_const.KEY_REF   : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_VCF_REF],
                        cbv_const.KEY_ALT   : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_VCF_ALT],
                        }
            prediction = {cbv_const.KEY_TARGETS : None}
            self.dataset.append({cbv_const.KEY_SNP_INFO_SECTION   : snp_data,
                                 cbv_const.KEY_PREDICTION_SECTION : prediction})

    def __load_cbv_data(self, file_name):
        self.clear_data()
        cbv_reader = CbvReader()
        cbv_reader.read(file_name)
        for rec in cbv_reader.fetch_hash_snps():
            snp_data = {cbv_const.KEY_CHROM : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CBV_CHROM],
                        cbv_const.KEY_POS   : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CBV_POS],
                        cbv_const.KEY_REF   : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CBV_REF],
                        cbv_const.KEY_ALT   : rec[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CBV_ALT],
                        }
            prediction = {cbv_const.KEY_TARGETS : rec[cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_CBV_TARGETS]}
            self.dataset.append({cbv_const.KEY_SNP_INFO_SECTION : snp_data,
                                 cbv_const.KEY_PREDICTION_SECTION : prediction})

    def validate_data(self):
        #to prevent misintepret due to different version between each data point by 
        #removing items from self.dataset if they are not exist in certain UCSC database
        self.dataset[:] = [item for item in self.dataset if self.referer.validate_snp(item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CHROM],
                                                                                      item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_POS],
                                                                                      item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_REF],
                                                                                      item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_ALT]
                                                                                      )]

    def calculate_scores(self):
        #get scores from LJB database
        for item in self.dataset:
            item[cbv_const.KEY_SCORES_SECTION] = self.referer.get_scores(item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CHROM],
                                                                                 item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_POS],
                                                                                 item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_REF],
                                                                                 item[cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_ALT]
                                                                                 )
        #remove items from self.dataset if they don't have scores
        self.dataset[:] = [item for item in self.dataset if item[cbv_const.KEY_SCORES_SECTION] is not None]

    def partition_data(self,
                       proportion_training_data   = cbv_const.PROPORTION_TRAINING_DATA,
                       proportion_validation_data = cbv_const.PROPORTION_VALIDATION_DATA,
                       ):
        total_proportion = proportion_training_data + proportion_validation_data
        self.training_data_size   = int(math.floor(len(self.dataset) * proportion_training_data / total_proportion))
        self.validation_data_size = len(self.dataset) - self.training_data_size

    def get_training_data(self):
        dataset = DataSet()
        for i in  xrange(0, self.training_data_size):
            dataset.append(self.dataset[i])
        return dataset

    def get_validation_data(self):
        dataset = DataSet()
        for i in xrange(self.training_data_size, len(self.dataset)):
            dataset.append(self.dataset[i])
        return dataset

    def set_shuffle_seed(self, shuffle_seed):
        self.dataset.set_shuffle_seed(shuffle_seed)

    def shuffle_data(self):
        self.dataset.shuffle()



