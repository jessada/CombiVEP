import math
import random
import numpy as np
import combivep.settings as cbv_const
from combivep.template import CombiVEPBase
from combivep.preproc.reader import VcfReader
from combivep.preproc.reader import CbvReader
from combivep.preproc.referer import Referer


class SnpDataRecord(CombiVEPBase):
    """ to store precalculated scores """

    def __init__(self, rec):
        self.__rec = rec

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str({'chrom': self.chrom,
                    'pos': self.pos,
                    'ref': self.ref,
                    'alt': self.alt,
                    'target': self.target,
                    })

    @property
    def chrom(self):
        return self.__rec.chrom

    @property
    def pos(self):
        return self.__rec.pos

    @property
    def ref(self):
        return self.__rec.ref

    @property
    def alt(self):
        return self.__rec.alt

    @property
    def target(self):
        return self.__rec.target


class DataSet(list):
    """

    one dataset represents one set of data:
    - training
    - validation
    - test
    - positive
    - negative
    - etc..

    """

    def __init__(self, *args):
        list.__init__(self, *args)
        self.shuffle_seed = None

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
            scores = item[cbv_const.KW_SCORES]
            tmp_array = []
            tmp_array.append(float(scores.phylop_score))
            tmp_array.append(float(scores.sift_score))
            tmp_array.append(float(scores.pp2_score))
            tmp_array.append(float(scores.lrt_score))
            tmp_array.append(float(scores.mt_score))
            tmp_array.append(float(scores.gerp_score))
            feature_vector_arrays.append(tmp_array)
        return np.matrix(feature_vector_arrays).T

    @property
    def targets(self):
        return self.__get_targets()

    def __get_targets(self):
        if len(self) == 0:
            return None
        if self[0][cbv_const.KW_SNP_DATA].target is None:
            return None
        target_array = []
        for item in self:
            target_array.append(item[cbv_const.KW_SNP_DATA].target)
        return np.array(target_array).astype(np.int)

    @property
    def n_features(self):
        return self.__get_n_features()

    def __get_n_features(self):
        for item in self:
            return len(item[cbv_const.KW_SCORES])
            break

    @property
    def n_data(self):
        return self.__get_n_data()

    def __get_n_data(self):
        return len(self)


class DataSetManager(CombiVEPBase):

    def __init__(self, cfg_file=cbv_const.CBV_CFG_FILE):
        CombiVEPBase.__init__(self)

        self.referer = Referer()
        self.referer.cfg_file = cfg_file
        self.referer.load_cfg()
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
        for rec in vcf_reader.fetch_snps():
            snp_data = SnpDataRecord(rec)
            self.dataset.append({cbv_const.KW_SNP_DATA: snp_data})

    def __load_cbv_data(self, file_name):
        self.clear_data()
        cbv_reader = CbvReader()
        cbv_reader.read(file_name)
        for rec in cbv_reader.fetch_snps():
            snp_data = SnpDataRecord(rec)
            self.dataset.append({cbv_const.KW_SNP_DATA: snp_data})

    def validate_data(self):
        #to prevent misintepretion due to different version
        #between each data point by removing items from self.dataset
        #if they are not exist in certain UCSC database
        tmp_dataset = DataSet()
        for item in self.dataset:
            snp_data = item[cbv_const.KW_SNP_DATA]
            if self.referer.validate_snp(snp_data.chrom,
                                         snp_data.pos,
                                         snp_data.ref,
                                         snp_data.alt,
                                         ):
                tmp_dataset.append(item)
        del self.dataset[:]
        self.dataset = tmp_dataset

    def calculate_scores(self):
        #get scores from LJB database
        tmp_dataset = DataSet()
        for item in self.dataset:
            snp_data = item[cbv_const.KW_SNP_DATA]
            scores = self.referer.get_scores(snp_data.chrom,
                                             snp_data.pos,
                                             snp_data.ref,
                                             snp_data.alt,
                                             )
            item[cbv_const.KW_SCORES] = scores
        #remove items from self.dataset if they don't have scores
        self.dataset[:] = [item for item in self.dataset if item[cbv_const.KW_SCORES] is not None]

    def partition_data(self,
                       prop_training=cbv_const.PROPORTION_TRAINING_DATA,
                       prop_validation=cbv_const.PROPORTION_VALIDATION_DATA,
                       ):
        self.__prop_training   = prop_training
        self.__prop_validation = prop_validation

    @property
    def n_training_data(self):
        total_prop = float(self.__prop_training + self.__prop_validation)
        ratio_training = self.__prop_training / total_prop
        return int(math.floor(len(self.dataset) * ratio_training))

    @property
    def n_validation_data(self):
        return len(self.dataset) - self.n_training_data

    def get_training_data(self):
        dataset = DataSet()
        for i in xrange(0, self.n_training_data):
            dataset.append(self.dataset[i])
        return dataset

    def get_validation_data(self):
        dataset = DataSet()
        for i in xrange(self.n_training_data, len(self.dataset)):
            dataset.append(self.dataset[i])
        return dataset

    def set_shuffle_seed(self, shuffle_seed):
        self.dataset.set_shuffle_seed(shuffle_seed)

    def shuffle_data(self):
        self.dataset.shuffle()
