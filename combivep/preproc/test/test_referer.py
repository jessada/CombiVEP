import unittest
import os
import combivep.settings as cbv_const
from combivep.preproc.test.template import SafePreProcTester
from combivep.preproc.referer import Referer


class TestReferer(SafePreProcTester):


    def __init__(self, test_name):
        SafePreProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'referer'

    def __create_referer_instance(self):
        referer = Referer()
        return referer

    def test_validate_snp(self):
        self.init_test('test_validate_snp')
        referer = self.__create_referer_instance()
        referer.cfg_file = cbv_const.CBV_SAMPLE_CFG_FILE
        referer.load_cfg()
        self.assertTrue(referer.validate_snp('1', 887560, 'A', 'C'),
                        "Incorrect SNP validating")
        self.assertTrue(referer.validate_snp('chr3', 25836088, 'C', 'A'),
                        "Incorrect SNP validating")
        self.assertTrue(referer.validate_snp('20' , 17474690, 'T', 'G'),
                        "Incorrect SNP validating")
        self.assertTrue(referer.validate_snp('chrX', 56296488, 'G', 'C'),
                        "Incorrect SNP validating")
        self.assertTrue(referer.validate_snp('Y', 15581983, 'G', 'A'),
                        "Incorrect SNP validating")
        self.assertFalse(referer.validate_snp('chr16', 21086416, 'T', 'A'),
                         "Incorrect SNP validating")

    def test_get_scores(self):
        self.init_test('test_get_scores')
        referer = self.__create_referer_instance()
        referer.cfg_file = cbv_const.CBV_SAMPLE_CFG_FILE
        referer.load_cfg()
        rec = referer.get_scores('3', 108541778, 'T', 'C')
        self.assertEqual(rec.phylop_score,
                         '0.102322',
                         "Incorrect LJB formatting")
        self.assertEqual(rec.sift_score,
                         '0.91',
                         "Incorrect LJB formatting")
        self.assertEqual(rec.pp2_score,
                         '0',
                         "Incorrect LJB formatting")
        self.assertEqual(rec.lrt_score,
                         '0.312516',
                         "Incorrect LJB formatting")
        self.assertEqual(rec.mt_score,
                         '0.000000',
                         "Incorrect LJB formatting")
        self.assertEqual(rec.gerp_score,
                         '-3.16',
                         "Incorrect LJB formatting")

    def tearDown(self):
        self.remove_working_dir()






