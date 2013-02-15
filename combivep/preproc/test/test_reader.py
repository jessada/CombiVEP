import unittest
import os
import combivep.settings as cbv_const
import combivep.preproc.reader as combivep_reader
from combivep.preproc.test.template import SafePreProcTester


class TestUcscReader(SafePreProcTester):


    def __init__(self, test_name):
        SafePreProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'ucsc_reader'

    def __create_ucsc_reader_instance(self):
        reader = combivep_reader.UcscReader()
        return reader

    def test_fetch_snps1(self):
        self.init_test('test_fetch_snps1')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('chr3', 110030150, 110030300))
        self.assertEqual(len(list(records)),
                         3,
                         "Incorrect number of records are being fetched")

    def test_fetch_snps2(self):
        self.init_test('test_fetch_snps2')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('3', 110030150, 110030300))
        self.assertEqual(len(records),
                         3,
                         "Incorrect number of records are being fetched")

    def test_formatting(self):
        self.init_test('test_formatting')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        readable = False
        for rec in reader.fetch_hash_snps('3', 110030157, 110030158):
            readable = True
            snp_info = rec[cbv_const.KEY_SNP_INFO_SECTION]
            self.assertEqual(snp_info[cbv_const.KEY_UCSC_CHROM],
                             'chr3',
                             "Incorrect UCSC formatting")
            self.assertEqual(snp_info[cbv_const.KEY_UCSC_START_POS],
                             '110030157',
                             "Incorrect UCSC formatting")
            self.assertEqual(snp_info[cbv_const.KEY_UCSC_END_POS],
                             '110030158',
                             "Incorrect UCSC formatting")
            self.assertEqual(snp_info[cbv_const.KEY_UCSC_STRAND],
                             '+',
                             "Incorrect UCSC formatting")
            self.assertEqual(snp_info[cbv_const.KEY_UCSC_REF],
                             'C',
                             "Incorrect UCSC formatting")
            self.assertEqual(snp_info[cbv_const.KEY_UCSC_OBSERVED],
                             'C/T',
                             "Incorrect UCSC formatting")
            break
        self.assertTrue(readable,
                        "UCSC database is not readable")

    def test_valid_indexing1(self):
        self.init_test('test_valid_indexing1')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('3', 110030330, 110030332))
        self.assertEqual(len(records),
                         2,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing2(self):
        self.init_test('test_valid_indexing2')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('3', 110030330, 110030331))
        self.assertEqual(len(records),
                         1,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing3(self):
        self.init_test('test_valid_indexing3')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('3', 110030331, 110030332))
        self.assertEqual(len(records),
                         1,
                         "Incorrect number of records are being fetched")

    def tearDown(self):
        self.remove_working_dir()


class TestLjbReader(SafePreProcTester):


    def __init__(self, test_name):
        SafePreProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'ljb_reader'

    def __create_ljb_reader_instance(self):
        reader = combivep_reader.LjbReader()
        return reader

    def test_valid_indexing1(self):
        self.init_test('test_valid_indexing1')
        reader = self.__create_ljb_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ljb_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('3', 108541778, 108541779))
        self.assertEqual(len(records),
                         3,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing2(self):
        self.init_test('test_valid_indexing1')
        reader = self.__create_ljb_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ljb_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('3', 108541778, 108541778))
        self.assertEqual(len(records),
                         2,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing3(self):
        self.init_test('test_valid_indexing1')
        reader = self.__create_ljb_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ljb_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_hash_snps('3', 108541779, 108541779))
        self.assertEqual(len(records),
                         1,
                         "Incorrect number of records are being fetched")

    def test_formatting(self):
        self.init_test('test_formatting')
        reader = self.__create_ljb_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ljb_reader.txt.gz')
        reader.read(test_file)
        readable = False
        for rec in reader.fetch_hash_snps('3', 108541778, 108541779):
            readable = True
            snp_info = rec[cbv_const.KEY_SNP_INFO_SECTION]
            self.assertEqual(snp_info[cbv_const.KEY_LJB_CHROM],
                             '3',
                             "Incorrect LJB formatting")
            self.assertEqual(snp_info[cbv_const.KEY_LJB_POS],
                             '108541778',
                             "Incorrect LJB formatting")
            self.assertEqual(snp_info[cbv_const.KEY_LJB_REF],
                             'T',
                             "Incorrect LJB formatting")
            self.assertEqual(snp_info[cbv_const.KEY_LJB_ALT],
                             'C',
                             "Incorrect LJB formatting")
            scores = rec[cbv_const.KEY_SCORES_SECTION]
            self.assertEqual(scores[cbv_const.KEY_PHYLOP_SCORE],
                             '0.102322',
                             "Incorrect LJB formatting")
            self.assertEqual(scores[cbv_const.KEY_SIFT_SCORE],
                             '0.91',
                             "Incorrect LJB formatting")
            self.assertEqual(scores[cbv_const.KEY_PP2_SCORE],
                             '0',
                             "Incorrect LJB formatting")
            self.assertEqual(scores[cbv_const.KEY_LRT_SCORE],
                             '0.312516',
                             "Incorrect LJB formatting")
            self.assertEqual(scores[cbv_const.KEY_MT_SCORE],
                             '0.000000',
                             "Incorrect LJB formatting")
            self.assertEqual(scores[cbv_const.KEY_GERP_SCORE],
                             '-3.16',
                             "Incorrect LJB formatting")
            break
        self.assertTrue(readable,
                        "LJB database is not readable")

    def tearDown(self):
        self.remove_working_dir()


class TestVcfReader(SafePreProcTester):


    def __init__(self, test_name):
        SafePreProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'vcf_reader'

    def __create_vcf_reader_instance(self):
        reader = combivep_reader.VcfReader()
        return reader

    def test_formatting(self):
        self.init_test('test_formatting')
        reader = self.__create_vcf_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_formatting.vcf')
        reader.read(test_file)
        readable = False
        for rec in reader.fetch_hash_snps():
            readable = True
            snp_info = rec[cbv_const.KEY_SNP_INFO_SECTION]
            self.assertEqual(snp_info[cbv_const.KEY_VCF_CHROM],
                             '1',
                             "Incorrect VCF formatting")
            self.assertEqual(snp_info[cbv_const.KEY_VCF_POS],
                             '887560',
                             "Incorrect VCF formatting")
            self.assertEqual(snp_info[cbv_const.KEY_VCF_REF],
                             'A',
                             "Incorrect VCF formatting")
            self.assertEqual(snp_info[cbv_const.KEY_VCF_ALT],
                             'C',
                             "Incorrect VCF formatting")
            break
        self.assertTrue(readable,
                        "VCF reader does not work properly")

    def tearDown(self):
        self.remove_working_dir()


class TestCbvReader(SafePreProcTester):


    def __init__(self, test_name):
        SafePreProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'cbv_reader'

    def __create_cbv_reader_instance(self):
        reader = combivep_reader.CbvReader()
        return reader

    def test_formatting(self):
        self.init_test('test_formatting')
        reader = self.__create_cbv_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_formatting.cbv')
        reader.read(test_file)
        readable = False
        for rec in reader.fetch_hash_snps():
            readable = True
            snp_info = rec[cbv_const.KEY_SNP_INFO_SECTION]
            self.assertEqual(snp_info[cbv_const.KEY_CBV_CHROM],
                             '1',
                             "Incorrect CBV formatting")
            self.assertEqual(snp_info[cbv_const.KEY_CBV_POS],
                             '19566382',
                             "Incorrect CBV formatting")
            self.assertEqual(snp_info[cbv_const.KEY_CBV_REF],
                             'T',
                             "Incorrect CBV formatting")
            self.assertEqual(snp_info[cbv_const.KEY_CBV_ALT],
                             'C',
                             "Incorrect CBV formatting")
            prediction = rec[cbv_const.KEY_PREDICTION_SECTION]
            self.assertEqual(prediction[cbv_const.KEY_CBV_TARGETS],
                             '1',
                             "Incorrect CBV formatting")
            break
        self.assertTrue(readable,
                        "CBV reader does not work properly")

    def test_counting(self):
        self.init_test('test_formatting')
        reader = self.__create_cbv_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_formatting.cbv')
        reader.read(test_file)
        self.assertEqual(len(list(reader.fetch_hash_snps())),
                         11,
                         'CbvReader does not read input file properly')

    def tearDown(self):
        self.remove_working_dir()




