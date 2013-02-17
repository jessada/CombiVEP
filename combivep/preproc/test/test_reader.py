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
        records = list(reader.fetch_snps('chr3', 110030150, 110030300))
        self.assertEqual(len(list(records)),
                         3,
                         "Incorrect number of records are being fetched")

    def test_fetch_snps2(self):
        self.init_test('test_fetch_snps2')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_snps('3', 110030150, 110030300))
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
        for rec in reader.fetch_snps('3', 110030157, 110030158):
            readable = True
            self.assertEqual(rec.chrom,
                             'chr3',
                             "Incorrect UCSC formatting")
            self.assertEqual(rec.start_pos,
                             '110030157',
                             "Incorrect UCSC formatting")
            self.assertEqual(rec.end_pos,
                             '110030158',
                             "Incorrect UCSC formatting")
            self.assertEqual(rec.strand,
                             '+',
                             "Incorrect UCSC formatting")
            self.assertEqual(rec.ref,
                             'C',
                             "Incorrect UCSC formatting")
            self.assertEqual(rec.observed,
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
        records = list(reader.fetch_snps('3', 110030330, 110030332))
        self.assertEqual(len(records),
                         2,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing2(self):
        self.init_test('test_valid_indexing2')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_snps('3', 110030330, 110030331))
        self.assertEqual(len(records),
                         1,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing3(self):
        self.init_test('test_valid_indexing3')
        reader = self.__create_ucsc_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ucsc_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_snps('3', 110030331, 110030332))
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
        records = list(reader.fetch_snps('3', 108541778, 108541779))
        self.assertEqual(len(records),
                         3,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing2(self):
        self.init_test('test_valid_indexing1')
        reader = self.__create_ljb_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ljb_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_snps('3', 108541778, 108541778))
        self.assertEqual(len(records),
                         2,
                         "Incorrect number of records are being fetched")

    def test_valid_indexing3(self):
        self.init_test('test_valid_indexing1')
        reader = self.__create_ljb_reader_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_ljb_reader.txt.gz')
        reader.read(test_file)
        records = list(reader.fetch_snps('3', 108541779, 108541779))
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
        for rec in reader.fetch_snps('3', 108541778, 108541779):
            readable = True
            self.assertEqual(rec.chrom,
                             '3',
                             "Incorrect LJB formatting")
            self.assertEqual(rec.pos,
                             '108541778',
                             "Incorrect LJB formatting")
            self.assertEqual(rec.ref,
                             'T',
                             "Incorrect LJB formatting")
            self.assertEqual(rec.alt,
                             'C',
                             "Incorrect LJB formatting")
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
        for rec in reader.fetch_snps():
            readable = True
            self.assertEqual(rec.chrom,
                             '1',
                             "Incorrect VCF formatting")
            self.assertEqual(rec.pos,
                             '887560',
                             "Incorrect VCF formatting")
            self.assertEqual(rec.ref,
                             'A',
                             "Incorrect VCF formatting")
            self.assertEqual(rec.alt,
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
        for rec in reader.fetch_snps():
            readable = True
            self.assertEqual(rec.chrom,
                             '1',
                             "Incorrect CBV formatting")
            self.assertEqual(rec.pos,
                             '19566382',
                             "Incorrect CBV formatting")
            self.assertEqual(rec.ref,
                             'T',
                             "Incorrect CBV formatting")
            self.assertEqual(rec.alt,
                             'C',
                             "Incorrect CBV formatting")
            self.assertEqual(rec.targets,
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
        self.assertEqual(len(list(reader.fetch_snps())),
                         11,
                         'CbvReader does not read input file properly')

    def tearDown(self):
        self.remove_working_dir()




