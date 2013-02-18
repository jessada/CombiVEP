import unittest
import os
import filecmp
from combivep.refdb.test.template import SafeRefDBTester
import combivep.settings as cbv_const
from combivep.refdb.control import UcscController
from combivep.refdb.control import LjbController
from combivep.preproc.reader import UcscReader
from combivep.preproc.reader import LjbReader


class TestUcscController(SafeRefDBTester):

    def __init__(self, test_name):
        SafeRefDBTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'ucsc_controller'

    def __create_ucsc_ctrller_instance(self):
        ucsc_ctrller = UcscController()
        return ucsc_ctrller

    def test_tabix_db(self):
        #init
        self.init_test('test_tabix_database')
        ucsc_ctrller = self.__create_ucsc_ctrller_instance()
        test_file         = os.path.join(self.data_dir,
                                         'test_tabix_database.txt')
        working_file      = os.path.join(self.working_dir,
                                         'test_tabix_database.txt')
        expected_out_file = os.path.join(self.working_dir,
                                         'test_tabix_database.txt.gz')
        self.copy_file(test_file, working_file)

        #test if the 'tabix' files are produced
        out_file = ucsc_ctrller.tabix_db(working_file)
        self.assertEqual(expected_out_file, out_file,
                         "Tabix doesn't work correctly")
        self.assertTrue(os.path.exists(out_file),
                        "Tabix doesn't work correctly")
        self.assertTrue(os.path.exists(out_file+'.tbi'),
                        "Tabix doesn't work correctly")

        #test if it is readable
        ucsc_reader = UcscReader()
        ucsc_reader.read(out_file)
        readable = False
        for rec in ucsc_reader.fetch_snps('chr3', 108572604, 108572605):
            readable = True
            self.assertEqual(rec.start_pos,
                             '108572604', 
                             "Database tabixing doesn't work correctly")
            break
        self.assertTrue(readable, 
                        "Tabixed ucsc database is not readable")

    @unittest.skip("temporary disable due to high bandwidth usage")
    def test_not_update(self):
        #init
        self.init_test('test_not_update')
        ucsc_ctrller = self.__create_ucsc_ctrller_instance()
        ucsc_ctrller.cfg_file = os.path.join(self.data_dir,
                                             'test_not_update_cfg_file.txt')
        #run test
        self.assertFalse(ucsc_ctrller.update(),
                         "UCSC controller can't identify valid update status")

    def tearDown(self):
        self.remove_working_dir()


class TestLjbController(SafeRefDBTester):

    def __init__(self, test_name):
        SafeRefDBTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'ljb_controller'

    def __create_ljb_ctrller_instance(self):
        ljb_ctrller = LjbController()
        return ljb_ctrller

    def test_clean_raw_db(self):
        #initialize variables
#        self.individual_debug = True
        self.init_test('test_clean_raw_database')
        ljb_ctrller = self.__create_ljb_ctrller_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_clean_raw_database.txt')
        out_file = os.path.join(self.working_dir,
                                'clean_database.txt')
        expected_out_file = os.path.join(self.data_dir,
                                         'expected_clean_database.txt')

        #call function
        ljb_ctrller.clean_raw_db(test_file, out_file)
        self.assertTrue(filecmp.cmp(out_file, expected_out_file),
                        "Raw LJB database haven't been clean properly")

    def test_tabix_db(self):
        #init
        self.init_test('test_tabix_database')
        ljb_ctrller = self.__create_ljb_ctrller_instance()
        test_file    = os.path.join(self.data_dir,
                                    'test_tabix_database.txt')
        working_file = os.path.join(self.working_dir,
                                    'test_tabix_database.txt')
        out_file     = os.path.join(self.working_dir,
                                    'test_tabix_database.txt.gz')
        self.copy_file(test_file, working_file)

        #test if the 'tabix' files are produced
        ljb_ctrller.tabix_db(working_file)
        self.assertTrue(os.path.exists(out_file),
                        "Tabix doesn't work correctly")
        self.assertTrue(os.path.exists(out_file+'.tbi'),
                        "Tabix doesn't work correctly")

        #test if it is readable
        ljb_reader = LjbReader()
        ljb_reader.read(out_file)
        readable = False
        for rec in ljb_reader.fetch_snps('3', 108549516, 108549517):
            readable = True
            self.assertEqual(rec.pos,
                             '108549517',
                             "Database tabixing doesn't work correctly")
            break
        self.assertTrue(readable,
                        "Tabixed ljb database is not readable")

    def test_concat_chromosone_files(self):
        #init
        self.individual_debug = True
        self.init_test('test_concat_chromosome_files')
        ljb_ctrller = self.__create_ljb_ctrller_instance()
        file_prefix = os.path.join(self.data_dir,
                                   'dummy_dbNSFP_light1.3')
        file_suffix = '.txt'
        out_file = os.path.join(self.working_dir,
                                'out_concat.txt')
        expected_out_file = os.path.join(self.data_dir,
                                         'expected_concat_file.txt')

        #runtest
        ljb_ctrller.concat_chromosome_files(file_prefix,
                                            file_suffix,
                                            out_file)
        self.assertTrue(filecmp.cmp(out_file, expected_out_file),
                                    "error in chromosome files concatenation")

    def test_not_update(self):
        #init
        self.init_test('test_not_update')
        ljb_ctrller = self.__create_ljb_ctrller_instance()
        ljb_ctrller.cfg_file = os.path.join(self.data_dir,
                                            'test_not_update_cfg_file.txt')
        #run test
        self.assertFalse(ljb_ctrller.update(),
                         "LJB controller can't identify valid update status")

    def tearDown(self):
        self.remove_working_dir()
