import unittest
import os
import filecmp
from combivep.test.template import SafeGeneralTester
import combivep.settings as cbv_const
from combivep.config import Configure


class TestConfigure(SafeGeneralTester):


    def __init__(self, test_name):
        SafeGeneralTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'configure'

    def init_configure_instance(self):
        self.__configure = Configure()

    def test_initial(self):
        #init
        self.individual_debug = True
        self.init_test('test_initial')
        self.init_configure_instance()
        output_file = os.path.join(self.working_dir, 'out_config.txt')
        self.__configure.config_file = output_file
        expected_out_file = os.path.join(self.data_dir, 'expected_update_initial_config_output.txt')
        #run test
        self.__configure.write_ljb_config('4.5', '/home/jessada/development/scilifelab/projects/CombiVEP/combivep/refdb/test/data/ljb_controller/dummy_dbNSFP_light1.4')
        self.assertTrue(filecmp.cmp(output_file, expected_out_file), "Configure cannot update LJB config correctly")

    def test_load(self):
        #init
        self.init_test('test_load')
        self.init_configure_instance()
        test_file = os.path.join(self.data_dir, 'test_load.txt')
        self.__configure.config_file = test_file
        #run test
        out = self.__configure.load_config()
        self.assertEqual(out[cbv_const.LATEST_UCSC_DB_VERSION], '7.5', "Configure cannot load configuration correctly")
        self.assertEqual(out[cbv_const.LATEST_UCSC_FILE_NAME], 'ucsc_file7.5.txt', "Configure cannot load configuration correctly")
        self.assertEqual(out[cbv_const.LATEST_LJB_DB_VERSION], '4.4', "Configure cannot load configuration correctly")
        self.assertEqual(out[cbv_const.LATEST_LJB_FILE_PREFIX], '/home/jessada/development/scilifelab/projects/CombiVEP/combivep/refdb/test/data/ljb_controller/dummy_dbNSFP_light1.3', "Configure cannot load configuration correctly")

    def test_update_ljb_config(self):
        #init
        self.init_test('test_update_ljb_config')
        self.init_configure_instance()
        test_file = os.path.join(self.data_dir, 'test_load.txt')
        output_file = os.path.join(self.working_dir, 'out_config.txt')
        self.copy_file(test_file, output_file)
        self.__configure.config_file = output_file
        expected_out_file = os.path.join(self.data_dir, 'expected_update_ljb_config_output.txt')
        #run test
        self.__configure.write_ljb_config('4.5', '/home/jessada/development/scilifelab/projects/CombiVEP/combivep/refdb/test/data/ljb_controller/dummy_dbNSFP_light1.4')
        self.assertTrue(filecmp.cmp(output_file, expected_out_file), "Configure cannot update LJB config correctly")

    def test_update_ucsc_config(self):
        #init
        self.init_test('test_update_ucsc_config')
        self.init_configure_instance()
        test_file = os.path.join(self.data_dir, 'test_load.txt')
        output_file = os.path.join(self.working_dir, 'out_config.txt')
        self.copy_file(test_file, output_file)
        self.__configure.config_file = output_file
        expected_out_file = os.path.join(self.data_dir, 'expected_update_ucsc_config_output.txt')
        #run test
        self.__configure.write_ucsc_config('7.6', 'ucsc_file7.6.txt')
        self.assertTrue(filecmp.cmp(output_file, expected_out_file), "Configure cannot update UCSC config correctly")

    def tearDown(self):
        self.remove_working_dir()


