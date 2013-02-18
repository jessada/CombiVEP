import unittest
import os
import filecmp
import combivep.settings as cbv_const
from combivep.test.template import SafeGeneralTester
from combivep.config import Configure


class TestConfigure(SafeGeneralTester):

    def __init__(self, test_name):
        SafeGeneralTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'configure'

    def __create_configure_instance(self):
        cfg = Configure()
        return cfg

    def test_initial(self):
        #init
        self.init_test('test_initial')
        cfg = self.__create_configure_instance()
        out_file = os.path.join(self.working_dir,
                                'out_config.txt')
        cfg.cfg_file = out_file
        expected_out_file = os.path.join(self.data_dir,
                                         'expected_initial_config_out.txt')
        #run test
        cfg.write_ljb_cfg('4.5', '/home/jessada/development/scilifelab/projects/CombiVEP/combivep/refdb/test/data/ljb_controller/dummy_dbNSFP_light1.4')
        self.assertTrue(filecmp.cmp(out_file, expected_out_file),
                        "Configure cannot update initial config correctly")

    def test_load(self):
        #init
        self.init_test('test_load')
        cfg = self.__create_configure_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_load.txt')
        cfg.cfg_file = test_file
        #run test
        out = cfg.load_cfg()
        self.assertEqual(out[cbv_const.LATEST_UCSC_DB_VER],
                         '7.5',
                         "Configure cannot load configuration correctly")
        self.assertEqual(out[cbv_const.LATEST_UCSC_FILE_NAME],
                         'ucsc_file7.5.txt',
                         "Configure cannot load configuration correctly")
        self.assertEqual(out[cbv_const.LATEST_LJB_DB_VER],
                         '4.4',
                         "Configure cannot load configuration correctly")
        self.assertEqual(out[cbv_const.LATEST_LJB_FILE_PREFIX],
                         '/home/jessada/development/scilifelab/projects/CombiVEP/combivep/refdb/test/data/ljb_controller/dummy_dbNSFP_light1.3',
                         "Configure cannot load configuration correctly")

    def test_update_ljb_cfg(self):
        #init
        self.init_test('test_update_ljb_config')
        cfg = self.__create_configure_instance()
        test_file    = os.path.join(self.data_dir,
                                    'test_load.txt')
        out_file     = os.path.join(self.working_dir,
                                    'out_config.txt')
        cfg.cfg_file = out_file
        expected_out_file = os.path.join(self.data_dir,
                                         'expected_ljb_config_out.txt')
        self.copy_file(test_file, out_file)
        #run test
        cfg.write_ljb_cfg('4.5', '/home/jessada/development/scilifelab/projects/CombiVEP/combivep/refdb/test/data/ljb_controller/dummy_dbNSFP_light1.4')
        self.assertTrue(filecmp.cmp(out_file, expected_out_file),
                        "Configure cannot update LJB config correctly")

    def test_update_ucsc_cfg(self):
        #init
        self.init_test('test_update_ucsc_config')
        cfg = self.__create_configure_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_load.txt')
        out_file = os.path.join(self.working_dir,
                                'out_config.txt')
        cfg.cfg_file = out_file
        expected_out_file = os.path.join(self.data_dir,
                                         'expected_ucsc_config_out.txt')
        self.copy_file(test_file, out_file)
        #run test
        cfg.write_ucsc_cfg('7.6', 'ucsc_file7.6.txt')
        self.assertTrue(filecmp.cmp(out_file, expected_out_file),
                        "Configure cannot update UCSC config correctly")

    def tearDown(self):
        self.remove_working_dir()
