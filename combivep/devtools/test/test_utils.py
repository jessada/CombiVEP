import os
import unittest
import combivep.settings as cbv_const
from combivep.devtools.test.template import SafeMiscTester
from combivep.devtools.utils import filter_cbv_data

class Test_Gadget(SafeMiscTester):


    def __init__(self, test_name):
        SafeMiscTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'utils'

    def init_configure_instance(self):
        pass

    def test_filter_cbv_data(self):
        #init
        self.individual_debug = True
        self.init_test('test_filter_cbv_data')
        test_file    = os.path.join(self.data_dir, 'test_filter_cbv_data.cbv')
        working_file = os.path.join(self.working_dir, 'test_filter_cbv_data.cbv')
        self.copy_file(test_file, working_file)
        filter_cbv_data(working_file, config_file=cbv_const.CBV_CENTRAL_TEST_CONFIG_FILE)

    def tearDown(self):
        self.remove_working_dir()


