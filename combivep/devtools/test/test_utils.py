import os
import unittest
from combivep.devtools.test.template import SafeMiscTester
import combivep.settings as combivep_settings
from combivep.devtools.utils import filter_cbv_data
#from combivep.devtools.utils import ucsc_join

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
        filter_cbv_data(working_file, config_file=combivep_settings.COMBIVEP_CENTRAL_TEST_CONFIGURATION_FILE)

#    def test_ucsc_join(self):
#        #init
#        self.individual_debug = True
#        self.init_test('test_ucsc_join')
#        test_file    = os.path.join(self.data_dir, 'test_ucsc_join.txt')
#        ucsc_join(test_file)
#
    def tearDown(self):
        self.remove_working_dir()


