import unittest
import os
import combivep.settings as cbv_const
from combivep.template import SafeTester
from combivep.template import RiskyTester


class SafeGeneralTester(SafeTester):
    """ General template for safe "base" modules testing """

    def __init__(self, test_name):
        SafeTester.__init__(self, test_name)

    def set_dir(self):
        self.working_dir = os.path.join(os.path.join(os.path.join(os.path.dirname(__file__),
                                                                  'tmp'),
                                                     self.test_class),
                                        self.test_function)
        self.data_dir = os.path.join(os.path.join(os.path.dirname(__file__),
                                                  'data'),
                                     self.test_class)


class RiskyGeneralTester(RiskyTester):
    """ General template for risky "base" modules testing """

    def __init__(self, test_name):
        RiskyTester.__init__(self, test_name)

    def set_dir(self):
        self.working_dir = cbv_const.COMBIVEP_WORKING_DIR
        self.data_dir = os.path.join(os.path.join(os.path.dirname(__file__),
                                                  'big_data'),
                                     self.test_class)
