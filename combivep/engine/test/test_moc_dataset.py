import unittest
import os.path
from combivep.engine.moc_dataset import DataSet
from combivep.template import Tester
import combivep.settings as cbv_const

class TestDataSet(Tester):
    """ to test dataset.py"""


    def setUp(self):
        pass

    def test_dataset_property(self):
        """

        check if number of rows, columns, and features are correctly counted.

        """
        dataset = DataSet(os.path.join(cbv_const.CBV_CENTRAL_TEST_DATASET_DIR,
                                                        'test_dataset'))
        self.assertEqual(dataset.n_data, 1718, msg='Dataset does not functional properly')
        self.assertEqual(dataset.n_cols, 8, msg='Dataset does not functional properly')
        self.assertEqual(dataset.n_features, 6, msg='Dataset does not functional properly')


