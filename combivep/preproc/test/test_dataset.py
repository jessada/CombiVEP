import unittest
import os
from combivep.preproc.test.template import SafePreProcTester
import combivep.settings as cbv_const
from combivep.preproc.dataset import DataSetManager


class TestDataSetManager(SafePreProcTester):


    def __init__(self, test_name):
        SafePreProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'dataset_manager'

    def init_dataset_instance(self):
        self.__dataset_manager = DataSetManager(config_file=cbv_const.CBV_CENTRAL_TEST_CFG_FILE)

    def test_vcf_load(self):
        self.init_test('test_vcf_load')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_vcf_load.vcf')
        self.__dataset_manager.load_data(test_file)
        self.assertEqual(len(self.__dataset_manager.dataset), 10, 'DataSetManager does not load VCF data correctly')

    def test_cbv_load(self):
        self.init_test('test_cbv_load')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_cbv_load.cbv')
        self.__dataset_manager.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        self.assertEqual(len(self.__dataset_manager.dataset), 11, 'DataSetManager does not load CBV data correctly')

    def test_validate_data(self):
        self.init_test('test_validate_data')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_vcf_load.vcf')
        self.__dataset_manager.load_data(test_file)
        self.__dataset_manager.validate_data()
        self.assertEqual(len(self.__dataset_manager.dataset), 7, 'DataSetManager does not clean data correctly')

    def test_calculate_scores(self):
        self.init_test('test_calculate_scores')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_calculate_scores.vcf')
        self.__dataset_manager.load_data(test_file)
        self.__dataset_manager.validate_data()
        self.__dataset_manager.calculate_scores()
        self.assertEqual(len(self.__dataset_manager.dataset), 3, 'DataSetManager does not calculate scores properly')

    def test_shuffle_data(self):
        self.init_test('test_shuffle_data')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_shuffle.vcf')
        self.__dataset_manager.load_data(test_file)
        self.__dataset_manager.validate_data()
        self.__dataset_manager.calculate_scores()
        self.assertEqual(self.__dataset_manager.dataset[5][cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_POS], '190999917', 'DataSetManager does not calculate scores properly')
        self.__dataset_manager.set_shuffle_seed(cbv_const.DEMO_SEED)
        self.__dataset_manager.shuffle_data()
        self.assertNotEqual(self.__dataset_manager.dataset[5][cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_POS], '190999917', 'DataSetManager may not shuffle data correctly')

    def test_partition_data(self):
        self.init_test('test_partition_data')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_partition.cbv')
        self.__dataset_manager.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        self.__dataset_manager.validate_data()
        self.__dataset_manager.calculate_scores()
        self.__dataset_manager.set_shuffle_seed(cbv_const.DEMO_SEED)
        self.__dataset_manager.shuffle_data()
        self.__dataset_manager.partition_data()
        self.assertEqual(len(self.__dataset_manager.get_training_data()), 11, 'DataSetManager does not correctly partition data')
        self.assertEqual(len(self.__dataset_manager.get_validation_data()), 3, 'DataSetManager does not correctly partition data')

    def test_vcf_dataset(self):
        self.init_test('test_dataset')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_shuffle.vcf')
        self.__dataset_manager.load_data(test_file)
        self.__dataset_manager.validate_data()
        self.__dataset_manager.calculate_scores()
        self.__dataset_manager.set_shuffle_seed(cbv_const.DEMO_SEED)
        self.__dataset_manager.shuffle_data()
        self.__dataset_manager.partition_data()
        training_dataset = self.__dataset_manager.get_training_data()
        self.assertEqual(training_dataset.n_features, 6, msg='Dataset does not functional properly')
        self.assertEqual(training_dataset.n_data, 8, msg='Dataset does not functional properly')

    def test_cbv_dataset(self):
        self.init_test('test_dataset')
        self.init_dataset_instance()
        test_file = os.path.join(self.data_dir, 'test_cbv_dataset.cbv')
        self.__dataset_manager.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        self.__dataset_manager.validate_data
        self.__dataset_manager.calculate_scores()
        self.__dataset_manager.set_shuffle_seed(cbv_const.DEMO_SEED)
        self.__dataset_manager.shuffle_data()
        self.__dataset_manager.partition_data()
        training_dataset = self.__dataset_manager.get_training_data()
        self.assertEqual(training_dataset.n_features, 6, msg='Dataset does not functional properly')
        self.assertEqual(training_dataset.n_data, 15, msg='Dataset does not functional properly')

    def test_add_dataset(self):
        self.init_test('test_add_dataset')
        test_file = os.path.join(self.data_dir, 'test_add_dataset1.cbv')
        self.dataset_manager1 = DataSetManager(config_file=cbv_const.CBV_CENTRAL_TEST_CFG_FILE)
        self.dataset_manager1.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        self.dataset_manager1.validate_data()
        self.dataset_manager1.calculate_scores()
        self.dataset_manager1.shuffle_data()
        self.dataset_manager1.partition_data()
        training_dataset1 = self.dataset_manager1.get_training_data()
        test_file = os.path.join(self.data_dir, 'test_add_dataset2.cbv')
        self.dataset_manager2 = DataSetManager(config_file=cbv_const.CBV_CENTRAL_TEST_CFG_FILE)
        self.dataset_manager2.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        self.dataset_manager2.validate_data()
        self.dataset_manager2.calculate_scores()
        self.dataset_manager2.shuffle_data()
        self.dataset_manager2.partition_data()
        training_dataset2 = self.dataset_manager2.get_training_data()
        combine_dataset = training_dataset1 + training_dataset2
        self.assertEqual(len(combine_dataset), 10, 'DataSetManager does not load VCF data correctly')

    def tearDown(self):
        self.remove_working_dir()






