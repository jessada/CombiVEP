import unittest
import os
import combivep.settings as cbv_const
from combivep.preproc.test.template import SafePreProcTester
from combivep.preproc.dataset import DataSetManager


class TestDataSetManager(SafePreProcTester):


    def __init__(self, test_name):
        SafePreProcTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'dataset_manager'

    def __create_dataset_manager_instance(self):
        dm = DataSetManager(cfg_file=cbv_const.CBV_CENTRAL_TEST_CFG_FILE)
        return dm

    def test_vcf_load(self):
        self.init_test('test_vcf_load')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_vcf_load.vcf')
        dm.load_data(test_file)
        self.assertEqual(len(dm.dataset),
                         10,
                         'DataSetManager does not load VCF data correctly')

    def test_cbv_load(self):
        self.init_test('test_cbv_load')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_cbv_load.cbv')
        dm.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        self.assertEqual(len(dm.dataset),
                         11,
                         'DataSetManager does not load CBV data correctly')

    def test_validate_data(self):
        self.init_test('test_validate_data')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_vcf_load.vcf')
        dm.load_data(test_file)
        dm.validate_data()
        self.assertEqual(len(dm.dataset),
                         7,
                         'DataSetManager does not clean data correctly')

    def test_calculate_scores(self):
        self.init_test('test_calculate_scores')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_calculate_scores.vcf')
        dm.load_data(test_file)
        dm.validate_data()
        dm.calculate_scores()
        self.assertEqual(len(dm.dataset),
                         3,
                         'DataSetManager does not calculate scores properly')

    def test_shuffle_data(self):
        self.init_test('test_shuffle_data')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_shuffle.vcf')
        dm.load_data(test_file)
        dm.validate_data()
        dm.calculate_scores()
        snp_info_5 = dm.dataset[5][cbv_const.KEY_SNP_INFO_SECTION]
        self.assertEqual(snp_info_5[cbv_const.KEY_POS],
                         '190999917',
                         'DataSetManager does not calculate scores properly')
        dm.set_shuffle_seed(cbv_const.DEMO_SEED)
        dm.shuffle_data()
        snp_info_5 = dm.dataset[5][cbv_const.KEY_SNP_INFO_SECTION]
        self.assertNotEqual(snp_info_5[cbv_const.KEY_POS],
                            '190999917',
                            'DataSetManager may not shuffle data correctly')

    def test_partition_data(self):
        self.init_test('test_partition_data')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_partition.cbv')
        dm.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        dm.validate_data()
        dm.calculate_scores()
        dm.set_shuffle_seed(cbv_const.DEMO_SEED)
        dm.shuffle_data()
        dm.partition_data()
        self.assertEqual(len(dm.get_training_data()),
                         11,
                         'DataSetManager does not correctly partition data')
        self.assertEqual(len(dm.get_validation_data()),
                         3,
                         'DataSetManager does not correctly partition data')

    def test_vcf_dataset(self):
        self.init_test('test_dataset')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_shuffle.vcf')
        dm.load_data(test_file)
        dm.validate_data()
        dm.calculate_scores()
        dm.set_shuffle_seed(cbv_const.DEMO_SEED)
        dm.shuffle_data()
        dm.partition_data()
        training_dataset = dm.get_training_data()
        self.assertEqual(training_dataset.n_features,
                         6,
                         msg='Dataset does not functional properly')
        self.assertEqual(training_dataset.n_data,
                         8,
                         msg='Dataset does not functional properly')

    def test_cbv_dataset(self):
        self.init_test('test_dataset')
        dm = self.__create_dataset_manager_instance()
        test_file = os.path.join(self.data_dir,
                                 'test_cbv_dataset.cbv')
        dm.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        dm.validate_data
        dm.calculate_scores()
        dm.set_shuffle_seed(cbv_const.DEMO_SEED)
        dm.shuffle_data()
        dm.partition_data()
        training_dataset = dm.get_training_data()
        self.assertEqual(training_dataset.n_features,
                         6,
                         msg='Dataset does not functional properly')
        self.assertEqual(training_dataset.n_data,
                         15,
                         msg='Dataset does not functional properly')

    def test_add_dataset(self):
        self.init_test('test_add_dataset')
        test_file = os.path.join(self.data_dir,
                                 'test_add_dataset1.cbv')
        dm1 = DataSetManager(cfg_file=cbv_const.CBV_CENTRAL_TEST_CFG_FILE)
        dm1.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        dm1.validate_data()
        dm1.calculate_scores()
        dm1.shuffle_data()
        dm1.partition_data()
        training_dataset1 = dm1.get_training_data()
        test_file = os.path.join(self.data_dir,
                                 'test_add_dataset2.cbv')
        dm2 = DataSetManager(cfg_file=cbv_const.CBV_CENTRAL_TEST_CFG_FILE)
        dm2.load_data(test_file, file_type=cbv_const.FILE_TYPE_CBV)
        dm2.validate_data()
        dm2.calculate_scores()
        dm2.shuffle_data()
        dm2.partition_data()
        training_dataset2 = dm2.get_training_data()
        combine_dataset = training_dataset1 + training_dataset2
        self.assertEqual(len(combine_dataset),
                         10,
                         'DataSetManager does not load VCF data correctly')

    def tearDown(self):
        self.remove_working_dir()






