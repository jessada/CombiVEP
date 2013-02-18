import unittest
import shutil
import os
import combivep.settings as cbv_const
from combivep.engine.moc_dataset import DataSet
from combivep.engine.test.template import SafeEngineTester
from combivep.engine.wrapper import Trainer
from combivep.engine.wrapper import Predictor


class TestTrainer(SafeEngineTester):

    def __init__(self, test_name):
        SafeEngineTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'trainer'

    def init_trainer_instance(self):
        pass

    def test_trainer(self):
        """

        to see if it can produce parameters file and produce figure

        """
        self.individual_debug = True
        self.init_test('test_trainer')
        self.init_trainer_instance()
        training_data   = DataSet(os.path.join(cbv_const.CBV_SAMPLE_DATASET_DIR,
                                               'training_dataset'))
        validation_data = DataSet(os.path.join(cbv_const.CBV_SAMPLE_DATASET_DIR,
                                               'validation_dataset'))
        trainer = Trainer(training_data,
                          validation_data,
                          seed=20,
                          n_hidden_nodes=7,
                          figure_dir=self.working_dir)
        trainer.train(iterations=50)

        params_file = os.path.join(self.working_dir,
                                   'params.npz')
        trainer.export_best_parameters(params_file=params_file)
        self.assertTrue(os.path.exists(params_file),
                        msg='Trainer does not functional properly')
        figure_file = os.path.join(self.working_dir,
                                   '07.eps')
        self.assertTrue(os.path.exists(figure_file),
                        msg='Trainer does not functional properly')

    def tearDown(self):
        self.remove_working_dir()


class TestPredictor(SafeEngineTester):

    def __init__(self, test_name):
        SafeEngineTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'predictor'

    def init_predictor_instance(self):
        pass

    def test_predictor(self):
        """

        to see if it can correctly predict SNPs in feature-vector format

        """
        self.individual_debug = True
        self.init_test('test_predictor')
        self.init_predictor_instance()

        predictor = Predictor()
        test_data = DataSet(os.path.join(cbv_const.CBV_SAMPLE_DATASET_DIR,
                                         'test_dataset'))
        params_file = os.path.join(self.data_dir,
                                   'params.npz')
        predictor.import_parameters(params_file=params_file)
        out = predictor.predict(test_data)
        self.assertEqual(round(out[0][0], 4),
                         0.2729,
                         msg='Predictor does not functional properly')
