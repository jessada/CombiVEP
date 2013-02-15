import unittest
import os
import numpy as np
import combivep.settings as cbv_const
from combivep.engine.mlp import Mlp
from combivep.engine.moc_dataset import DataSet
from combivep.engine.test.template import SafeEngineTester

class TestMlp(SafeEngineTester):


    def __init__(self, test_name):
        SafeEngineTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'trainer'

    def init_trainer_instance(self):
        pass

    def test_fix_random_weight(self):
        """

        check if the initial random values of weight matrixes are likely to be
        random.

        """
        training_data = DataSet(os.path.join(cbv_const.CBV_SAMPLE_DATASET_DIR,
                                             'dummy_training_dataset'))
        mlp = Mlp(training_data.n_features, seed=20)
        self.assertEqual(round(mlp.weights1[0][1], 4),
                         0.0090,
                         msg='MLP is not ready for test because the random value is not fix')
        self.assertEqual(round(mlp.weights1[0][0], 4),
                         0.0059,
                         msg='MLP is not ready for test because the random value is not fix')

    def test_forward_propagation(self):
        """

        check if the matrix multiplications in forward propagationin are
        working properly.

        """
        training_data = DataSet(os.path.join(cbv_const.CBV_SAMPLE_DATASET_DIR,
                                             'dummy_training_dataset'))
        mlp = Mlp(training_data.n_features, seed=20)
        out = mlp.forward_propagation(training_data)
        self.assertEqual(round(out[0][0], 4),
                         0.5022,
                         msg='forward propagation does not functional properly')

    def test_one_round_forward_backward_weight_update(self):
        """

        to see if can correctly run one round of "forward", "backward" and
        "weight update"

        """
        training_data = DataSet(os.path.join(cbv_const.CBV_SAMPLE_DATASET_DIR,
                                             'dummy_training_dataset'))
        mlp = Mlp(training_data.n_features, seed=20)
        mlp.forward_propagation(training_data)
        mlp.backward_propagation(training_data)
        weights1, weights2 = mlp.weight_update(training_data)
        self.assertEqual(round(weights1[0][0], 4),
                         0.0059,
                         msg='one round of forward propagation, backward propagation and weight update, does not functional properly')


