import os
import unittest
import combivep.settings as cbv_const
from combivep.test.template import SafeGeneralTester
from combivep.app import train_combivep_using_cbv_data
from combivep.app import predict_deleterious_probability


class TestApp(SafeGeneralTester):


    def __init__(self, test_name):
        SafeGeneralTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'app'

    def init_configure_instance(self):
        pass

    def test_train_combivep_using_cbv_data(self):
        #init
        self.init_test('train_combivep_using_cbv_data')
        test_file    = os.path.join(self.data_dir,
                                    'test_train_combivep.cbv')
        params_file  = os.path.join(self.working_dir,
                                    'params.npz')
        #run test
        train_combivep_using_cbv_data(test_file,
                                      random_seed=1,
                                      n_hidden_nodes=7,
                                      iterations=50,
                                      figure_dir=self.working_dir,
                                      params_out_file=params_file,
                                      cfg_file=cbv_const.CBV_SAMPLE_CFG_FILE)
        self.assertTrue(os.path.exists(params_file),
                        msg='Trainer does not functional properly')
        figure_file  = os.path.join(self.working_dir,
                                    '07.eps')
        self.assertTrue(os.path.exists(figure_file),
                        msg='Trainer does not functional properly')

    def test_predict_deleterious_probability_cbv(self):
        #init
        self.init_test('predict_deleterious_probability_cbv')
        test_file   = os.path.join(self.data_dir,
                                   'test_test_combivep.cbv')
        output_file = os.path.join(self.working_dir,
                                   'cbv_output.txt')
        #run test
        kwargs = {}
        kwargs["params_file"] = cbv_const.CBV_SAMPLE_PARAM_FILE
        kwargs["file_type"]   = cbv_const.FILE_TYPE_CBV
        kwargs["cfg_file"]    = cbv_const.CBV_SAMPLE_CFG_FILE
        kwargs["output_file"] = output_file
        predict_deleterious_probability(test_file,
                                        **kwargs)
        self.assertTrue(os.path.exists(output_file),
                        msg='Predictor does not functional properly')
        f = open(output_file, 'r')
        self.assertEqual(f.readline().strip(),
                         "#"+"\t".join(cbv_const.PREDICTION_OUT_COLS_HEADER))
        self.assertEqual(f.readline().strip(),
                         '1\t35227264\tT\tC\t1\t0.2605\t0.968087\t0.96\t0.031\t1.000000\t0.838867\t4.45')
        f.close()

    def test_predict_deleterious_probability_vcf(self):
        #init
        self.init_test('predict_deleterious_probability_vcf')
        test_file   = os.path.join(self.data_dir,
                                   'test_test_combivep.vcf')
        output_file = os.path.join(self.working_dir,
                                   'vcf_output.txt')
        #run test
        kwargs = {}
        kwargs["params_file"] = cbv_const.CBV_SAMPLE_PARAM_FILE
        kwargs["file_type"]   = cbv_const.FILE_TYPE_VCF
        kwargs["cfg_file"]    = cbv_const.CBV_SAMPLE_CFG_FILE
        kwargs["output_file"] = output_file
        predict_deleterious_probability(test_file,
                                        **kwargs)
        self.assertTrue(os.path.exists(output_file),
                        msg='Predictor does not functional properly')
        f = open(output_file, 'r')
        self.assertEqual(f.readline().strip(),
                         "#"+"\t".join(cbv_const.PREDICTION_OUT_COLS_HEADER))
        self.assertEqual(f.readline().strip(),
                         '3\t361508\tC\tT\tNone\t0.0001\t0.024209\t0.0\t0\t0.950380\t0.000019\t-2.66')
        f.close()

    def tearDown(self):
        self.remove_working_dir()


