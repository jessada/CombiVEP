import os
import unittest
from combivep.test.template import SafeGeneralTester
import combivep.settings as combivep_settings
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
#        self.individual_debug = True
        self.init_test('train_combivep_using_cbv_data')
        test_file    = os.path.join(self.data_dir, 'test_train_combivep.cbv')
        params_file  = os.path.join(self.working_dir, 'params.npz')
        #run test
        train_combivep_using_cbv_data(test_file,
                                      random_seed=1,
                                      n_hidden_nodes=7,
                                      iterations=50,
                                      figure_dir=self.working_dir,
                                      params_out_file=params_file,
                                      config_file=combivep_settings.COMBIVEP_CENTRAL_TEST_CONFIGURATION_FILE)
        self.assertTrue(os.path.exists(params_file), msg='Trainer does not functional properly')
        figure_file  = os.path.join(self.working_dir, '07.eps')
        self.assertTrue(os.path.exists(figure_file), msg='Trainer does not functional properly')

    def test_predict_deleterious_probability_cbv(self):
        #init
#        self.individual_debug = True
        self.init_test('predict_deleterious_probability_cbv')
        test_file = os.path.join(self.data_dir, 'test_test_combivep.cbv')
        output_file  = os.path.join(self.working_dir, 'cbv_output.txt')
        #run test
        args = {}
        args["params_file"] = combivep_settings.COMBIVEP_CENTRAL_TEST_PARAMETER_FILE
        args["file_type"]   = combivep_settings.FILE_TYPE_CBV
        args["config_file"] = combivep_settings.COMBIVEP_CENTRAL_TEST_CONFIGURATION_FILE
        args["output_file"] = output_file
        predict_deleterious_probability(test_file,
                                        **args
                                        )
        self.assertTrue(os.path.exists(output_file), msg='Predictor does not functional properly')
        f = open(output_file, 'r')
#        self.assertEqual(f.readline().strip(), '#CHROM\tPOS\tREF\tALT\tACTUAL_DELETERIOUS_EFFECT\tPREDICTED_DELETERIOUS_PROBABILITY')
        self.assertEqual(f.readline().strip(), '#CHROM\tPOS\tREF\tALT\tACTUAL_DELETERIOUS_EFFECT\tPREDICTED_DELETERIOUS_PROBABILITY\tPHYLOP_SCORE\tSIFT_SCORE\tPP2_SCORE\tLRT_SCORT\tMT_SCORE\tGERP_SCORE')
        self.assertEqual(f.readline().strip(), '1\t35227264\tT\tC\t1\t0.2605\t0.968087\t0.96\t0.031\t1.000000\t0.838867\t4.45')
        f.close()

    def test_predict_deleterious_probability_vcf(self):
        #init
        self.init_test('predict_deleterious_probability_vcf')
        test_file = os.path.join(self.data_dir, 'test_test_combivep.vcf')
        output_file  = os.path.join(self.working_dir, 'vcf_output.txt')
        #run test
        predict_deleterious_probability(test_file,
                                        params_file=combivep_settings.COMBIVEP_CENTRAL_TEST_PARAMETER_FILE,
                                        file_type=combivep_settings.FILE_TYPE_VCF,
                                        config_file=combivep_settings.COMBIVEP_CENTRAL_TEST_CONFIGURATION_FILE,
                                        output_file=output_file,
                                        )
        self.assertTrue(os.path.exists(output_file), msg='Predictor does not functional properly')
        f = open(output_file, 'r')
#        self.assertEqual(f.readline().strip(), '#CHROM\tPOS\tREF\tALT\tACTUAL_DELETERIOUS_EFFECT\tPREDICTED_DELETERIOUS_PROBABILITY')
        self.assertEqual(f.readline().strip(), '#CHROM\tPOS\tREF\tALT\tACTUAL_DELETERIOUS_EFFECT\tPREDICTED_DELETERIOUS_PROBABILITY\tPHYLOP_SCORE\tSIFT_SCORE\tPP2_SCORE\tLRT_SCORT\tMT_SCORE\tGERP_SCORE')
        self.assertEqual(f.readline().strip(), '3\t361508\tC\tT\tNone\t0.0001\t0.024209\t0.0\t0\t0.950380\t0.000019\t-2.66')
        f.close()

    def tearDown(self):
        self.remove_working_dir()


