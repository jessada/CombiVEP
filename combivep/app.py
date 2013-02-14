import argparse
import sys
import os
import numpy as np
import combivep.settings as cbv_const
from combivep.preproc.dataset import DataSetManager
from combivep.engine.wrapper import Trainer
from combivep.engine.wrapper import Predictor
from combivep.refdb.control import UcscController
from combivep.refdb.control import LjbController


def app_combivep_reference_updater():
    #update LJB reference database
    ljb_controller = LjbController()
    ljb_controller.update()
    #update UCSC reference database
    ucsc_controller = UcscController()
    ucsc_controller.update()

def app_combivep_trainer():
    argp = argparse.ArgumentParser(description="An application to train the consensus predictor using Multilayer Perceptron Learning technique.")
    tmp_help=[]
    tmp_help.append("A file with list of SNPs in CBV (CombiVEP) format.")
    tmp_help.append("CBV format consists of five fields: CHROM, POS, REF, ALT, ACTUAL_DELETERIOUS_EFFECT.")
    tmp_help.append("Each field is separated by a tab.")
    tmp_help.append("SNP Position(POS) is 1-based index.")
    argp.add_argument('training_data_file', help=' '.join(tmp_help))
    args = argp.parse_args()
    train_combivep_using_cbv_data(args.training_data_file)

def app_combivep_predictor():
    argp = argparse.ArgumentParser(description="An application to predict how likely the given SNPs will have deleterious effect")
    tmp_help=[]
    tmp_help.append("A file with list of SNPs.")
    tmp_help.append("It can be either in standard VCF format or CBV (CombiVEP) format.")
    tmp_help.append("CBV format consists of five fields: CHROM, POS, REF, ALT, ACTUAL_DELETERIOUS_EFFECT.")
    tmp_help.append("Each field is separated by a tab.")
    tmp_help.append("SNP Position(POS) is 1-based index.")
    argp.add_argument('input_file', help=' '.join(tmp_help))
    argp.add_argument('-F', help='Input format (%(choices)s). Default is in %(default)s format.', 
                            choices=[cbv_const.FILE_TYPE_VCF, cbv_const.FILE_TYPE_CBV], 
                            metavar='FORMAT', 
                            default=cbv_const.FILE_TYPE_VCF,
                            dest='input_format'
                            )
    args = argp.parse_args()
    predict_deleterious_probability(args.input_file, file_type=args.input_format)

def train_combivep_using_cbv_data(training_data_file, 
                                  params_out_file=cbv_const.USER_PARAMS_FILE,
                                  random_seed=cbv_const.DEFAULT_SEED,
                                  n_hidden_nodes=cbv_const.DEFAULT_HIDDEN_NODES,
                                  figure_dir=cbv_const.DEFAULT_FIGURE_DIR,
                                  iterations=cbv_const.DEFAULT_ITERATIONS,
                                  cfg_file=cbv_const.CBV_CFG_FILE,
                                  ):
    """

    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
    CBV has 5 fields, CHROM, POS, REF, ALT, EFFECT (1=deleterious, 0=neutral). All are tab separated
    Required arguments
    - neutral_data_file : list of SNPs with no harmful effect, CBV format
    - pathognice_data_file : list of SNPs with deleterious effect, CBV format

    """
    #pre-processing dataset
    print >> sys.stderr, 'pre-processing dataset, this may take a while (around 750 SNPs/mins). . . . '
    dm = DataSetManager(cfg_file=cfg_file)
    dm.load_data(training_data_file, file_type=cbv_const.FILE_TYPE_CBV)
    dm.validate_data()
    dm.calculate_scores()
    dm.set_shuffle_seed(random_seed)
    dm.shuffle_data()
    dm.partition_data()

    #partition data
    training_dataset      = dm.get_training_data() 
    validation_dataset    = dm.get_validation_data() 

    #train !!!
    print >> sys.stderr, 'Training CombiVEP, please wait (around 500 SNPs/mins) . . . . '
    trainer = Trainer(training_dataset, validation_dataset, random_seed, n_hidden_nodes, figure_dir)
    trainer.train(iterations)
    if not os.path.exists(cbv_const.USER_PARAMS_DIR):
        os.makedirs(cbv_const.USER_PARAMS_DIR)
    trainer.export_best_parameters(params_out_file)

def predict_deleterious_probability(SNPs_file,
                                    params_file=cbv_const.USER_PARAMS_FILE,
                                    file_type=cbv_const.FILE_TYPE_VCF,
                                    output_file=None,
                                    cfg_file=cbv_const.CBV_CFG_FILE,
                                    ):
    """

    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
    CBV has 5 fields, CHROM, POS, REF, ALT, EFFECT (1=deleterious, 0=neutral). All are tab separated
    Required arguments
    - SNPs_file : list of SNPs to be predicted, can be either VCF or CBV (default is VCF)

    """
    #pre-processing test dataset
    print >> sys.stderr, 'pre-processing dataset, this may take a while (around 750 SNPs/mins). . . . '
    dm = DataSetManager(cfg_file=cfg_file)
    dm.load_data(SNPs_file, file_type=file_type)
    dm.validate_data()
    dm.calculate_scores()

    #predict
    predictor = Predictor()
    predictor.import_parameters(params_file=params_file)
    out = (np.array(predictor.predict(dm.dataset)).reshape(-1,))

    #print output
    if output_file is not None:
        sys.stdout = open(output_file, 'w')
    print >> sys.stdout, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ("CHROM", "POS", "REF", "ALT", "ACTUAL_DELETERIOUS_EFFECT", "PREDICTED_DELETERIOUS_PROBABILITY", "PHYLOP_SCORE", "SIFT_SCORE", "PP2_SCORE", "LRT_SCORT", "MT_SCORE", "GERP_SCORE")
    for i in xrange(len(dm.dataset)):
        print >> sys.stdout, "%s\t%s\t%s\t%s\t%s\t%6.4f\t%s\t%s\t%s\t%s\t%s\t%s" % (dm.dataset[i][cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_CHROM],
                                                        dm.dataset[i][cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_POS],
                                                        dm.dataset[i][cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_REF],
                                                        dm.dataset[i][cbv_const.KEY_SNP_INFO_SECTION][cbv_const.KEY_ALT],
                                                        dm.dataset[i][cbv_const.KEY_PREDICTION_SECTION][cbv_const.KEY_TARGETS],
                                                        out[i],
                                                        dm.dataset[i][cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_PHYLOP_SCORE],
                                                        dm.dataset[i][cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_SIFT_SCORE],
                                                        dm.dataset[i][cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_PP2_SCORE],
                                                        dm.dataset[i][cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_LRT_SCORE],
                                                        dm.dataset[i][cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_MT_SCORE],
                                                        dm.dataset[i][cbv_const.KEY_SCORES_SECTION][cbv_const.KEY_GERP_SCORE],
                                                        )
    sys.stdout = sys.__stdout__

