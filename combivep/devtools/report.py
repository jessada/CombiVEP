import numpy as np
import os
import combivep.settings as cbv_const
import combivep.devtools.settings as dev_const
import matplotlib.pyplot as plt
from combivep.devtools.utils import filter_cbv_data
from combivep.devtools.utils import calculate_roc
from combivep.devtools.utils import print_preproc
from combivep.devtools.utils import fast_training
from combivep.devtools.utils import fast_predict
from combivep.devtools.utils import generate_roc_auc_figures
from combivep.devtools.utils import generate_scores_dist_figure
from combivep.devtools.utils import generate_precision_measurements_figures
from combivep.devtools.utils import measure_precision
from combivep.devtools.settings import PRECISION_MEASURES
from combivep.devtools.settings import PREDICTOR_NAMES
from combivep.app import predict_deleterious_probability


def filter_all_cbv():
    filter_cbv_data(os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                                 'training.cbv'))
    filter_cbv_data(os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                                 'test.cbv'))


def demo_training():
    fast_training(os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                               'training.cbv.scores'),
                  params_out_file=dev_const.PUB_PARAM_FILE,
                  random_seed=cbv_const.DEMO_SEED,
                  figure_dir=dev_const.PUB_FIGS_DIR,
                  )


def demo_predicting():
    fast_predict(os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                              'test.cbv.scores'),
                 params_file=dev_const.PUB_PARAM_FILE,
                 file_type=cbv_const.FILE_TYPE_CBV,
                 output_file=dev_const.PUB_RAW_PREDICTION_RESULT,
                 cfg_file=cbv_const.CBV_CFG_FILE,
                 )


def generate_figures():
    params = {'legend.fontsize': 9}
    plt.rcParams.update(params)
    header_data = np.loadtxt(dev_const.PUB_CONDEL_PREDICTION_RESULT,
                             dtype='S20'
                             )[:, :5]
    scores_data = np.loadtxt(dev_const.PUB_CONDEL_PREDICTION_RESULT,
                             dtype='S20'
                             )[:, 5:13].astype(np.float)
    patho_scores = scores_data[header_data[:, 4] == '1']
    neutr_scores = scores_data[header_data[:, 4] == '0']
    predictor_colors = ('k',
                        'm',
                        'c',
                        'g',
                        'b',
                        'coral',
                        'darkred',
                        'r',
                        )
    #plot roc and auc curve
    figs = generate_roc_auc_figures(plt,
                                    patho_scores,
                                    neutr_scores,
                                    PREDICTOR_NAMES,
                                    predictor_colors)
    auc_fig, roc_fig = figs

    #plot scores distribution
    scores_dist_fig = generate_scores_dist_figure(plt,
                                                  patho_scores,
                                                  neutr_scores)

    #plot precision performance
    figs = generate_precision_measurements_figures(plt,
                                                   patho_scores,
                                                   neutr_scores,
                                                   PREDICTOR_NAMES,
                                                   predictor_colors)
    result_classes_fig, precision_fig = figs

    return (roc_fig,
            auc_fig,
            scores_dist_fig,
            result_classes_fig,
            precision_fig,
            )


def generate_preproc_report():
    n_original_neutr = 21561
    n_original_patho = 21401
    print_preproc("", "Neutral SNPs", "Pathogenic SNPs")
    print_preproc("Parsed VariBench samples",
                  str(n_original_neutr),
                  str(n_original_patho))

    n_uncertain_neutr = 47
    n_uncertain_patho = 52
    print_preproc("Uncertain",
                  str(n_uncertain_neutr),
                  str(n_uncertain_patho))

    clean_training_file = os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                                       'training.cbv.clean')
    clean_test_file = os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                                   'test.cbv.clean')
    data = np.loadtxt(clean_training_file, dtype='S20')
    n_clean_patho_training = data[data[:, 4] == '1'].shape[0]
    n_clean_neutr_training = data[data[:, 4] == '0'].shape[0]
    data = np.loadtxt(clean_test_file, dtype='S20')
    n_clean_patho_test = data[data[:, 4] == '1'].shape[0]
    n_clean_neutr_test = data[data[:, 4] == '0'].shape[0]
    print_preproc("Unknown reference",
                  str(n_original_neutr - (n_clean_neutr_test+n_clean_neutr_training+n_uncertain_neutr)),
                  str(n_original_patho - (n_clean_patho_test+n_clean_patho_training+n_uncertain_patho)),
                  )

    scores_training_file = os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                                        'training.cbv.scores')
    scores_test_file = os.path.join(cbv_const.CBV_SAMPLE_CBV_DIR,
                                    'test.cbv.scores')
    data = np.loadtxt(scores_training_file, dtype='S20')
    n_scores_patho_training = data[data[:, 4] == '1'].shape[0]
    n_scores_neutr_training = data[data[:, 4] == '0'].shape[0]
    data = np.loadtxt(scores_test_file, dtype='S20')
    n_scores_patho_test = data[data[:, 4] == '1'].shape[0]
    n_scores_neutr_test = data[data[:, 4] == '0'].shape[0]
    print_preproc("Unidentified by effect predictors",
                  str((n_clean_neutr_test+n_clean_neutr_training) - (n_scores_neutr_test+n_scores_neutr_training)),
                  str((n_clean_patho_test+n_clean_patho_training) - (n_scores_patho_test+n_scores_patho_training)),
                  )
    print_preproc("Ready to be used by CombiVEP & Condel",
                  str(n_scores_neutr_test+n_scores_neutr_training),
                  str(n_scores_patho_test+n_scores_patho_training),
                  )
    print_preproc("Training dataset",
                  str(n_scores_neutr_training),
                  str(n_scores_patho_training),
                  )
    print_preproc("Test dataset",
                  str(n_scores_neutr_test),
                  str(n_scores_patho_test),
                  )


def generate_precision_performance_report():
    prediction_result = np.loadtxt(dev_const.PUB_CONDEL_PREDICTION_RESULT,
                                   dtype='S20')
    pos_data = prediction_result[prediction_result[:, 4] == '1']
    neg_data = prediction_result[prediction_result[:, 4] == '0']
    pos_samples = pos_data[:, 5:13].astype(np.float)
    neg_samples = neg_data[:, 5:13].astype(np.float)
    precision_performances = []
    for i in xrange(neg_samples.shape[1]):
        precision_performances.append(measure_precision(0.5,
                                                        pos_samples[:, i],
                                                        neg_samples[:, i]))

    print "%-20s" % "", "".join(map(lambda x: "%11s" % x, PREDICTOR_NAMES))
    for i in xrange(len(PRECISION_MEASURES)):
        print "%-20s" % PRECISION_MEASURES[i], "".join(map(lambda x: "%11s" % ("%.4f" % x[i]), precision_performances))
