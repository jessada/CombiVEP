import matplotlib.pyplot as plt
import numpy as np
import os
from combivep.devtools.utils import filter_cbv_data
from combivep.devtools.utils import calculate_roc
from combivep.app import train_combivep_using_cbv_data
import combivep.settings as combivep_settings
import combivep.devtools.settings as devtools_settings
from combivep.app import predict_deleterious_probability

def filter_all_cbv():
    filter_cbv_data(os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'training.cbv'))
    filter_cbv_data(os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'test.cbv'))

def demo_training():
    train_combivep_using_cbv_data(os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'training.cbv.scores'),
                                  params_out_file=devtools_settings.PUBLICATION_PARAMETER_FILE,
                                  random_seed=combivep_settings.DEMO_SEED,
                                  figure_dir=devtools_settings.PUBLICATION_FIGURES_DIR,
                                  )

def demo_predicting():
    predict_deleterious_probability(os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'test.cbv.scores'),
                                    params_file=devtools_settings.PUBLICATION_PARAMETER_FILE,
                                    file_type=combivep_settings.FILE_TYPE_CBV,
                                    output_file=devtools_settings.PUBLICATION_RAW_PREDICTION_RESULT,
                                    config_file=combivep_settings.COMBIVEP_CONFIGURATION_FILE,
                                    )

def generate_performance_report():
    data      = np.loadtxt(devtools_settings.PUBLICATION_CONDEL_PREDICTION_RESULT, dtype='S20')
    min_value = np.amin(data[:, 5:13].astype(np.float))
    max_value = np.amax(data[:, 5:13].astype(np.float))

    #produce roc data from CombiVEP, Phylop, SIFT, PP2, LRT, MT, GERP, Condel
    fp_rates, tp_rates = calculate_roc(data[data[:, 4] == '1'][:, 5:13].astype(np.float),
                                       data[data[:, 4] == '0'][:, 5:13].astype(np.float),
                                       np.linspace(min_value, max_value, 5001))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(fp_rates[:, 0], tp_rates[:, 0], 'k',
            fp_rates[:, 1], tp_rates[:, 1], 'm--',
            fp_rates[:, 2], tp_rates[:, 2], 'c--',
            fp_rates[:, 3], tp_rates[:, 3], 'g--',
            fp_rates[:, 4], tp_rates[:, 4], 'k--',
            fp_rates[:, 5], tp_rates[:, 5], 'r--',
            fp_rates[:, 6], tp_rates[:, 6], 'b--',
            fp_rates[:, 7], tp_rates[:, 7], 'r',
            )
    ax.set_ylabel('true positive rate')
    ax.set_xlabel('false positive rate')
    fig.savefig(devtools_settings.PUBLICATION_ROC_CURVE_FIGURE, bbox_inches='tight', pad_inches=0.01)
    return devtools_settings.PUBLICATION_ROC_CURVE_FIGURE







