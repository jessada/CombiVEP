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
                                    output_file=devtools_settings.PUBLICATION_PREDICTION_RESULT,
                                    config_file=combivep_settings.COMBIVEP_CONFIGURATION_FILE,
                                    )

def generate_performance_report():
    data = np.loadtxt(devtools_settings.PUBLICATION_PREDICTION_RESULT, dtype='S20')
    false_positive_rates, true_positive_rates = calculate_roc(data[data[:, 4] == '1'][:, 5:11].astype(np.float), 
                                                              data[data[:, 4] == '0'][:, 5:11].astype(np.float), 
                                                              np.linspace(0, 1, 101))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(false_positive_rates[:, 0], true_positive_rates[:, 0], 'k',
            false_positive_rates[:, 1], true_positive_rates[:, 1], 'm--',
            false_positive_rates[:, 2], true_positive_rates[:, 2], 'c--',
            false_positive_rates[:, 3], true_positive_rates[:, 3], 'g--',
            false_positive_rates[:, 4], true_positive_rates[:, 4], 'k--',
            false_positive_rates[:, 5], true_positive_rates[:, 5], 'r--',
            )
    ax.set_ylabel('true positive rate')
    ax.set_xlabel('false positive rate')
    fig.savefig(devtools_settings.PUBLICATION_ROC_CURVE_FIGURE, bbox_inches='tight', pad_inches=0.01)
    return devtools_settings.PUBLICATION_ROC_CURVE_FIGURE







