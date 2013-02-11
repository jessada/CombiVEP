import matplotlib.pyplot as plt
import numpy as np
import os
from combivep.devtools.utils import filter_cbv_data
from combivep.devtools.utils import calculate_roc
from combivep.app import train_combivep_using_cbv_data
import combivep.settings as combivep_settings
import combivep.devtools.settings as devtools_settings
from combivep.app import predict_deleterious_probability
from sklearn.metrics import auc

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
    header_data = np.loadtxt(devtools_settings.PUBLICATION_CONDEL_PREDICTION_RESULT, dtype='S20')[:, :5]
    scores_data = np.loadtxt(devtools_settings.PUBLICATION_CONDEL_PREDICTION_RESULT, dtype='S20')[:, 5:13].astype(np.float)
    min_value   = np.amin(scores_data)
    max_value   = np.amax(scores_data)
    predictor_names = ('CombiVEP',
                       'Phylop',
                       'SIFT',
                       'PP2',
                       'LRT',
                       'MT',
                       'GERP',
                       'Condel',
                       )
    predictor_colors = ('k',
                        'm',
                        'c',
                        'g',
                        'b',
                        'coral',
                        'darkred',
                        'r',
                        )

    #produce roc data from CombiVEP, Phylop, SIFT, PP2, LRT, MT, GERP, Condel
    fp_rates, tp_rates = calculate_roc(scores_data[header_data[:, 4] == '1'],
                                       scores_data[header_data[:, 4] == '0'],
                                       np.linspace(min_value, max_value, 5001))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in xrange(len(predictor_names)):
        ax.plot(fp_rates[:, i], tp_rates[:, i], predictor_colors[i], label=predictor_names[i])
    ax.set_ylabel('true positive rate')
    ax.set_xlabel('false positive rate')
    ax.legend(bbox_to_anchor=(0.9999, 0.0001), loc=4)
    fig.savefig(devtools_settings.PUBLICATION_ROC_FIGURE, bbox_inches='tight', pad_inches=0.01)

    #produce auc data from roc data
    fig  = plt.figure()
    aucs = []
    ind  = []
    ax   = fig.add_subplot(111)
    for i in xrange(len(predictor_names)):
        aucs.append(auc(fp_rates[:, i], tp_rates[:, i]))
        ind.append(0.5*(i+1)-0.4)
    ax.bar(ind, aucs, 0.3, color=predictor_colors)
    for i in xrange(len(aucs)):
        ax.text(ind[i], aucs[i] + 0.01, "%0.3f" % aucs[i])
    ax.set_ylim([0.7, 0.9])
    ax.set_xticks(np.array(ind) + 0.15)
    ax.set_xticklabels(predictor_names, rotation=30)
    fig.savefig(devtools_settings.PUBLICATION_AUC_FIGURE, bbox_inches='tight', pad_inches=0.01)

#    #plot scores distribution
#    fig        = plt.figure()
#    ax         = fig.add_subplot(211)
#    hist_range = (-0.005, 1.005)
#    pathogenic_hist, bins = np.histogram(scores_data[header_data[:, 4] == '1'][:, 0], bins = 100, range = hist_range)
#    neutral_hist, bins    = np.histogram(scores_data[header_data[:, 4] == '0'][:, 0], bins = 100, range = hist_range)
#    center = (bins[:-1]+bins[1:])/2
#    ax.plot(center, pathogenic_hist, 'r--')
#    ax.plot(center, neutral_hist, 'b--')
#    ax = fig.add_subplot(212)
#    pathogenic_hist, bins = np.histogram(scores_data[header_data[:, 4] == '1'][:, 7], bins = 100, range = hist_range)
#    neutral_hist, bins    = np.histogram(scores_data[header_data[:, 4] == '0'][:, 7], bins = 100, range = hist_range)
#    center = (bins[:-1]+bins[1:])/2
#    ax.plot(center, pathogenic_hist, 'r--')
#    ax.plot(center, neutral_hist, 'b--')
#    fig.savefig(devtools_settings.PUBLICATION_SCORE_DISTRIBUTION_FIGURE, bbox_inches='tight', pad_inches=0.01)
#
    return (devtools_settings.PUBLICATION_ROC_FIGURE,
            devtools_settings.PUBLICATION_AUC_FIGURE,
#            devtools_settings.PUBLICATION_SCORE_DISTRIBUTION_FIGURE,
            )







