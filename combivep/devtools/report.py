import matplotlib.pyplot as plt
import numpy as np
import os
from combivep.devtools.utils import filter_cbv_data
from combivep.devtools.utils import calculate_roc
from combivep.devtools.utils import print_preproc
from combivep.devtools.utils import print_precision
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

def generate_figures():
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
    fig.savefig(devtools_settings.PUBLICATION_ROC_FIGURE, bbox_inches='tight', pad_inches=0.05)

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
    fig.savefig(devtools_settings.PUBLICATION_AUC_FIGURE, bbox_inches='tight', pad_inches=0.05)

    #plot scores distribution
    fig        = plt.figure()
    ax         = fig.add_subplot(211)
    hist_range = (-0.005, 1.005)
    pathogenic_hist, bins = np.histogram(scores_data[header_data[:, 4] == '1'][:, 0], bins = 100, range = hist_range)
    neutral_hist, bins    = np.histogram(scores_data[header_data[:, 4] == '0'][:, 0], bins = 100, range = hist_range)
    center = (bins[:-1]+bins[1:])/2
    ax.plot(center, pathogenic_hist, 'r--', label='pathogenic variants')
    ax.plot(center, neutral_hist, 'b--', label='neutral variants')
    ax.set_title('CombiVEP score distributuion')
    ax.set_ylabel('samples')
    ax.set_xlabel('score')
    ax.legend(bbox_to_anchor=(0.999, 0.999), loc=1)
    ax = fig.add_subplot(212)
    pathogenic_hist, bins = np.histogram(scores_data[header_data[:, 4] == '1'][:, 7], bins = 100, range = hist_range)
    neutral_hist, bins    = np.histogram(scores_data[header_data[:, 4] == '0'][:, 7], bins = 100, range = hist_range)
    center = (bins[:-1]+bins[1:])/2
    ax.plot(center, pathogenic_hist, 'r--', label='pathogenic variants')
    ax.plot(center, neutral_hist, 'b--', label='neutral variants')
    ax.set_title('Condel score distributuion')
    ax.set_ylabel('samples')
    ax.set_xlabel('score')
    ax.legend(bbox_to_anchor=(0.999, 0.999), loc=1)
    fig.tight_layout()
    fig.savefig(devtools_settings.PUBLICATION_SCORE_DISTRIBUTION_FIGURE, bbox_inches='tight', pad_inches=0.05)

    return (devtools_settings.PUBLICATION_ROC_FIGURE,
            devtools_settings.PUBLICATION_AUC_FIGURE,
            devtools_settings.PUBLICATION_SCORE_DISTRIBUTION_FIGURE,
            )

def generate_preproc_report():
    original_neutral_records    = 21561
    original_pathogenic_records = 21401
    print_preproc("", "Neutral SNPs", "Pathogenic SNPs")
    print_preproc("Parsed VariBench samples", str(original_neutral_records), str(original_pathogenic_records))

    uncertain_neutral_records    = 47
    uncertain_pathogenic_records = 52
    print_preproc("Uncertain", str(uncertain_neutral_records), str(uncertain_pathogenic_records))

    clean_training_file = os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'training.cbv.clean')
    clean_test_file     = os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'test.cbv.clean')
    data = np.loadtxt(clean_training_file, dtype='S20')
    clean_pathogenic_training_records = data[data[:, 4] == '1'].shape[0]
    clean_neutral_training_records    = data[data[:, 4] == '0'].shape[0]
    data = np.loadtxt(clean_test_file, dtype='S20')
    clean_pathogenic_test_records = data[data[:, 4] == '1'].shape[0]
    clean_neutral_test_records    = data[data[:, 4] == '0'].shape[0]
    print_preproc("Unknown reference",
                  str(original_neutral_records - (clean_neutral_test_records+clean_neutral_training_records+uncertain_neutral_records)),
                  str(original_pathogenic_records - (clean_pathogenic_test_records+clean_pathogenic_training_records+uncertain_pathogenic_records)),
                  )

    scores_training_file = os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'training.cbv.scores')
    scores_test_file     = os.path.join(combivep_settings.COMBIVEP_CENTRAL_TEST_CBV_DIR, 'test.cbv.scores')
    data = np.loadtxt(scores_training_file, dtype='S20')
    scores_pathogenic_training_records = data[data[:, 4] == '1'].shape[0]
    scores_neutral_training_records    = data[data[:, 4] == '0'].shape[0]
    data = np.loadtxt(scores_test_file, dtype='S20')
    scores_pathogenic_test_records = data[data[:, 4] == '1'].shape[0]
    scores_neutral_test_records    = data[data[:, 4] == '0'].shape[0]
    print_preproc("Unidentified by effect predictors",
                  str((clean_neutral_test_records+clean_neutral_training_records) - (scores_neutral_test_records+scores_neutral_training_records)),
                  str((clean_pathogenic_test_records+clean_pathogenic_training_records) - (scores_pathogenic_test_records+scores_pathogenic_training_records)),
                  )
    print_preproc("Ready to be used by CombiVEP & Condel",
                  str(scores_neutral_test_records+scores_neutral_training_records),
                  str(scores_pathogenic_test_records+scores_pathogenic_training_records),
                  )
    print_preproc("Training dataset",
                  str(scores_neutral_training_records),
                  str(scores_pathogenic_training_records),
                  )
    print_preproc("Test dataset",
                  str(scores_neutral_test_records),
                  str(scores_pathogenic_test_records),
                  )

def generate_precision_performance_report():
    discriminant_ratio = 0.5
    prediction_result  = np.loadtxt(devtools_settings.PUBLICATION_CONDEL_PREDICTION_RESULT, dtype='S20')

    cbv_positive_samples = prediction_result[prediction_result[:, 5].astype(np.float) > discriminant_ratio]
    cbv_true_positive_prediction  = cbv_positive_samples[cbv_positive_samples[:, 4] == '1'].shape[0]
    cbv_false_positive_prediction = cbv_positive_samples[cbv_positive_samples[:, 4] == '0'].shape[0]
    cbv_negative_samples = prediction_result[prediction_result[:, 5].astype(np.float) <= discriminant_ratio]
    cbv_true_negative_prediction  = cbv_negative_samples[cbv_negative_samples[:, 4] == '0'].shape[0]
    cbv_false_negative_prediction = cbv_negative_samples[cbv_negative_samples[:, 4] == '1'].shape[0]

    cbv_accuracy         = float(cbv_true_positive_prediction+cbv_true_negative_prediction)/(cbv_true_positive_prediction+cbv_true_negative_prediction+cbv_false_positive_prediction+cbv_false_negative_prediction)
    cbv_sensitivity      = float(cbv_true_positive_prediction)/(cbv_true_positive_prediction+cbv_false_negative_prediction)
    cbv_specificity      = float(cbv_true_negative_prediction)/(cbv_true_negative_prediction+cbv_false_positive_prediction)
    cbv_balance_accuracy = (cbv_sensitivity+cbv_specificity)/2

    cdl_positive_samples = prediction_result[prediction_result[:, 12].astype(np.float) > discriminant_ratio]
    cdl_true_positive_prediction  = cdl_positive_samples[cdl_positive_samples[:, 4] == '1'].shape[0]
    cdl_false_positive_prediction = cdl_positive_samples[cdl_positive_samples[:, 4] == '0'].shape[0]
    cdl_negative_samples = prediction_result[prediction_result[:, 12].astype(np.float) <= discriminant_ratio]
    cdl_true_negative_prediction  = cdl_negative_samples[cdl_negative_samples[:, 4] == '0'].shape[0]
    cdl_false_negative_prediction = cdl_negative_samples[cdl_negative_samples[:, 4] == '1'].shape[0]

    cdl_accuracy         = float(cdl_true_positive_prediction+cdl_true_negative_prediction)/(cdl_true_positive_prediction+cdl_true_negative_prediction+cdl_false_positive_prediction+cdl_false_negative_prediction)
    cdl_sensitivity      = float(cdl_true_positive_prediction)/(cdl_true_positive_prediction+cdl_false_negative_prediction)
    cdl_specificity      = float(cdl_true_negative_prediction)/(cdl_true_negative_prediction+cdl_false_positive_prediction)
    cdl_balance_accuracy = (cdl_sensitivity+cdl_specificity)/2

    print_precision("", "Condel", "CombiVEP")
    print_precision("False positive samples", str(cdl_false_positive_prediction), str(cbv_false_positive_prediction))
    print_precision("True negative samples",  str(cdl_true_negative_prediction),  str(cbv_true_negative_prediction))
    print_precision("False negative samples", str(cdl_false_negative_prediction), str(cbv_false_negative_prediction))
    print_precision("True positive samples",  str(cdl_true_positive_prediction),  str(cbv_true_positive_prediction))
    print_precision("Accuracy",               "%.4f" % cdl_accuracy,              "%.4f" % cbv_accuracy)
    print_precision("Sensitivity",            "%.4f" % cdl_sensitivity,           "%.4f" % cbv_sensitivity)
    print_precision("Specificity",            "%.4f" % cdl_specificity,           "%.4f" % cbv_specificity)
    print_precision("Balance accuracy",       "%.4f" % cdl_balance_accuracy,      "%.4f" % cbv_balance_accuracy)

