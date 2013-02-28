import os
import sys
import numpy as np
import combivep.devtools.settings as dev_const
import combivep.settings as cbv_const
from numpy import histogram as hist
from combivep.preproc.dataset import DataSetManager
from combivep.preproc.dataset import SnpDataRecord
from combivep.preproc.reader import ScoresRecord
from combivep.config import Configure
from combivep.template import CombiVEPBase
from combivep.engine.wrapper import Trainer
from combivep.engine.wrapper import Predictor
from collections import namedtuple
from combivep.devtools.settings import PRECISION_MEASURES
from sklearn.metrics import auc


PrecisionPerformance = namedtuple("precision_performance", PRECISION_MEASURES)


class DevScoresRecord(CombiVEPBase):
    """ to automatically parse data with scores """

    def __init__(self, array_snp):
        self.__array_snp = array_snp

    @property
    def chrom(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_CHROM]

    @property
    def pos(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_POS]

    @property
    def ref(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_REF]

    @property
    def alt(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_ALT]

    @property
    def target(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_TARGET]

    @property
    def phylop_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_PHYLOP_SCORE]

    @property
    def sift_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_SIFT_SCORE]

    @property
    def pp2_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_PP2_SCORE]

    @property
    def lrt_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_LRT_SCORE]

    @property
    def mt_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_MT_SCORE]

    @property
    def gerp_score(self):
        return self.__array_snp[dev_const.SCORES_0_IDX_GERP_SCORE]


class ScoresReader(CombiVEPBase):
    """

    to read parsed SNPS file
    The format are CHROM, POS, REF, ALT, EFFECT.
    All fields are tab separated.

    """

    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, scores_file):
        self.scores_file_name = scores_file

    def fetch_snps(self):
        scores_file = open(self.scores_file_name)
        for line in scores_file:
            if line[0] == '#':
                continue
            yield DevScoresRecord(line.rstrip('\n').split('\t'))


class FastDataSetManager(DataSetManager):

    def __init__(self, cfg_file=cbv_const.CBV_CFG_FILE):
        DataSetManager.__init__(self, cfg_file=cbv_const.CBV_CFG_FILE)

    def load_data(self, file_name, file_type=dev_const.FILE_TYPE_SCORES):
        if file_type == dev_const.FILE_TYPE_SCORES:
            return self.__load_scores_data(file_name)
        else:
            return DataSetManager.load_data(self, file_name, file_type)

    def __load_scores_data(self, file_name):
        self.clear_data()
        scores_reader = ScoresReader()
        scores_reader.read(file_name)
        for rec in scores_reader.fetch_snps():
            snp_data = SnpDataRecord(rec)
            scores   = ScoresRecord(rec)
            self.dataset.append({cbv_const.KW_SNP_DATA: snp_data,
                                 cbv_const.KW_SCORES: scores,
                                 })


def filter_cbv_data(cbv_file,
                    cfg_file=cbv_const.CBV_CFG_FILE):
    (dir_name, file_name) = os.path.split(cbv_file)
    report_fmt = "{caption:<25}:{value:>6d}"
    clean_out_fmt = "{chrom}\t{pos}\t{ref}\t{alt}\t{target}\n"
    scores_out_fmt = "{chrom}\t{pos}\t{ref}\t{alt}\t{target}\t"
    scores_out_fmt += "{phylop_score}\t{sift_score}\t{pp2_score}\t"
    scores_out_fmt += "{lrt_score}\t{mt_score}\t{gerp_score}\n"

    print
    print "> > >  " + file_name
    dm = DataSetManager()
    dm.load_data(cbv_file, file_type=cbv_const.FILE_TYPE_CBV)
    print report_fmt.format(caption="original",
                            value=len(dm.dataset))
    print

    dm.validate_data()
    f_clean = open(cbv_file + '.clean', 'w')
    for item in dm.dataset:
        snp_data = item[cbv_const.KW_SNP_DATA]
        f_clean.write(clean_out_fmt.format(chrom=snp_data.chrom,
                                           pos=snp_data.pos,
                                           ref=snp_data.ref,
                                           alt=snp_data.alt,
                                           target=snp_data.target,
                                           ))
    f_clean.close()
    print report_fmt.format(caption="Clean pathogenic",
                            value=len([item for item in dm.dataset if item[cbv_const.KW_SNP_DATA].target == '1']))
    print report_fmt.format(caption="Clean neutral",
                            value=len([item for item in dm.dataset if item[cbv_const.KW_SNP_DATA].target == '0']))
    print report_fmt.format(caption="Total",
                            value=len(dm.dataset))
    print

    dm.calculate_scores()
    f_scores = open(cbv_file + '.scores', 'w')
    for item in dm.dataset:
        snp_data = item[cbv_const.KW_SNP_DATA]
        scores   = item[cbv_const.KW_SCORES]
        f_scores.write(scores_out_fmt.format(chrom=snp_data.chrom,
                                             pos=snp_data.pos,
                                             ref=snp_data.ref,
                                             alt=snp_data.alt,
                                             target=snp_data.target,
                                             phylop_score=scores.phylop_score,
                                             sift_score=scores.sift_score,
                                             pp2_score=scores.pp2_score,
                                             lrt_score=scores.lrt_score,
                                             mt_score=scores.mt_score,
                                             gerp_score=scores.gerp_score
                                             ))
    f_scores.close()
    print report_fmt.format(caption="Scored pathogenic",
                            value=len([item for item in dm.dataset if item[cbv_const.KW_SNP_DATA].target == '1']))
    print report_fmt.format(caption="Scored neutral",
                            value=len([item for item in dm.dataset if item[cbv_const.KW_SNP_DATA].target == '0']))
    print report_fmt.format(caption="Total",
                            value=len(dm.dataset))
    print

    dm.set_shuffle_seed(cbv_const.DEMO_SEED)
    dm.shuffle_data()
    dm.partition_data()

    #partition data
    training_dataset   = dm.get_training_data()
    validation_dataset = dm.get_validation_data()

    print report_fmt.format(caption="Training pathogenic",
                            value=len([item for item in training_dataset if item[cbv_const.KW_SNP_DATA].target == '1']))
    print report_fmt.format(caption="Training neutral",
                            value=len([item for item in training_dataset if item[cbv_const.KW_SNP_DATA].target == '0']))
    print report_fmt.format(caption="Total",
                            value=len(training_dataset))
    print

    print report_fmt.format(caption="Validation pathogenic",
                            value=len([item for item in validation_dataset if item[cbv_const.KW_SNP_DATA].target == '1']))
    print report_fmt.format(caption="Validation neutral",
                            value=len([item for item in validation_dataset if item[cbv_const.KW_SNP_DATA].target == '0']))
    print report_fmt.format(caption="Total",
                            value=len(validation_dataset))


def calculate_roc(patho_dataset, neutr_dataset, roc_range):
    fp_rates = np.zeros([len(roc_range),
                         patho_dataset.shape[1]])
    tp_rates = np.zeros([len(roc_range),
                         neutr_dataset.shape[1]])
    patho_data_size = patho_dataset.shape[0]
    neutr_data_size = neutr_dataset.shape[0]
    for i in xrange(len(roc_range)):
        fp_rates[i, :] = np.matrix(np.sum(neutr_dataset > roc_range[i],
                                          axis=0
                                          ).astype(np.float))/neutr_data_size
        tp_rates[i, :] = np.matrix(np.sum(patho_dataset > roc_range[i],
                                          axis=0
                                          ).astype(np.float))/patho_data_size
    return (fp_rates, tp_rates)


def print_preproc(col1, col2, col3, col4=None):
    report_fmt = "{measurement:<40}{neutr_val:>18}{patho_val:>18}{total_val:>18}"
    if col4 is None:
        print report_fmt.format(measurement=col1,
                                neutr_val=col2,
                                patho_val=col3,
                                total_val=int(col2)+int(col3))
    else:
        print report_fmt.format(measurement=col1,
                                neutr_val=col2,
                                patho_val=col3,
                                total_val=col4)


def info(msg):
    print >> sys.stderr, msg


def fast_training(training_data_file,
                  params_out_file=cbv_const.USER_PARAMS_FILE,
                  random_seed=cbv_const.DFLT_SEED,
                  n_hidden_nodes=cbv_const.DFLT_HIDDEN_NODES,
                  figure_dir=cbv_const.DFLT_FIGURE_DIR,
                  iterations=cbv_const.DFLT_ITERATIONS,
                  cfg_file=cbv_const.CBV_CFG_FILE,
                  ):
    """

    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
    CBV has 5 fields, CHROM, POS, REF, ALT, EFFECT (1=deleterious, 0=neutr).
    All are tab separated
    Required arguments
    - neutr_data_file : list of SNPs with no harmful effect, CBV format
    - pathognice_data_file : list of SNPs with deleterious effect, CBV format

    """
    #pre-processing dataset
    info('pre-processing dataset, this may take a while (around 750 SNPs/mins). . .')
    dm = FastDataSetManager(cfg_file=cfg_file)
    dm.load_data(training_data_file, file_type=dev_const.FILE_TYPE_SCORES)
    dm.set_shuffle_seed(random_seed)
    dm.shuffle_data()
    dm.partition_data()

    #partition data
    training_dataset   = dm.get_training_data()
    validation_dataset = dm.get_validation_data()

    #train !!!
    info('Training CombiVEP, please wait (around 500 SNPs/mins) . . .')
    trainer = Trainer(training_dataset,
                      validation_dataset,
                      random_seed,
                      n_hidden_nodes,
                      figure_dir)
    trainer.train(iterations)
    if not os.path.exists(cbv_const.USER_PARAMS_DIR):
        os.makedirs(cbv_const.USER_PARAMS_DIR)
    trainer.export_best_parameters(params_out_file)


def fast_predict(SNPs_file,
                 params_file=cbv_const.USER_PARAMS_FILE,
                 file_type=cbv_const.FILE_TYPE_VCF,
                 output_file=None,
                 cfg_file=cbv_const.CBV_CFG_FILE,
                 ):
    """

    CBV (CombiVEP format) is a parsed format intended to be used by CombiVEP.
    CBV has 5 fields, CHROM, POS, REF, ALT, EFFECT (1=deleterious, 0=neutr).
    All are tab separated
    Required arguments
    - SNPs_file : list of SNPs to be predicted, can be either VCF or CBV
                  (default is VCF)

    """
    #pre-processing test dataset
    info('pre-processing dataset, this may take a while (around 750 SNPs/mins). . .')
    dm = FastDataSetManager(cfg_file=cfg_file)
    dm.load_data(SNPs_file, file_type=dev_const.FILE_TYPE_SCORES)

    #predict
    predictor = Predictor()
    predictor.import_parameters(params_file=params_file)
    out = (np.array(predictor.predict(dm.dataset)).reshape(-1,))

    #print output
    if output_file is not None:
        sys.stdout = open(output_file, 'w')
    tmp_rec = []
    tmp_rec.append("CHROM")
    tmp_rec.append("POS")
    tmp_rec.append("REF")
    tmp_rec.append("ALT")
    tmp_rec.append("ACTUAL_DELETERIOUS_EFFECT")
    tmp_rec.append("PREDICTED_DELETERIOUS_PROBABILITY")
    tmp_rec.append("PHYLOP_SCORE")
    tmp_rec.append("SIFT_SCORE")
    tmp_rec.append("PP2_SCORE")
    tmp_rec.append("LRT_SCORT")
    tmp_rec.append("MT_SCORE")
    tmp_rec.append("GERP_SCORE")
    print "#" + "\t".join(tmp_rec)
    for i in xrange(len(dm.dataset)):
        del tmp_rec[:]
        snp_data = dm.dataset[i][cbv_const.KW_SNP_DATA]
        scores   = dm.dataset[i][cbv_const.KW_SCORES]
        tmp_rec.append(snp_data.chrom)
        tmp_rec.append(snp_data.pos)
        tmp_rec.append(snp_data.ref)
        tmp_rec.append(snp_data.alt)
        tmp_rec.append(snp_data.target)
        tmp_rec.append("%6.4f" % out[i])
        tmp_rec.append(scores.phylop_score)
        tmp_rec.append(scores.sift_score)
        tmp_rec.append(scores.pp2_score)
        tmp_rec.append(scores.lrt_score)
        tmp_rec.append(scores.mt_score)
        tmp_rec.append(scores.gerp_score)
        print "\t".join(tmp_rec)
    sys.stdout = sys.__stdout__


def generate_roc_auc_figures(plt,
                             patho_scores,
                             neutr_scores,
                             predictor_names,
                             predictor_colors):
    min_value = min(np.amin(patho_scores), np.amin(neutr_scores))
    max_value = max(np.amax(patho_scores), np.amax(neutr_scores))

    #produce roc data from CombiVEP, Phylop, SIFT, PP2, LRT, MT, GERP, Condel
    fp_rates, tp_rates = calculate_roc(patho_scores,
                                       neutr_scores,
                                       np.linspace(min_value, max_value, 5001))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in xrange(len(predictor_names)):
        ax.plot(fp_rates[:, i],
                tp_rates[:, i],
                predictor_colors[i],
                label=predictor_names[i])
    ax.set_ylabel('true positive rate')
    ax.set_xlabel('false positive rate')
    ax.legend(bbox_to_anchor=(0.99, 0.01), loc=4)
    fig.savefig(dev_const.PUB_ROC_FIG, bbox_inches='tight', pad_inches=0.05)

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
        ax.text(ind[i], aucs[i] + 0.005, "%0.4f" % aucs[i], size=10)
    ax.set_ylim([0.7, 0.9])
    ax.set_xticks(np.array(ind) + 0.15)
    ax.set_xticklabels(predictor_names, rotation=30)
    ax.set_ylabel('Area Under ROC Curve')
    ax.set_xlabel('Predictors')
    fig.savefig(dev_const.PUB_AUC_FIG, bbox_inches='tight', pad_inches=0.05)
    return dev_const.PUB_ROC_FIG, dev_const.PUB_AUC_FIG


def scores_dist_plot(figure,
                     subplot,
                     patho_scores,
                     neutr_scores,
                     title):
    #calculate histogram
    hist_range = (-0.005, 1.005)
    patho_hist, bins = hist(patho_scores,
                            bins=100,
                            range=hist_range)
    neutr_hist, bins = hist(neutr_scores,
                            bins=100,
                            range=hist_range)
    center = (bins[:-1]+bins[1:]) / 2
    #add the first and the last data point to smooth the plot
    center = np.concatenate((np.array([-0.005]), center), axis=0)
    center = np.concatenate((center, np.array([1.005])), axis=0)
    patho_hist = np.concatenate((np.array([0]), patho_hist), axis=0)
    patho_hist = np.concatenate((patho_hist, np.array([0])), axis=0)
    neutr_hist = np.concatenate((np.array([0]), neutr_hist), axis=0)
    neutr_hist = np.concatenate((neutr_hist, np.array([0])), axis=0)
    ax = figure.add_subplot(subplot)
    ax.plot(center, patho_hist, 'r--', label='pathogenic variants')
    ax.plot(center, neutr_hist, 'b--', label='neutral variants')
    ax.set_title(title)
    ax.set_ylim([0, 450])
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylabel('samples')
    ax.set_xlabel('score')
    return ax


def generate_scores_dist_figure(plt,
                                patho_scores,
                                neutr_scores):
    fig = plt.figure()
    #CombiVEP
    ax = scores_dist_plot(fig,
                          221,
                          patho_scores[:, 0],
                          neutr_scores[:, 0],
                          'CombiVEP scores distribution')
    ax.legend(bbox_to_anchor=(1.435, 1.35), loc=1)
    #Condel
    ax = scores_dist_plot(fig,
                          222,
                          patho_scores[:, 7],
                          neutr_scores[:, 7],
                          'Condel scores distribution')
    #PP2
    ax = scores_dist_plot(fig,
                          223,
                          patho_scores[:, 3],
                          neutr_scores[:, 3],
                          'PolyPhen2 scores distribution')
    #SIFT
    ax = scores_dist_plot(fig,
                          224,
                          patho_scores[:, 2],
                          neutr_scores[:, 2],
                          'SIFT scores distribution')
    fig.tight_layout()
    fig.savefig(dev_const.PUB_SCORES_DISTR_FIG,
                bbox_inches='tight',
                pad_inches=0.05)
    return dev_const.PUB_SCORES_DISTR_FIG


def precision_int_plot(ax,
                       ind,
                       results,
                       width,
                       color):
    bars = ax.bar(ind,
                  results,
                  width,
                  color=color)
    for i in xrange(len(results)):
        ax.text(ind[i]+(width/2),
                results[i]+5,
                results[i],
                size=8,
                ha='center')
    return bars


def precision_float_plot(ax,
                         ind,
                         results,
                         width,
                         color):
    bars = ax.bar(ind,
                  results,
                  width,
                  color=color)
    for i in xrange(len(results)):
        ax.text(ind[i]+(width/2)-0.02,
                results[i]+0.01,
                "%.2f" % results[i],
                size=8,
                ha='center')
    return bars


def generate_precision_measurements_figures(plt,
                                            patho_scores,
                                            neutr_scores,
                                            predictor_names,
                                            predictor_colors):
    #measure precision by using 0.5 as a discriminant retio
    precision_performances = []
    for i in xrange(patho_scores.shape[1]):
        precision_performances.append(measure_precision(0.5,
                                                        patho_scores[:, i],
                                                        neutr_scores[:, i]))

    width = 0.1
    ind = np.subtract(np.multiply(np.add(np.arange(len(predictor_names)),
                                         1),
                                  0.5),
                      0.45)
    #plot number of results in each class for each predictor
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bars1 = precision_int_plot(ax,
                               ind,
                               map(lambda x: x.true_positive,
                                   precision_performances),
                               width,
                               'g')
    bars2 = precision_int_plot(ax,
                               ind+width,
                               map(lambda x: x.false_positive,
                                   precision_performances),
                               width,
                               'darkred')
    bars3 = precision_int_plot(ax,
                               ind+(2*width),
                               map(lambda x: x.true_negative,
                                   precision_performances),
                               width,
                               'b')
    bars4 = precision_int_plot(ax,
                               ind+(3*width),
                               map(lambda x: x.false_negative,
                                   precision_performances),
                               width,
                               'r')
    ax.set_ylim([0, 1400])
    ax.set_xticks(ind + 0.2)
    ax.set_xticklabels(predictor_names, rotation=30)
    ax.set_ylabel('Precision performance')
    ax.set_xlabel('Predictors')
    ax.legend((bars1[0], bars2[0], bars3[0], bars4[0]),
              ('True positive results',
               'False positive results',
               'True negative results',
               'False negative results'),
              bbox_to_anchor=(0.99, 0.99),
              loc=1)
    fig.savefig(dev_const.PUB_RESULT_CLASSES_FIG,
                bbox_inches='tight',
                pad_inches=0.05)
    #plot precision measurement
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bars1 = precision_float_plot(ax,
                                 ind,
                                 map(lambda x: x.accuracy,
                                     precision_performances),
                                 width,
                                 'g')
    bars2 = precision_float_plot(ax,
                                 ind+width,
                                 map(lambda x: x.sensitivity,
                                     precision_performances),
                                 width,
                                 'darkred')
    bars3 = precision_float_plot(ax,
                                 ind+(2*width),
                                 map(lambda x: x.specificity,
                                     precision_performances),
                                 width,
                                 'b')
    bars4 = precision_float_plot(ax,
                                 ind+(3*width),
                                 map(lambda x: x.balance_accuracy,
                                     precision_performances),
                                 width,
                                 'r')
    ax.set_ylim([0, 1.19])
    ax.set_xticks(ind + 0.2)
    ax.set_xticklabels(predictor_names, rotation=30)
    ax.set_ylabel('Precision performance')
    ax.set_xlabel('Predictors')
    ax.legend((bars1[0], bars2[0], bars3[0], bars4[0]),
              ('Accuracy',
               'Sensitivity',
               'Specificity',
               'Balance accuracy'),
              bbox_to_anchor=(0.9999, 0.9999),
              loc=1)
    fig.savefig(dev_const.PUB_PRECISION_FIG,
                bbox_inches='tight',
                pad_inches=0.05)
    return dev_const.PUB_RESULT_CLASSES_FIG, dev_const.PUB_PRECISION_FIG


def measure_precision(ratio,
                      pos_samples_scores,
                      neg_samples_scores,
                      ):
    true_pos  = len(pos_samples_scores[pos_samples_scores > ratio])
    false_neg = len(pos_samples_scores[pos_samples_scores <= ratio])
    true_neg  = len(neg_samples_scores[neg_samples_scores < ratio])
    false_pos = len(neg_samples_scores[neg_samples_scores >= ratio])

    pos_prediction = true_pos + false_pos
    neg_prediction = true_neg + false_neg

    total_samples = true_pos + true_neg + false_pos + false_neg

    accuracy = float(true_pos+true_neg) / total_samples
    sensitivity = float(true_pos) / (true_pos+false_neg)
    specificity = float(true_neg) / (true_neg+false_pos)
    balance_accuracy = (sensitivity+specificity) / 2

    return PrecisionPerformance(pos_prediction,
                                true_pos,
                                false_pos,
                                neg_prediction,
                                true_neg,
                                false_neg,
                                accuracy,
                                sensitivity,
                                specificity,
                                balance_accuracy)
