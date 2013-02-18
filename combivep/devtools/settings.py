import os


# > > > > > > > > > > > > > development files & folders < < < < < < < < < <
PUB_RES_ROOT = '/home/jessada/development/scilifelab/assignments/20121119_CombiVEP_publication/resources'

#general resources
PUB_RAW_PREDICTION_RESULT    = os.path.join(PUB_RES_ROOT,
                                            'prediction_result')
PUB_CONDEL_PREDICTION_RESULT = os.path.join(PUB_RES_ROOT,
                                            'full_prediction_result')
PUB_PARAM_FILE = os.path.join(PUB_RES_ROOT, 'params.npz')

#figures
PUB_FIGS_DIR = os.path.join(PUB_RES_ROOT, 'figures')
PUB_ROC_FIG  = os.path.join(PUB_FIGS_DIR, 'roc_curve.eps')
PUB_AUC_FIG  = os.path.join(PUB_FIGS_DIR, 'auc_curve.eps')
PUB_SCORES_DISTR_FIG = os.path.join(PUB_FIGS_DIR, 'scores_dist.eps')


# > > > > > > > > > > > > > file type < < < < < < < < < <
#file type
FILE_TYPE_SCORES = 'scores'

#SNP information key
KEY_SCORES_CHROM  = 'SCORES_chrom'
KEY_SCORES_POS    = 'SCORES_pos'
KEY_SCORES_REF    = 'SCORES_ref'
KEY_SCORES_ALT    = 'SCORES_alt'
KEY_SCORES_TARGET = 'SCORES_target'

#SCORES index
#0-based index, used by python
SCORES_0_IDX_CHROM        = 0
SCORES_0_IDX_POS          = 1
SCORES_0_IDX_REF          = 2
SCORES_0_IDX_ALT          = 3
SCORES_0_IDX_TARGET       = 4
SCORES_0_IDX_PHYLOP_SCORE = 5
SCORES_0_IDX_SIFT_SCORE   = 6
SCORES_0_IDX_PP2_SCORE    = 7
SCORES_0_IDX_LRT_SCORE    = 8
SCORES_0_IDX_MT_SCORE     = 9
SCORES_0_IDX_GERP_SCORE   = 10

# > > > > > > > > > > > > > data structure < < < < < < < < < <
PREDICTOR_NAMES = ('CombiVEP',
                   'Phylop',
                   'SIFT',
                   'PP2',
                   'LRT',
                   'MT',
                   'GERP',
                   'Condel')
PRECISION_MEASURES = ('true_postive',
                      'false_negative',
                      'true_negative',
                      'false_positive',
                      'accuracy',
                      'sensitivity',
                      'specificity',
                      'balance_accuracy')
