import os
import combivep.settings as combivep_settings


# > > > > > > > > > > > > > development files & folders < < < < < < < < < <
PUBLICATION_RESOURCES_ROOT             = '/home/jessada/development/scilifelab/assignments/20121119_CombiVEP_publication/resources'
PUBLICATION_PARAMETER_FILE             = os.path.join(PUBLICATION_RESOURCES_ROOT, 'params.npz')
PUBLICATION_RAW_PREDICTION_RESULT      = os.path.join(PUBLICATION_RESOURCES_ROOT, 'prediction_result')
PUBLICATION_CONDEL_PREDICTION_RESULT   = os.path.join(PUBLICATION_RESOURCES_ROOT, 'full_prediction_result')
PUBLICATION_FIGURES_DIR                = os.path.join(PUBLICATION_RESOURCES_ROOT, 'figures')
PUBLICATION_SCORES_DISTRIBUTION_FIGURE = os.path.join(PUBLICATION_FIGURES_DIR, 'scores_dist.eps')
PUBLICATION_ROC_FIGURE                 = os.path.join(PUBLICATION_FIGURES_DIR, 'roc_curve.eps')
PUBLICATION_AUC_FIGURE                 = os.path.join(PUBLICATION_FIGURES_DIR, 'auc_curve.eps')
PUBLICATION_SCORE_DISTRIBUTION_FIGURE  = os.path.join(PUBLICATION_FIGURES_DIR, 'score_distribution.eps')


# > > > > > > > > > > > > > file type < < < < < < < < < <
#file type
FILE_TYPE_SCORES = 'scores' #CombiVEP format

#SNP information key
KEY_SCORES_CHROM   = 'SCORES_chrom'
KEY_SCORES_POS     = 'SCORES_pos'
KEY_SCORES_REF     = 'SCORES_ref'
KEY_SCORES_ALT     = 'SCORES_alt'
KEY_SCORES_TARGETS = 'SCORES_targets'

#SCORES index
#0-based index, used by python
SCORES_0_INDEX_CHROM        = 0
SCORES_0_INDEX_POS          = 1
SCORES_0_INDEX_REF          = 2
SCORES_0_INDEX_ALT          = 3
SCORES_0_INDEX_TARGETS      = 4
SCORES_0_INDEX_PHYLOP_SCORE = 5
SCORES_0_INDEX_SIFT_SCORE   = 6
SCORES_0_INDEX_PP2_SCORE    = 7
SCORES_0_INDEX_LRT_SCORE    = 8
SCORES_0_INDEX_MT_SCORE     = 9
SCORES_0_INDEX_GERP_SCORE   = 10



