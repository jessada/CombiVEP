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






