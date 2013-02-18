import os

"""

This module is all about constant value that are used in this application

"""

DEBUG_MODE = 0

# > > > > > > > > > > > > > development files & folders < < < < < < < < < <
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))

#'central' data for testing and for demo
CBV_SAMPLE_DATA_ROOT   = os.path.join(PROJECT_ROOT, 'combivep/data')
CBV_SAMPLE_DATASET_DIR = os.path.join(CBV_SAMPLE_DATA_ROOT, 'dataset')
CBV_SAMPLE_CBV_DIR     = os.path.join(CBV_SAMPLE_DATA_ROOT, 'CBV')
CBV_SAMPLE_VCF_DIR     = os.path.join(CBV_SAMPLE_DATA_ROOT, 'VCF')
CBV_SAMPLE_UCSC_DIR    = os.path.join(CBV_SAMPLE_DATA_ROOT, 'UCSC')
CBV_SAMPLE_LJB_DIR     = os.path.join(CBV_SAMPLE_DATA_ROOT, 'LJB')
CBV_SAMPLE_CFG_FILE    = os.path.join(CBV_SAMPLE_DATA_ROOT, 'config.txt')
CBV_SAMPLE_PARAM_DIR   = os.path.join(CBV_SAMPLE_DATA_ROOT, 'params')
CBV_SAMPLE_PARAM_FILE  = os.path.join(CBV_SAMPLE_PARAM_DIR, 'params.npz')


# > > > > > > > > > > > > > User files & folders < < < < < < < < < <
USER_DATA_ROOT = os.path.expanduser('~/.CombiVEP')

#to keep user data produced by CombiVEP engine
USER_PARAMS_DIR  = os.path.join(USER_DATA_ROOT, 'params')
USER_PARAMS_FILE = os.path.join(USER_PARAMS_DIR, 'params.npz')

#to keep the reference database from UCSC and LJB
USER_UCSC_REF_DB_DIR = os.path.join(USER_DATA_ROOT, 'ref/UCSC')
USER_LJB_REF_DB_DIR  = os.path.join(USER_DATA_ROOT, 'ref/LJB')


# > > > > > > > > > > > > > status file < < < < < < < < < <
CBV_CFG_FILE = os.path.join(USER_DATA_ROOT, 'config.txt')

#key
LATEST_UCSC_DB_VER     = 'latest_ucsc_database_version'
LATEST_UCSC_FILE_NAME  = 'latest_ucsc_file_name'
LATEST_LJB_DB_VER      = 'latest_ljb_database_version'
LATEST_LJB_FILE_PREFIX = 'latest_ljb_file_prefix'


# > > > > > > > > > > > > > temporay files and folder < < < < < < < < < <
#the only temporay working folder used for processing data
CBV_WORKING_DIR = os.path.join(USER_DATA_ROOT, 'tmp')

#temporary files for reference database processing
TMP_UCSC_CLEAN_DB_FILE = os.path.join(CBV_WORKING_DIR, 'tmp_ucsc_clean_db.txt')
TMP_LJB_CLEAN_DB_FILE  = os.path.join(CBV_WORKING_DIR, 'tmp_ljb_clean_db.txt')


# > > > > > > > > > > > > > MLP configuration < < < < < < < < < <
#general model configuration
MLP_COEFFICIENT   = 0.9
STEP_SIZE         = 0.0005
MAX_ALLOWED_ERROR = 0.35
MIN_IMPROVEMENT   = 0.000000001

#default model argument values
DFLT_ITERATIONS   = 30000
DFLT_HIDDEN_NODES = 6
DFLT_SEED         = None
DFLT_FIGURE_DIR   = None

#proportion of data partitioning
PROPORTION_TRAINING_DATA   = 70
PROPORTION_VALIDATION_DATA = 15


# > > > > > > > > > > > > > Demo configuration < < < < < < < < < <
DEMO_SEED = 20


# > > > > > > > reference database configuration < < < < < < <
#UCSC
UCSC_FOLDER_URL      = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database'
UCSC_LIST_FILE_NAME  = 'ucsc_list_file'
UCSC_FILES_PATTERN   = r"""href="(?P<file_name>snp\d{3}.txt.gz)">.*>.*(?P<date>\d{2}-[a-zA-Z]{3}-\d{4})"""
UCSC_VERSION_PATTERN = r"""[a-zA-Z]*(?P<version>[\d]*)[a-zA-Z.]*"""

#LJB
LJB_FOLDER_URL      = 'http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/'
LJB_LIST_FILE_NAME  = 'ljb_list_file'
LJB_FILES_PATTERN   = r"""href="(?P<file_name>dbNSFP_light.*.zip)">"""
LJB_VERSION_PATTERN = r"""[a-zA-Z_]*(?P<version>[\d.]*)[.][a-zA-Z.]*"""


# > > > > > > > > > > > > > Dataset data structure < < < < < < < < < <
#section key
KW_SNP_DATA = 'snp_data'
KW_SCORES   = 'scores'

##global SNP data keywords
#KW_CHROM = 'chrom'
#KW_POS   = 'pos'
#KW_REF   = 'ref'
#KW_ALT   = 'alt'
#KW_TARGET = 'target'

#file type
FILE_TYPE_CBV = 'CBV'
FILE_TYPE_VCF = 'VCF'


# > > > > > > > > > > > > > UCSC format configuration < < < < < < < < < <
##general key
#KW_UCSC_CHROM     = 'ucsc_chrom'
#KW_UCSC_START_POS = 'ucsc_start_pos'
#KW_UCSC_END_POS   = 'ucsc_end_pos'
#KW_UCSC_STRAND    = 'ucsc_strand'
#KW_UCSC_REF       = 'ucsc_ref'
#KW_UCSC_OBSERVED  = 'ucsc_observed'

#UCSC index
#0-based index, used by python
UCSC_0_IDX_CHROM     = 1
UCSC_0_IDX_START_POS = 2
UCSC_0_IDX_END_POS   = 3
UCSC_0_IDX_STRAND    = 6
UCSC_0_IDX_REF       = 8
UCSC_0_IDX_OBSERVED  = 9
UCSC_EXPECTED_LENGTH = 26


# > > > > > > > > > > > > > LJB format configuration < < < < < < < < < <
##SNP information key
#KW_LJB_CHROM = 'ljb_chrom'
#KW_LJB_POS   = 'ljb_hg19_pos'
#KW_LJB_REF   = 'ljb_ref'
#KW_LJB_ALT   = 'ljb_alt'

##score key
#KW_PHYLOP_SCORE = 'phylop_score'
#KW_SIFT_SCORE   = 'sift_score'
#KW_PP2_SCORE    = 'pp2_score'
#KW_LRT_SCORE    = 'lrt_score'
#KW_MT_SCORE     = 'mt_score'
#KW_GERP_SCORE   = 'gerp_score'

#LJB index
#1-based index, used by awk for parsing
LJB_RAW_1_IDX_CHROM        = 1
LJB_RAW_1_IDX_HG19_POS     = 7
LJB_RAW_1_IDX_REF          = 3
LJB_RAW_1_IDX_ALT          = 4
LJB_RAW_1_IDX_PHYLOP_SCORE = 8
LJB_RAW_1_IDX_SIFT_SCORE   = 9
LJB_RAW_1_IDX_PP2_SCORE    = 10
LJB_RAW_1_IDX_LRT_SCORE    = 11
LJB_RAW_1_IDX_MT_SCORE     = 13
LJB_RAW_1_IDX_GERP_SCORE   = 17

#0-based index, used by python for reading
LJB_PARSED_0_IDX_CHROM        = 0
LJB_PARSED_0_IDX_POS          = 1
LJB_PARSED_0_IDX_REF          = 2
LJB_PARSED_0_IDX_ALT          = 3
LJB_PARSED_0_IDX_PHYLOP_SCORE = 4
LJB_PARSED_0_IDX_SIFT_SCORE   = 5
LJB_PARSED_0_IDX_PP2_SCORE    = 6
LJB_PARSED_0_IDX_LRT_SCORE    = 7
LJB_PARSED_0_IDX_MT_SCORE     = 8
LJB_PARSED_0_IDX_GERP_SCORE   = 9
LJB_PARSED_EXPECTED_LENGTH    = 10


# > > > > > > > > > > > > > VCF format configuration < < < < < < < < < <
##SNP information key
#KW_VCF_CHROM = 'vcf_chrom'
#KW_VCF_POS   = 'vcf_pos'
#KW_VCF_REF   = 'vcf_ref'
#KW_VCF_ALT   = 'vcf_alt'

#VCF index
#0-based index, used by python
VCF_0_IDX_CHROM = 0
VCF_0_IDX_POS   = 1
VCF_0_IDX_REF   = 3
VCF_0_IDX_ALT   = 4


# > > > > > > > > > > CBV format (CombiVEP format) configuration < < < < < < <
##SNP information key
#KW_CBV_CHROM  = 'CBV_chrom'
#KW_CBV_POS    = 'CBV_pos'
#KW_CBV_REF    = 'CBV_ref'
#KW_CBV_ALT    = 'CBV_alt'
#KW_CBV_TARGET = 'CBV_target'

#CBV index
#0-based index, used by python
CBV_0_IDX_CHROM  = 0
CBV_0_IDX_POS    = 1
CBV_0_IDX_REF    = 2
CBV_0_IDX_ALT    = 3
CBV_0_IDX_TARGET = 4


# > > > > > > > > > > prediction output format < < < < < < < < < <
PREDICTION_OUT_COLS_HEADER = ['CHROM',
                              'POS',
                              'REF',
                              'ALT',
                              'ACTUAL_DELETERIOUS_EFFECT',
                              'PREDICTED_DELETERIOUS_PROBABILITY',
                              'PHYLOP_SCORE',
                              'SIFT_SCORE',
                              'PP2_SCORE',
                              'LRT_SCORT',
                              'MT_SCORE',
                              'GERP_SCORE']
