import pysam
import combivep.settings as cbv_const
from combivep.template import CombiVEPBase


class UcscReader(CombiVEPBase):
    """to read UCSC parsed file in tabix format"""

    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, tabix_file):
        self.db_file_name = tabix_file

    def fetch_array_snps(self, chromosome, start_pos, end_pos):
        tabix_file = pysam.Tabixfile(self.db_file_name)
        if chromosome.startswith('chr'):
            chrom = chromosome
        else:
            chrom = 'chr' + chromosome
        for line in tabix_file.fetch(chrom, start_pos, end_pos):
            yield line.rstrip('\n').split('\t')

    def fetch_hash_snps(self, chromosome, start_pos, end_pos):
        for rec in self.fetch_array_snps(chromosome, start_pos, end_pos):
            if len(rec) != cbv_const.UCSC_EXPECTED_LENGTH :
                raise Exception("Invalid formatting is found in file '%s'>> Chrom : %s\tStart pos : %s\tEnd pos : %s" % (self.db_file_name, rec[cbv_const.UCSC_0_IDX_CHROM], rec[cbv_const.UCSC_0_IDX_START_POS], rec[cbv_const.UCSC_0_IDX_END_POS]))
            snp_info = {cbv_const.KEY_UCSC_CHROM     : rec[cbv_const.UCSC_0_IDX_CHROM], 
                        cbv_const.KEY_UCSC_START_POS : rec[cbv_const.UCSC_0_IDX_START_POS],
                        cbv_const.KEY_UCSC_END_POS   : rec[cbv_const.UCSC_0_IDX_END_POS],
                        cbv_const.KEY_UCSC_STRAND    : rec[cbv_const.UCSC_0_IDX_STRAND],
                        cbv_const.KEY_UCSC_REF       : rec[cbv_const.UCSC_0_IDX_REF],
                        cbv_const.KEY_UCSC_OBSERVED  : rec[cbv_const.UCSC_0_IDX_OBSERVED],
                        }
            yield {cbv_const.KEY_SNP_INFO_SECTION : snp_info}


class LjbReader(CombiVEPBase):
    """to read parsed LJB file"""


    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, ljb_file):
        self.db_file_name = ljb_file
        self.__tabix_file = pysam.Tabixfile(ljb_file)

    def fetch_array_snps(self, chromosome, start_pos, end_pos):
        for line in self.__tabix_file.fetch(chromosome, int(start_pos)-1, int(end_pos)):
            yield line.rstrip('\n').split('\t')

    def fetch_hash_snps(self, chromosome, start_pos, end_pos):
        for rec in self.fetch_array_snps(chromosome, start_pos, end_pos):
            if len(rec) != cbv_const.LJB_PARSED_EXPECTED_LENGTH :
                raise Exception("Invalid formatting is found in file '%s'>> Chrom : %s\tpos : %s" % (self.db_file_name, rec[cbv_const.LJB_PARSED_0_IDX_CHROM], rec[cbv_const.LJB_PARSED_0_IDX_POS]))
            snp_info = {cbv_const.KEY_LJB_CHROM : rec[cbv_const.LJB_PARSED_0_IDX_CHROM],
                        cbv_const.KEY_LJB_POS   : rec[cbv_const.LJB_PARSED_0_IDX_POS],
                        cbv_const.KEY_LJB_REF   : rec[cbv_const.LJB_PARSED_0_IDX_REF],
                        cbv_const.KEY_LJB_ALT   : rec[cbv_const.LJB_PARSED_0_IDX_ALT],
                        }
            scores   = {cbv_const.KEY_PHYLOP_SCORE  : rec[cbv_const.LJB_PARSED_0_IDX_PHYLOP_SCORE],
                        cbv_const.KEY_SIFT_SCORE    : rec[cbv_const.LJB_PARSED_0_IDX_SIFT_SCORE],
                        cbv_const.KEY_PP2_SCORE     : rec[cbv_const.LJB_PARSED_0_IDX_PP2_SCORE],
                        cbv_const.KEY_LRT_SCORE     : rec[cbv_const.LJB_PARSED_0_IDX_LRT_SCORE],
                        cbv_const.KEY_MT_SCORE      : rec[cbv_const.LJB_PARSED_0_IDX_MT_SCORE],
                        cbv_const.KEY_GERP_SCORE    : rec[cbv_const.LJB_PARSED_0_IDX_GERP_SCORE],
                        }
            yield {cbv_const.KEY_SNP_INFO_SECTION : snp_info,
                   cbv_const.KEY_SCORES_SECTION   : scores,
                   }

    def get_scores(self, chromosome, pos, ref, alt):
        for rec in self.fetch_array_snps(chromosome, pos, pos):
            if rec[cbv_const.LJB_PARSED_0_IDX_REF] != ref:
                continue
            if rec[cbv_const.LJB_PARSED_0_IDX_ALT] != alt:
                continue
            return {cbv_const.KEY_PHYLOP_SCORE  : rec[cbv_const.LJB_PARSED_0_IDX_PHYLOP_SCORE],
                    cbv_const.KEY_SIFT_SCORE    : rec[cbv_const.LJB_PARSED_0_IDX_SIFT_SCORE],
                    cbv_const.KEY_PP2_SCORE     : rec[cbv_const.LJB_PARSED_0_IDX_PP2_SCORE],
                    cbv_const.KEY_LRT_SCORE     : rec[cbv_const.LJB_PARSED_0_IDX_LRT_SCORE],
                    cbv_const.KEY_MT_SCORE      : rec[cbv_const.LJB_PARSED_0_IDX_MT_SCORE],
                    cbv_const.KEY_GERP_SCORE    : rec[cbv_const.LJB_PARSED_0_IDX_GERP_SCORE],
                    }
        return None

class VcfReader(CombiVEPBase):
    """to read parsed VCF file"""


    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, vcf_file):
        self.vcf_file_name = vcf_file

    def fetch_array_snps(self):
        vcf_file = open(self.vcf_file_name)
        for line in vcf_file:
            if line[0] == '#':
                continue
            yield line.rstrip('\n').split('\t')

    def fetch_hash_snps(self):
        for rec in self.fetch_array_snps():
            snp_info = {cbv_const.KEY_VCF_CHROM : rec[cbv_const.VCF_0_IDX_CHROM],
                        cbv_const.KEY_VCF_POS   : rec[cbv_const.VCF_0_IDX_POS],
                        cbv_const.KEY_VCF_REF   : rec[cbv_const.VCF_0_IDX_REF],
                        cbv_const.KEY_VCF_ALT   : rec[cbv_const.VCF_0_IDX_ALT],
                        }
            yield {cbv_const.KEY_SNP_INFO_SECTION : snp_info}


class CbvReader(CombiVEPBase):
    """

    to read parsed SNPS file
    The format are CHROM, POS, REF, ALT, EFFECT.
    All fields are tab separated.

    """


    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, cbv_file):
        self.cbv_file_name = cbv_file

    def fetch_array_snps(self):
        cbv_file = open(self.cbv_file_name)
        for line in cbv_file:
            if line[0] == '#':
                continue
            yield line.rstrip('\n').split('\t')

    def fetch_hash_snps(self):
        for rec in self.fetch_array_snps():
            snp_info = {cbv_const.KEY_CBV_CHROM : rec[cbv_const.CBV_0_IDX_CHROM],
                        cbv_const.KEY_CBV_POS   : rec[cbv_const.CBV_0_IDX_POS],
                        cbv_const.KEY_CBV_REF   : rec[cbv_const.CBV_0_IDX_REF],
                        cbv_const.KEY_CBV_ALT   : rec[cbv_const.CBV_0_IDX_ALT],
                        }
            prediction = {cbv_const.KEY_CBV_TARGETS : rec[cbv_const.CBV_0_IDX_TARGETS]}
            yield {cbv_const.KEY_SNP_INFO_SECTION   : snp_info,
                   cbv_const.KEY_PREDICTION_SECTION : prediction,
                   }









