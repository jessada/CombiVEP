import pysam
import combivep.settings as cbv_const
from combivep.template import CombiVEPBase
from collections import namedtuple


class UcscRecord(CombiVEPBase):
    """ to automatically parse ucsc data """


    def __init__(self, array_snp):
        self.__array_snp = array_snp

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.__array_snp)

    @property
    def chrom(self):
        return self.__array_snp[cbv_const.UCSC_0_IDX_CHROM]

    @property
    def start_pos(self):
        return self.__array_snp[cbv_const.UCSC_0_IDX_START_POS]

    @property
    def end_pos(self):
        return self.__array_snp[cbv_const.UCSC_0_IDX_END_POS]

    @property
    def strand(self):
        return self.__array_snp[cbv_const.UCSC_0_IDX_STRAND]

    @property
    def ref(self):
        return self.__array_snp[cbv_const.UCSC_0_IDX_REF]

    @property
    def observed(self):
        return self.__array_snp[cbv_const.UCSC_0_IDX_OBSERVED]

class UcscReader(CombiVEPBase):
    """to read UCSC parsed file in tabix format"""

    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, tabix_file):
        self.db_file_name = tabix_file

    def fetch_snps(self, chromosome, start_pos, end_pos):
        tabix_file = pysam.Tabixfile(self.db_file_name)
        if chromosome.startswith('chr'):
            chrom = chromosome
        else:
            chrom = 'chr' + chromosome
        for line in tabix_file.fetch(chrom, start_pos, end_pos):
            raw_rec = line.rstrip('\n').split('\t')
            if len(raw_rec) != cbv_const.UCSC_EXPECTED_LENGTH:
                self.throw(self.format_error.format(data=line.rstrip('\n')))
            yield UcscRecord(raw_rec)

    @property
    def format_error(self):
        return "Invalid formatting is found at:\n{data}"

class LjbRecord(CombiVEPBase):
    """ to automatically parse ljb data """


    def __init__(self, array_snp):
        self.__array_snp = array_snp

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.__array_snp)

    @property
    def chrom(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_CHROM]

    @property
    def pos(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_POS]

    @property
    def ref(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_REF]

    @property
    def alt(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_ALT]

    @property
    def phylop_score(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_PHYLOP_SCORE]

    @property
    def sift_score(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_SIFT_SCORE]

    @property
    def pp2_score(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_PP2_SCORE]

    @property
    def lrt_score(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_LRT_SCORE]

    @property
    def mt_score(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_MT_SCORE]

    @property
    def gerp_score(self):
        return self.__array_snp[cbv_const.LJB_PARSED_0_IDX_GERP_SCORE]

class ScoresRecord(CombiVEPBase):
    """ to store precalculated scores """


    def __init__(self, rec):
        self.__rec = rec

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str({'phylop_score': self.phylop_score,
                    'sift_score': self.sift_score,
                    'pp2_score': self.pp2_score,
                    'lrt_score': self.lrt_score,
                    'mt_score': self.mt_score,
                    'gerp_score': self.gerp_score,
                    })

    def __len__(self):
        return 6

    @property
    def phylop_score(self):
        return self.__rec.phylop_score

    @property
    def sift_score(self):
        return self.__rec.sift_score

    @property
    def pp2_score(self):
        return self.__rec.pp2_score

    @property
    def lrt_score(self):
        return self.__rec.lrt_score

    @property
    def mt_score(self):
        return self.__rec.mt_score

    @property
    def gerp_score(self):
        return self.__rec.gerp_score

class LjbReader(CombiVEPBase):
    """to read parsed LJB file"""


    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, ljb_file):
        self.db_file_name = ljb_file
        self.__tabix_file = pysam.Tabixfile(ljb_file)

    def fetch_snps(self, chromosome, start_pos, end_pos):
        for line in self.__tabix_file.fetch(chromosome,
                                            int(start_pos)-1,
                                            int(end_pos)):
            raw_rec = line.rstrip('\n').split('\t')
            if len(raw_rec) != cbv_const.LJB_PARSED_EXPECTED_LENGTH:
                self.throw(self.format_error.format(data=line.rstrip('\n')))
            yield LjbRecord(raw_rec)

    def get_scores(self, chromosome, pos, ref, alt):
        for rec in self.fetch_snps(chromosome, pos, pos):
            if rec.ref != ref:
                continue
            if rec.alt != alt:
                continue
            return ScoresRecord(rec)
        return None

    @property
    def format_error(self):
        return "Invalid formatting is found at:\n{data}"

class VcfRecord(CombiVEPBase):
    """ to automatically parse vcf data"""


    def __init__(self, array_snp):
        self.__array_snp = array_snp

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.__array_snp)

    @property
    def chrom(self):
        return self.__array_snp[cbv_const.VCF_0_IDX_CHROM]

    @property
    def pos(self):
        return self.__array_snp[cbv_const.VCF_0_IDX_POS]

    @property
    def ref(self):
        return self.__array_snp[cbv_const.VCF_0_IDX_REF]

    @property
    def alt(self):
        return self.__array_snp[cbv_const.VCF_0_IDX_ALT]

class VcfReader(CombiVEPBase):
    """to read parsed VCF file"""


    def __init__(self):
        CombiVEPBase.__init__(self)

    def read(self, vcf_file):
        self.vcf_file_name = vcf_file

    def fetch_snps(self):
        vcf_file = open(self.vcf_file_name)
        for line in vcf_file:
            if line[0] == '#':
                continue
            yield VcfRecord(line.rstrip('\n').split('\t'))

class CbvRecord(CombiVEPBase):
    """ to automatically parse vcf data"""


    def __init__(self, array_snp):
        self.__array_snp = array_snp

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.__array_snp)

    @property
    def chrom(self):
        return self.__array_snp[cbv_const.CBV_0_IDX_CHROM]

    @property
    def pos(self):
        return self.__array_snp[cbv_const.CBV_0_IDX_POS]

    @property
    def ref(self):
        return self.__array_snp[cbv_const.CBV_0_IDX_REF]

    @property
    def alt(self):
        return self.__array_snp[cbv_const.CBV_0_IDX_ALT]

    @property
    def targets(self):
        return self.__array_snp[cbv_const.CBV_0_IDX_TARGETS]

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

    def fetch_snps(self):
        cbv_file = open(self.cbv_file_name)
        for line in cbv_file:
            if line[0] == '#':
                continue
            yield CbvRecord(line.rstrip('\n').split('\t'))

