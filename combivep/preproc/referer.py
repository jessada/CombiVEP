import combivep.settings as cbv_const
from combivep.config import Configure
from combivep.preproc.reader import UcscReader
from combivep.preproc.reader import LjbReader


class Referer(Configure):
    """To connect to reference database"""

    def __init__(self):
        Configure.__init__(self)

    def load_cfg(self):
        Configure.load_cfg(self)
        ucsc_db = self.cfg_vals[cbv_const.LATEST_UCSC_FILE_NAME]
        ljb_db  = self.cfg_vals[cbv_const.LATEST_LJB_FILE_PREFIX] + '.txt.gz'
        self.__ucsc_reader = UcscReader()
        self.__ucsc_reader.read(ucsc_db)
        self.__ljb_reader = LjbReader()
        self.__ljb_reader.read(ljb_db)

    def validate_snp(self, chrom, pos, ref, alt):
        """

        This function checks if a given snp is valid by referencing
        with UCSC database

        The inputs of this function are in string format
        except "pos", which is integer.

        "chrom" can be either in format "chr1" or "1"
        "pos" is 1-based index

        return True if the snp is presented in UCSC reference database
        and False otherwise.

        """
        reader = self.__get_ucsc_reader()
        for rec in reader.fetch_snps(chrom, int(pos)-1, int(pos)):
            if rec.ref != ref:
                continue
            if ref == alt:
                continue
            ucsc_alts = rec.observed.split('/')
            for ucsc_alt in ucsc_alts:
                if ucsc_alt == alt:
                    return True
        return False

    def get_scores(self, chrom, pos, ref, alt):
        """

        This function returns precomputed prediction scores from LJB database

        The inputs of this function are in string format
        except "pos", which is integer.

        "chrom" can be either in format "chr1" or "1"
        "pos" is 1-based index

        return hash scores if the snp is precomputed and None otherwise

        """
        reader = self.__get_ljb_reader()
        return reader.get_scores(chrom, pos, ref, alt)

    def __get_ucsc_reader(self):
        return self.__ucsc_reader

    def __get_ljb_reader(self):
        return self.__ljb_reader
