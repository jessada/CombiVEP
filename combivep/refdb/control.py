import sys
import os
import pysam
import combivep.settings as cbv_const
from combivep.refdb.updater import UcscUpdater
from combivep.refdb.updater import LjbUpdater
from combivep.config import Configure


class UcscController(UcscUpdater, Configure):
    """UCSC database controller"""


    def __init__(self):
        UcscUpdater.__init__(self)
        Configure.__init__(self)

    def update(self):
        self.load_cfg()
        print >> sys.stderr, 'Checking new UCSC reference database version . . . . '
        new_file, new_version = self.check_new_file(self.cfg_values[cbv_const.LATEST_UCSC_DB_VERSION])
        if not new_version:
            print >> sys.stderr, 'UCSC reference database is already up-to-date (version %s) . . . . . ' % (self.cfg_values[cbv_const.LATEST_UCSC_DB_VERSION])
            return False
        self.download_new_file()
        new_database = self.__tabix_database()
        self.write_ucsc_cfg(new_version, new_database)
        print >> sys.stderr, 'Finish updating UCSC reference database . . . . '
        return True

    def tabix_database(self, file_name):
        """ interface for testing purpose """
        self.raw_db_file = file_name
        return self.__tabix_database()

    def __tabix_database(self):
        return self.__tabix(self.raw_db_file)

    def __tabix(self, file_name):
        """ tabix into gz and tbi file """
        print >> sys.stderr, 'indexing ucsc database . . . . . '
        return pysam.tabix_index(file_name,
                                 force     = True,
                                 seq_col   = cbv_const.UCSC_0_IDX_CHROM,
                                 start_col = cbv_const.UCSC_0_IDX_START_POS,
                                 end_col   = cbv_const.UCSC_0_IDX_END_POS,
                                 zerobased = True)


class LjbController(LjbUpdater, Configure):
    """LJB database controller"""


    def __init__(self):
        LjbUpdater.__init__(self)
        Configure.__init__(self)
        self.chromosome_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    def update(self):
        self.load_cfg()
        print >> sys.stderr, 'Checking new LJB reference database version . . . . '
        new_file, new_version = self.check_new_file(self.cfg_values[cbv_const.LATEST_LJB_DB_VERSION])
        if not new_version:
            print >> sys.stderr, 'LJB reference database is already up-to-date (version %s) . . . . . ' % (self.cfg_values[cbv_const.LATEST_LJB_DB_VERSION])
            return False
        self.download_new_file()
        file_prefix, dummy_ext = os.path.splitext(self.downloaded_file)
        self.delete_file(file_prefix + '.txt')
        #clean and concat then tabix
        print >> sys.stderr, 'cleaning database . . . . . '
        for chromosome_file in self.__get_chromosome_files(file_prefix):
            self.__clean_raw_database(chromosome_file, chromosome_file + '.clean')
            self.__concat_file(chromosome_file + '.clean', file_prefix + '.txt')
        print >> sys.stderr, 'indexing ljb database . . . . . '
        self.__tabix_database(file_prefix + '.txt')
        #save file information to the configuration file
        self.write_ljb_cfg(new_version, file_prefix)
        #remove downloaded and temporary files
        self.delete_file(self.downloaded_file)
        for chromosome_file in self.__get_chromosome_files(file_prefix):
            self.delete_file(chromosome_file)
            self.delete_file(chromosome_file + '.clean')
        print >> sys.stderr, 'Finish updating LJB reference database . . . . '
        return True

    def __get_chromosome_files(self, file_prefix):
        for chromosome in self.chromosome_list:
            yield file_prefix + '.chr' + chromosome

    def concat_chromosome_files(self, file_prefix, file_suffix, out_file):
        return self.__concat_chromosome_files(file_prefix, file_suffix, out_file)

    def __concat_chromosome_files(self, file_prefix, file_suffix, out_file):
        self.delete_file(out_file)
        for chromosome_file in self.__get_chromosome_files(file_prefix):
            self.__concat_file(chromosome_file + file_suffix, out_file)

    def __concat_file(self, source, target):
        cmd = []
        cmd.append(' cat ')
        cmd.append(source)
        cmd.append(' >> ')
        cmd.append(target)
        return self.exec_sh(''.join(cmd))

    def clean_raw_database(self, input_file, output_file):
        """ interface for testing purpose """
        return self.__clean_raw_database(input_file, output_file)

    def __clean_raw_database(self, input_file, output_file):
        cmd = []
        #remove records that any of the scores are 'NA'
        cmd.append('awk -F\'\\t\' \'($')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_PHYLOP_SCORE))
        cmd.append(' !~ /[A-Z]/) && ($')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_SIFT_SCORE))
        cmd.append(' !~ /[A-Z]/) && ($')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_PP2_SCORE))
        cmd.append(' !~ /[A-Z]/) && ($')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_LRT_SCORE))
        cmd.append(' !~ /[A-Z]/) && ($')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_MT_SCORE))
        cmd.append(' !~ /[A-Z]/) && ($')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_GERP_SCORE))
        cmd.append(' !~ /[A-Z]/)\' ')
        cmd.append(input_file)
        #reformat the columns
        cmd.append(' | awk -F\'\\t\' \'{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_CHROM))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_POS))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_REF))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_ALT))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_PHYLOP_SCORE))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_SIFT_SCORE))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_PP2_SCORE))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_LRT_SCORE))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_MT_SCORE))
        cmd.append(', $')
        cmd.append(str(cbv_const.LJB_RAW_1_IDX_GERP_SCORE))
        cmd.append('}\' ')
        #redirect output to a file
        cmd.append(' | sort -k2 -n ')
        #redirect output to a file
        cmd.append(' > ')
        cmd.append(output_file)
        return self.exec_sh(''.join(cmd))

    def tabix_database(self, file_name):
        """ interface for testing purpose """
        self.__tabix_database(file_name)

    def __tabix_database(self, file_name):
        return self.__tabix(file_name)

    def __tabix(self, file_name):
        """ tabix into gz and tbi file """
        return pysam.tabix_index(file_name,
                                 force     = True,
                                 seq_col   = cbv_const.LJB_PARSED_0_IDX_CHROM,
                                 start_col = cbv_const.LJB_PARSED_0_IDX_POS,
                                 end_col   = cbv_const.LJB_PARSED_0_IDX_POS,
                                 zerobased = False)







