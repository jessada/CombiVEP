import os
import combivep.settings as cbv_const


class Configure(object):
    """ CombiVEP base class """


    def __init__(self):
        self.cfg_file = cbv_const.CBV_CFG_FILE
        self.cfg_values  = {}
        self.cfg_values[cbv_const.LATEST_UCSC_DB_VERSION] = '0'
        self.cfg_values[cbv_const.LATEST_UCSC_FILE_NAME]  = ''
        self.cfg_values[cbv_const.LATEST_LJB_DB_VERSION]  = '0.1'
        self.cfg_values[cbv_const.LATEST_LJB_FILE_PREFIX] = ''

    def load_cfg(self):
        f = open(self.cfg_file, 'r')
        for line in f:
            rec = line.strip().split('=')
            if rec[0] == cbv_const.LATEST_UCSC_DB_VERSION:
                self.cfg_values[cbv_const.LATEST_UCSC_DB_VERSION] = rec[1]
            elif rec[0] == cbv_const.LATEST_UCSC_FILE_NAME:
                self.cfg_values[cbv_const.LATEST_UCSC_FILE_NAME]  = rec[1]
            elif rec[0] == cbv_const.LATEST_LJB_DB_VERSION:
                self.cfg_values[cbv_const.LATEST_LJB_DB_VERSION]  = rec[1]
            elif rec[0] == cbv_const.LATEST_LJB_FILE_PREFIX:
                self.cfg_values[cbv_const.LATEST_LJB_FILE_PREFIX] = rec[1]
        f.close()
        return self.cfg_values

    def __save(self):
        f = open(self.cfg_file, 'w')
        f.write("%s=%s\n" % (cbv_const.LATEST_UCSC_DB_VERSION, self.cfg_values[cbv_const.LATEST_UCSC_DB_VERSION]))
        f.write("%s=%s\n" % (cbv_const.LATEST_UCSC_FILE_NAME,  self.cfg_values[cbv_const.LATEST_UCSC_FILE_NAME]))
        f.write("%s=%s\n" % (cbv_const.LATEST_LJB_DB_VERSION,  self.cfg_values[cbv_const.LATEST_LJB_DB_VERSION]))
        f.write("%s=%s\n" % (cbv_const.LATEST_LJB_FILE_PREFIX, self.cfg_values[cbv_const.LATEST_LJB_FILE_PREFIX]))
        f.close()

    def write_ljb_cfg(self, version, file_prefix):
        if os.path.exists(self.cfg_file):
            self.load_cfg()
        self.cfg_values[cbv_const.LATEST_LJB_DB_VERSION]  = version
        self.cfg_values[cbv_const.LATEST_LJB_FILE_PREFIX] = file_prefix
        self.__save()

    def write_ucsc_cfg(self, version, file_name):
        self.load_cfg()
        self.cfg_values[cbv_const.LATEST_UCSC_DB_VERSION] = version
        self.cfg_values[cbv_const.LATEST_UCSC_FILE_NAME]  = file_name
        self.__save()






