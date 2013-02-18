import os
import combivep.settings as cbv_const


class Configure(object):
    """ CombiVEP base class """

    def __init__(self):
        self.cfg_file = cbv_const.CBV_CFG_FILE
        self.cfg_vals = {}
        self.cfg_vals[cbv_const.LATEST_UCSC_DB_VER] = '0'
        self.cfg_vals[cbv_const.LATEST_UCSC_FILE_NAME] = ''
        self.cfg_vals[cbv_const.LATEST_LJB_DB_VER] = '0.1'
        self.cfg_vals[cbv_const.LATEST_LJB_FILE_PREFIX] = ''

    def load_cfg(self):
        f = open(self.cfg_file, 'r')
        for line in f:
            rec = line.strip().split('=')
            if rec[0] == cbv_const.LATEST_UCSC_DB_VER:
                self.cfg_vals[cbv_const.LATEST_UCSC_DB_VER] = rec[1]
            elif rec[0] == cbv_const.LATEST_UCSC_FILE_NAME:
                self.cfg_vals[cbv_const.LATEST_UCSC_FILE_NAME] = rec[1]
            elif rec[0] == cbv_const.LATEST_LJB_DB_VER:
                self.cfg_vals[cbv_const.LATEST_LJB_DB_VER] = rec[1]
            elif rec[0] == cbv_const.LATEST_LJB_FILE_PREFIX:
                self.cfg_vals[cbv_const.LATEST_LJB_FILE_PREFIX] = rec[1]
        f.close()
        return self.cfg_vals

    def __save(self):
        f = open(self.cfg_file, 'w')
        cfg_names = []
        cfg_names.append(cbv_const.LATEST_UCSC_DB_VER)
        cfg_names.append(cbv_const.LATEST_UCSC_FILE_NAME)
        cfg_names.append(cbv_const.LATEST_LJB_DB_VER)
        cfg_names.append(cbv_const.LATEST_LJB_FILE_PREFIX)
        for cfg_name in cfg_names:
            cfg_val = self.cfg_vals[cfg_name]
            f.write("{name}={val}\n".format(name=cfg_name,
                                            val=cfg_val))
        f.close()

    def write_ljb_cfg(self, version, file_prefix):
        if os.path.exists(self.cfg_file):
            self.load_cfg()
        self.cfg_vals[cbv_const.LATEST_LJB_DB_VER]      = version
        self.cfg_vals[cbv_const.LATEST_LJB_FILE_PREFIX] = file_prefix
        self.__save()

    def write_ucsc_cfg(self, version, file_name):
        if os.path.exists(self.cfg_file):
            self.load_cfg()
        self.cfg_vals[cbv_const.LATEST_UCSC_DB_VER]    = version
        self.cfg_vals[cbv_const.LATEST_UCSC_FILE_NAME] = file_name
        self.__save()
