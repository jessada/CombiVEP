import subprocess
import os
import re
import sys
import zipfile
import combivep.settings as cbv_const
from combivep.template import CombiVEPBase


def ungz(gz_file):
    args = []
    args.append('gzip')
    args.append('-d')
    args.append(gz_file)
    error = subprocess.call(args)
    if error:
        return None
    file_name, dummy_gz = os.path.splitext(gz_file)
    return file_name


def unzip(zip_file, out_dir):
    """

    unzip file to the output directory
    then return the name of unzipped file together with their path

    """

    unzip_files = zipfile.ZipFile(zip_file)
    out_files = []
    for unzip_file in unzip_files.namelist():
        (dir_name, file_name) = os.path.split(unzip_file)
        unzip_out_dir         = os.path.join(out_dir, dir_name)
        if not os.path.exists(unzip_out_dir):
            os.makedirs(unzip_out_dir)
        out_file    = os.path.join(unzip_out_dir, file_name)
        zipped_file = os.path.join(dir_name, file_name)
        if not out_file.endswith('/'):
            fd = open(out_file,"w")
            fd.write(unzip_files.read(zipped_file))
            fd.close()
            out_files.append(out_file)
    return out_files


class Downloader(CombiVEPBase):
    """to download file"""


    def __init__(self):
        CombiVEPBase.__init__(self)

    def download(self, target_url, output_dir, output_file_name=None):
        """

        to download file from the target url and save it at the specified
        output directory

        """
        current_working_dir = os.getcwd()
        os.chdir(output_dir)
        args = []
        args.append('wget')
        args.append('-q')
        args.append('-N')
        if output_file_name:
            args.append('--output-document=%s' % (output_file_name))
        args.append(target_url)
        error_code = subprocess.call(args)
        os.chdir(current_working_dir)
        return error_code

class Updater(Downloader):


    def __init__(self):
        Downloader.__init__(self)
        self.working_dir = cbv_const.CBV_WORKING_DIR

        #define the property of updater
        #URL of the folder that contain target files
        self.folder_url       = None
        #pattern to find the target files
        self.files_pattern    = None
        #pattern to find the version in each file name
        self.version_pattern  = None
        #directory to store (new) reference DB
        self.local_ref_db_dir = None

        self.tmp_file         = 'tmp_list'

    def check_new_file(self, last_version):
        """

        Check if there is a new version at the target site
        if yes, return the new file together with its version
        other, return 'None'
        Please note that this function only check but not d/l anything

        """
        if not self.__ready():
            return None
        self.create_dir(self.working_dir)
        tmp_list_file  = os.path.join(self.working_dir, self.tmp_file)
        self.download(self.folder_url,
                      self.working_dir,
                      output_file_name=tmp_list_file)
        files_list  = self.__parse(tmp_list_file)
        max_version = max(sorted(files_list.keys()))
        if max_version <= last_version:
            return None, None
        else:
            self.new_version = max_version
            self.new_file    = os.path.join(self.folder_url,
                                            files_list[max_version])
            return self.new_file, self.new_version

    def parse(self, list_file):
        return self.__parse(list_file)

    def __parse(self, list_file):
        """

        extract only the file of interest together with their versions

        """
        out     = {}
        fp      = re.compile(self.files_pattern)
        matches = fp.finditer(open(list_file).read())
        for match in matches:
            vp           = re.compile(self.version_pattern)
            file_name    = match.group('file_name')
            version      = vp.match(file_name).group('version')
            out[version] = file_name
        return out

    def download_new_file(self):
        """download the 'new_file' and store it in 'local_db_dir'"""
        if not os.path.exists(self.local_ref_db_dir):
            os.makedirs(self.local_ref_db_dir)
        self.info('Downloading %s . . .' % (self.new_file))
        error = self.download(self.new_file, self.local_ref_db_dir)
        if error:
            return None
        self.downloaded_file = os.path.join(self.local_ref_db_dir,
                                            os.path.basename(self.new_file))
        return self.downloaded_file

    def __ready(self):
        return self.folder_url and self.files_pattern and self.version_pattern


class UcscUpdater(Updater):
    """ to check if local UCSC DB is up-to-date """


    def __init__(self):
        Updater.__init__(self)

        #specific configuration
        self.folder_url       = cbv_const.UCSC_FOLDER_URL
        self.files_pattern    = cbv_const.UCSC_FILES_PATTERN
        self.version_pattern  = cbv_const.UCSC_VERSION_PATTERN
        self.local_ref_db_dir = cbv_const.USER_UCSC_REF_DB_DIR

    def download_new_file(self):
        if not Updater.download_new_file(self):
            return False
        if self.new_file.endswith('.gz'):
            self.raw_db_file = ungz(self.downloaded_file)
        else:
            self.raw_db_file = self.downloaded_file
        return self.raw_db_file

class LjbUpdater(Updater):
    """ to check if local LJB DB is up-to-date """


    def __init__(self):
        Updater.__init__(self)

        #specific configuration
        self.folder_url       = cbv_const.LJB_FOLDER_URL
        self.files_pattern    = cbv_const.LJB_FILES_PATTERN
        self.version_pattern  = cbv_const.LJB_VERSION_PATTERN
        self.local_ref_db_dir = cbv_const.USER_LJB_REF_DB_DIR

    def download_new_file(self):
        if not Updater.download_new_file(self):
            return False
        if self.new_file.endswith('.zip'):
            self.raw_db_files = unzip(self.downloaded_file,
                                      self.local_ref_db_dir)
        else:
            self.raw_db_files = self.downloaded_file
        return self.raw_db_files






