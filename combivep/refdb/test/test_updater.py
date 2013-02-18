import unittest
import os
import combivep.refdb.test.template as test_template
import combivep.settings as cbv_const
import combivep.refdb.updater as combivep_updater


class TestDownloader(test_template.SafeRefDBTester):

    def __init__(self, test_name):
        test_template.SafeRefDBTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'downloader'

    def __create_downloader_instance(self):
        downloader = combivep_updater.Downloader()
        return downloader

    def test_download(self):
        self.init_test('test_download')
        downloader = self.__create_downloader_instance()
        target_url = 'http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP_light_v1.3.readme.txt'
        target_file = os.path.join(self.working_dir,
                                   os.path.basename(target_url))
        downloader.download(target_url, self.working_dir)
        self.assertGreater(os.stat(target_file).st_size,
                           3397,
                           msg='Download does not functional properly')

    def tearDown(self):
        self.remove_working_dir()


class TestUpdater(test_template.SafeRefDBTester):

    def __init__(self, test_name):
        test_template.SafeRefDBTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'updater'

    def __create_updater_instance(self):
        updater = combivep_updater.Updater()
        updater.working_dir = self.working_dir
        return updater

    def test_ready1(self):
        self.init_test('test_ready1')
        updater = self.__create_updater_instance()
        self.assertEqual(updater.check_new_file('135'),
                         None,
                         msg='error in checking Updater ready state')

    def test_ready2(self):
        self.init_test('test_ready2')
        updater = self.__create_updater_instance()
        updater.folder_url      = 'http://dummy_url'
        updater.version_pattern = 'dummy_version_pattern'
        self.assertFalse(updater.check_new_file('135'),
                         msg='error in checking Updater ready state')

    def test_ready3(self):
        self.init_test('test_ready3')
        updater = self.__create_updater_instance()
        updater.files_pattern   = r"""href="(?P<file_name>snp\d{3}.sql)">.*>.*(?P<date>\d{2}-[a-zA-Z]{3}-\d{4})"""
        updater.version_pattern = 'dummy_version_pattern'
        self.assertFalse(updater.check_new_file('135'),
                         msg='error in checking Updater ready state')

    def tearDown(self):
        self.remove_working_dir()


#@unittest.skip("temporary disable due to high bandwidth usage")
class TestUcscUpdater(test_template.SafeRefDBTester):

    def __init__(self, test_name):
        test_template.SafeRefDBTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'ucsc_updater'

    def __create_ucsc_updater_instance(self):
        ucsc_updater = combivep_updater.UcscUpdater()
        ucsc_updater.working_dir      = self.working_dir
        ucsc_updater.files_pattern    = r"""href="(?P<file_name>snp\d{3}.sql)">.*>.*(?P<date>\d{2}-[a-zA-Z]{3}-\d{4})"""
        ucsc_updater.tmp_file         = cbv_const.UCSC_LIST_FILE_NAME
        ucsc_updater.local_ref_db_dir = self.working_dir
        return ucsc_updater

    @unittest.skip("temporary disable due to high bandwidth usage")
    def test_update1(self):
#        self.individual_debug = True
        self.init_test('test_update1')
        ucsc_updater = self.__create_ucsc_updater_instance()
        new_file, new_version = ucsc_updater.check_new_file('135')
        self.assertTrue(new_file.endswith('.sql'),
                        msg='incorrectly identify updating result')

    @unittest.skip("temporary disable due to high bandwidth usage")
    def test_update2(self):
        self.init_test('test_update2')
        ucsc_updater = self.__create_ucsc_updater_instance()
        new_file, new_version = ucsc_updater.check_new_file('136')
        self.assertTrue(new_file.endswith('.sql'),
                        msg='incorrectly identify updating result')

    @unittest.skip("temporary disable due to high bandwidth usage")
    def test_not_update1(self):
        self.init_test('test_not_update1')
        ucsc_updater = self.__create_ucsc_updater_instance()
        new_file, new_version = ucsc_updater.check_new_file('137')
        self.assertEqual(new_version,
                         None,
                         msg='incorrectly identify updating result')

    @unittest.skip("temporary disable due to high bandwidth usage")
    def test_not_update2(self):
        self.init_test('test_not_update2')
        ucsc_updater = self.__create_ucsc_updater_instance()
        new_file, new_version = ucsc_updater.check_new_file('138')
        self.assertequal(new_version,
                         none,
                         msg='incorrectly identify updating result')

    @unittest.skip("temporary disable due to high bandwidth usage")
    def test_full_update1(self):
        #init
        self.init_test('test_full_update1')
        ucsc_updater = self.__create_ucsc_updater_instance()
        #run test
        new_file, new_version = ucsc_updater.check_new_file('135')
        self.asserttrue(new_file.endswith('.sql'),
                        msg='UCSC updater got an invalid new file')
        self.assertEqual(new_version,
                         '137',
                         msg='UCSC updater got an invalid new version number')
        downloaded_file = ucsc_updater.download_new_file()
        self.assertTrue(os.path.exists(downloaded_file),
                        msg="UCSC's downloaded file (%s) doesn't exist" % (downloaded_file))

    def test_parse(self):
        self.init_test('test_parse')
        ucsc_updater = self.__create_ucsc_updater_instance()
        test_file = os.path.join(self.data_dir, 'dummy_ucsc_list_file')
        out       = ucsc_updater.parse(test_file)
        self.assertEqual(out['135'],
                         'snp135.sql',
                         msg='UCSC parser cannot find UCSC files')

    def tearDown(self):
        self.remove_working_dir()


class TestLJBUpdater(test_template.SafeRefDBTester):

    def __init__(self, test_name):
        test_template.SafeRefDBTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'ljb_updater'

    def __create_ljb_updater_instance(self):
        ljb_updater = combivep_updater.LjbUpdater()
        ljb_updater.working_dir      = self.working_dir
        ljb_updater.files_pattern    = r"""href="(?P<file_name>dbNSFPv[\d.]*.readme.txt)">"""
        ljb_updater.tmp_file         = cbv_const.LJB_LIST_FILE_NAME
        ljb_updater.local_ref_db_dir = self.working_dir
        return ljb_updater

    def test_update1(self):
        self.init_test('test_update1')
        ljb_updater = self.__create_ljb_updater_instance()
        new_file, new_version = ljb_updater.check_new_file('1.2')
        self.assertNotEqual(new_version,
                            None,
                            msg='incorrectly identify updating result')

    def test_not_update1(self):
        self.init_test('test_not_update1')
        ljb_updater = self.__create_ljb_updater_instance()
        new_file, new_version = ljb_updater.check_new_file('1.3')
        self.assertEqual(new_version,
                         None,
                         msg='incorrectly identify updating result')

    def test_not_update2(self):
        self.init_test('test_not_update2')
        ljb_updater = self.__create_ljb_updater_instance()
        new_file, new_version = ljb_updater.check_new_file('1.4')
        self.assertEqual(new_version,
                         None,
                         msg='incorrectly identify updating result')

    def test_full_update1(self):
        #init
        self.init_test('test_full_update1')
        ljb_updater = self.__create_ljb_updater_instance()
        #run test
        new_file, new_version = ljb_updater.check_new_file('1.2')
        self.assertTrue(new_file.endswith('.readme.txt'),
                        msg='LJB updater got an invalid new file')
        self.assertEqual(new_version,
                         '1.3',
                         msg='LJB updater got an invalid new version number')
        downloaded_file = ljb_updater.download_new_file()
        self.assertTrue(os.path.exists(downloaded_file),
                        msg="LJB's downloaded file (%s) doesn't exist" % (downloaded_file))

    def test_parse(self):
        self.init_test('test_parse')
        ljb_updater = self.__create_ljb_updater_instance()
        test_file = os.path.join(self.data_dir, 'dummy_ljb_list_file')
        out       = ljb_updater.parse(test_file)
        self.assertEqual(out['1.3'],
                         'dbNSFPv1.3.readme.txt',
                         msg='LJB parser cannot find LJB files')

    def tearDown(self):
        self.remove_working_dir()


class TestMisc(test_template.SafeRefDBTester):
    """ test (a few) miscellaneous function(s) """

    def __init__(self, test_name):
        test_template.SafeRefDBTester.__init__(self, test_name)

    def setUp(self):
        self.test_class = 'misc'

    def test_unzip(self):
        self.init_test('test_unzip')
        test_file = os.path.join(self.data_dir, 'dummy_ljb.zip')
        out = combivep_updater.unzip(test_file, self.working_dir)
        self.assertEqual(out[0],
                         os.path.join(self.working_dir,
                                      'try.in'),
                         msg='some files are missing')
        self.assertEqual(out[1],
                         os.path.join(self.working_dir,
                                      'search_dbNSFP_light_v1.3.readme.pdf'),
                         msg='some files are missing')
        self.assertEqual(out[2],
                         os.path.join(self.working_dir,
                                      'search_dbNSFP_light_v1.3.readme.doc'),
                         msg='some files are missing')
        self.assertEqual(out[3],
                         os.path.join(self.working_dir,
                                      'search_dbNSFP_light13.class'),
                         msg='some files are missing')
        self.assertEqual(out[4],
                         os.path.join(self.working_dir,
                                      'dbNSFP_light_v1.3.readme.txt'),
                         msg='some files are missing')
        self.assertEqual(out[5],
                         os.path.join(self.working_dir,
                                      'abc/tryhg19.in'),
                         msg='some files are missing')
        for out_file in out:
            self.assertTrue(os.path.exists(out_file),
                            msg='"%s" file is missing' % (out_file))

    def test_gz(self):
        self.init_test('test_gz')
        test_file    = os.path.join(self.data_dir, 'test_gz.txt.gz')
        working_file = os.path.join(self.working_dir, 'test_gz.txt.gz')
        self.copy_file(test_file, working_file)
        out = combivep_updater.ungz(working_file)
        self.assertTrue(out,
                        'no output file from "ungz" function')
        self.assertTrue(os.path.exists(out),
                        'ungz function does not work properly')

    def tearDown(self):
        self.remove_working_dir()
