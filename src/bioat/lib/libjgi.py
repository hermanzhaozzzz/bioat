import gzip
import json
import os
import re
import sys
import tarfile
import textwrap
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from hashlib import md5

import requests
from requests import HTTPError
from requests.exceptions import ChunkedEncodingError
from requests.utils import cookiejar_from_dict
from tqdm import tqdm
from urllib3.exceptions import InvalidChunkLength, ProtocolError

from bioat.lib.libpath import HOME
from bioat.logger import get_logger

__module_name__ = 'bioat.lib.libjgi'


class JGIDoc:
    DEFAULT_CATEGORIES = ["ESTs",
                          "EST Clusters",
                          "Assembled scaffolds (unmasked)",
                          "Assembled scaffolds (masked)",
                          "Transcripts",
                          "Genes",
                          "CDS",
                          "Proteins",
                          "Additional Files"]
    usage_example_blurb = """\
        This tool will retrieve files from JGI
        It will return a list of possible files for downloading.
    
        To get <jgi_address>, go to: http://genome.jgi.doe.gov/ and search for your
        species of interest. Click through until you are at the "Info" page. For
        \x1B[3mNematostella vectensis\x1B[23m, the appropriate page is
        "http://genome.jgi.doe.gov/Nemve1/Nemve1.info.html".
    
        To query using only the name simply requires the specific JGI organism
        abbreviation, as referenced in the full url.
    
        For the above example, the proper input syntax for this script would be:
    
        $ bioat meta JGI_query -q http://genome.jgi.doe.gov/Nemve1/Nemve1.info.html
    
                                 -or-
    
        $ bioat meta JGI_query -q Nemve1
    
        If you already have the XML file for the query in the directory, you may use
        the --xml flag to avoid redownloading it (particularly useful if querying
        large, top-level groups with many sub-species, such as "fungi"):
    
        $ bioat meta JGI_query -x <your_xml_index>
        
        If the XML filename is omitted when using the -x/--xml flag, it is assumed that
        the XML file is named "jgi-xml-query.result_<organism-name>.xml". In such cases, the
        organism name is required.
        """
    select_blurb = """\
        # SYNTAX ///////////////////////////////////////////////////////////////////////
        Using the following format syntax to download selected file:
            <category number>:<i>[,<i>, <i>];<category number>:<i>-<i>;...
        
        Indices (<i>) may be a mixture of comma-separated values and hyphen-separated 
        ranges. 
        
        Example:
            '3:4,5; 7:1-10,13' will select elements 4 and 5 from category 3, and 1-10 
            plus 13 from category 7.
        # /SYNTAX ///////////////////////////////////////////////////////////////////////
        """


class JGIConfig:
    URL_JGI_MAIN = "https://genome.jgi.doe.gov"  # url home for jgi-img database
    URL_JGI_LOGIN = 'https://signon.jgi.doe.gov/signon/create'  # url login
    URL_JGI_FETCH_XML = "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory"  # url fetch xml
    FILENAME_TEMPLATE_COOKIE = "jgi-xml-query.cookie_{}.cookie"
    FILENAME_TEMPLATE_XML = "jgi-xml-query.result_{}.xml"
    FILENAME_TEMPLATE_LOG_FAIL = "jgi-xml-query.failed_{}.log"
    FILENAME_CONFIG_PATH = os.path.join(HOME, '.bioat', 'JGI', "account.conf")

    def __init__(self, overwrite_conf: bool = False, logger=None):
        if not logger:
            self.logger = get_logger(level='DEBUG', module_name=__module_name__,
                                     func_name=sys._getframe().f_code.co_name)
        else:
            self.logger = logger

        self.info: dict = {"user": None, "password": None, "categories": None}

        self.logger.debug('start to set / load user info...')
        if overwrite_conf:
            # 只要指定重写参数就重写，无需判断其他
            self.input_user_info()  # 从终端手动输入JGI用户信息
            self.save_config()  # 保存用户信息
            # 然后直接退出
            self.logger.info("Configuration complete.\n"
                             "Script may now be used to query JGI.\n"
                             "You need re-run your command without parameter `--overwrite_conf`\n")
            # self.load_config()  # 从文件加载用户信息
        else:
            if os.path.isfile(self.FILENAME_CONFIG_PATH) and os.path.getsize(self.FILENAME_CONFIG_PATH) > 0:
                # 直接加载配置信息
                self.load_config()
            else:
                # 配置信息保存文件异常
                self.logger.warning(f'errors occur when loading file: {self.FILENAME_CONFIG_PATH}!\n'
                                    f'try to get user info by manual inputting...')
                self.input_user_info()  # 从终端手动输入JGI用户信息
                self.save_config()  # 保存用户信息
                self.load_config()  # 从文件加载用户信息
        self.logger.debug('set / load user info success!')

    def load_config(self):
        """
        Reads "user", "password" and "categories" entries
        from config file.

        """
        self.logger.debug('loading JGI account info...')

        with open(self.FILENAME_CONFIG_PATH, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith("user"):
                    self.info["user"] = line.split("=")[1]
                if line.startswith("password"):
                    self.info["password"] = line.split("=")[1]
                if line.startswith("categories"):
                    cats = line.strip().split("=")[1]
                    self.info["categories"] = [e.strip() for e in cats.split(",")]

        if not all([self.info["user"], self.info["password"]]):
            self.logger.critical(
                f"Config file present ({self.FILENAME_CONFIG_PATH}), but user and/or password not found.")
            sys.exit(1)

        self.logger.debug('loading JGI account info done')

    def input_user_info(self):
        """
        Dialog with user to gather user information for
        use with the curl query. Returns a dict.

        """
        blurb = """
        === USER SETUP ===

        JGI access configuration:

        Before continuing, you will need to provide your JGI login credentials.
        These are required by JGI's curl api, and will be stored in a config
        file for future use (unless you choose to delete them).

        If you need to sign up for a JGI account, use the registration link at
        https://contacts.jgi.doe.gov/registration/new

        === CREDENTIALS ===
        """
        # Print left-justified triple-quoted text blocks
        print(textwrap.dedent(blurb))
        user_query = "JGI account username/email (or 'q' to quit): "

        # write user and passwd from interactive terminal
        self.info["user"] = input(user_query)

        if self.info["user"] == "q":
            sys.exit("Exiting now.")

        self.info["password"] = input("JGI account password (or 'q' to quit): ")

        if self.info["password"] == "q":
            sys.exit("Exiting now.")

        if self.info["user"].strip() == '' or self.info["password"] == '':
            self.logger.warning('user and password can not be None!')
            self.input_user_info()

        input_blurb = (
            f"Proceed with USER={self.info['user']}, PASSWORD={self.info['password']} to configure "
            "script?\n([y]es, [n]o, [r]estart): "
        )

        self.info["categories"]: list[str] = JGIDoc.DEFAULT_CATEGORIES

        while True:  # catch invalid responses
            choice = input(input_blurb)

            if choice.lower() == "y":
                break
            elif choice.lower() == "n":
                sys.exit("Exiting now.")
            elif choice.lower() == "r":
                self.input_user_info()

    def save_config(self):
        """
        Creates a config file <config_path> using
        credentials from dict <config_info>.

        """
        self.logger.debug('saving JGI account info...')
        u = self.info["user"]
        p = self.info["password"]
        c = self.info["categories"]
        c = ",".join(c)
        header = ("# bioat meta JGI_query: JGI account info and essentials {}\n".format("#" * 24))
        info = f"user={u}\npassword={p}\ncategories={c}"

        if not os.path.exists(os.path.dirname(self.FILENAME_CONFIG_PATH)):
            os.makedirs(os.path.dirname(self.FILENAME_CONFIG_PATH))

        with open(self.FILENAME_CONFIG_PATH, 'wt') as f:
            f.write(header)
            f.write(info)

        self.logger.debug('saving JGI account info done')


class JGIOperator:
    def __init__(
            self,
            # pick one from three
            query_info: str | None = None,
            xml: str | None = None,
            log_fails: str | None = None,
            # runtime params
            nretry: int = 5,
            timeout: int = -1,
            regex: str | None = None,
            all_get: bool = False,
            overwrite_conf: bool = False,
            filter_files: bool = False,
            # doc helper
            syntax_help: bool = False,
            usage: bool = False,
            # log
            log_level: str = 'INFO'
    ):
        # from cmd parameters
        self.query_info = query_info
        self.xml = xml
        self.log_fails = log_fails
        self.nretry = nretry
        self.timeout = None if timeout == -1 else timeout
        self.regex = regex
        self.all_get = all_get
        self.overwrite_conf = overwrite_conf
        self.filter_files = filter_files  # TODO 搞清楚用法
        self.syntax_help = syntax_help
        self.usage = usage
        self.log_level = log_level
        # self
        self._dict_to_get = dict()
        self._urls_to_get = set()
        self._url_to_validate = defaultdict(dict)
        self._selections = dict()
        self._downloaded_files = list()
        self._failed_urls = list()
        self._desired_categories = dict()
        # From other obj #
        # load configs; auto check if you need overwrite user info or not
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name='class: JGIConfig')
        self.config = JGIConfig(overwrite_conf=overwrite_conf, logger=logger)
        self.interactive = False
        # load docs;
        self.docs = JGIDoc
        # / From other obj

    # step 01 print and exit
    def run_doc(self):
        """Checker for doc mode.

        print syntax_help and/or usage if these parameters are defined. And then, exit progress.
        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        run_docs = any([self.syntax_help, self.usage])
        if run_docs:
            logger.debug('doc mode is selected')
            # Check if user wants query help
            if self.syntax_help:
                print(f'\n[syntax_help]:\n{self.docs.select_blurb}')
            if self.usage:
                print(f'\n[usage]:\n{self.docs.usage_example_blurb}')
            sys.exit('Done. exit.')
            # finally exit

    # step 02 update self.query_info and self.log_fails and followed by self.login
    def parse_input(self):
        """Checker for input.

        A checker for this rule: ONLY ONE of the parameters ('query_info', 'xml' and 'log_fails') can be specified;
        Parse query_info/xml/log_fails to update self.query_info and self.log_fails.
        After this method, you should call self.login method and then self.query method
        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        logger.debug("update JGIOperator obj.query_info and obj.log_fails using parameters:"
                     " 'query_info', 'xml' and 'log_fails'")

        if sum(map(bool, [self.query_info, self.xml, self.log_fails])) != 1:
            # 这三个参数只可以指定一个！
            logger.error("ONLY ONE of the parameters ('query_info', 'xml' and 'log_fails') can be specified!")
            sys.exit('Done. exit.')
        else:
            logger.debug("pass param checker: 'query_info', 'xml' and 'log_fails'")

        if self.query_info:
            try:
                # if query_info is a URL
                logger.debug('attempt to parse query_info as a url')
                query_regex = re.compile(r'\.jgi.+\.(?:gov|org).*\/(.+)\/(?!\/)')
                self.query_info = query_regex.search(self.query_info).group(1)
            except AttributeError:
                # if query_info is an organism name
                logger.debug('not a url, try to define query_info as an organism name abbreviation')
                # query_info = query_info
        elif self.log_fails:
            # load failed info from log file if provided
            # filename see: class JGIConfig.FILENAME_TEMPLATE_LOG_FAIL
            self.query_info = os.path.basename(self.log_fails).split('.')[1].split('_')[0]

            with open(self.log_fails, 'rt') as f:
                self.log_fails = f.read().splitlines()
        elif self.xml:
            # parse query_info from xml file content
            name_pattern = r"name=\"(.+)\""
            _org_line = None
            with open(self.xml, 'rt') as f:
                for line in f:
                    if "organismDownloads" in line:
                        # standardized name indicator
                        # <organismDownloads name="Nemve1">
                        _org_line = line.strip()
                        break  # don't keep looking, already found
            try:
                self.query_info = re.search(name_pattern, _org_line).group(1)
            except TypeError:  # org_line still None
                logger.critical('the xml file seems wrong')
                sys.exit('Exit with errors.')
        else:
            logger.error("one of the parameters 'query_info' and 'xml' should be specified")
            sys.exit('Done. exit.')

    # step 03 login
    def login(self):
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        # prepare info
        url = self.config.URL_JGI_LOGIN
        cookie_file = self.config.FILENAME_TEMPLATE_COOKIE.format(self.query_info)
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/72.0.3626.109 Safari/537.36'
        }
        data = {
            'login': self.config.info['user'],
            'password': self.config.info['password'],
            'commit': 'Sign In'
        }
        # set session
        s = requests.session()
        # get response
        response = s.post(url, headers=headers, data=data, timeout=self.timeout)

        # check response status code
        logger.debug(f'check response status code: {response.status_code}')
        if response.status_code == 200:
            logger.debug('login successes, get cookie...')
            # save cookie to file
            logger.debug(f'save cookie to {cookie_file}')
            cookies_dict = requests.utils.dict_from_cookiejar(s.cookies)
            cookies_str = json.dumps(cookies_dict)

            with open(cookie_file, 'wt') as f:
                f.write(cookies_str)
            logger.debug(f'successfully login, write cookie @ {cookie_file}')
        else:
            logger.critical("Couldn't connect with server. Please check Internet connection and nretry.")
            self._clean_exit(logger=logger)

        # step 04 query xml info

    # step 04 query
    def query(self):
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        # prepare info
        url = self.config.URL_JGI_FETCH_XML
        cookie_file = self.config.FILENAME_TEMPLATE_COOKIE.format(self.query_info)
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/72.0.3626.109 Safari/537.36'
        }
        params = {"organism": self.query_info}

        with open(cookie_file, 'rt') as f:
            cookies_dict = json.loads(f.read())
        cookies = requests.utils.cookiejar_from_dict(cookies_dict)
        # logger.debug(f'requests.get: url={url}, params={params}, header={headers}')
        response = requests.get(
            url,
            params=params,
            cookies=cookies,
            allow_redirects=True,
            stream=True,
            headers=headers,
            timeout=self.timeout,
        )
        try:
            response.raise_for_status()  # 如果响应的状态码不是200，将引发HTTPError异常
        except HTTPError:
            logger.critical(f'response status: {response.status_code}')
            logger.critical("Couldn't connect with server. Please check Internet connection and nretry.")

        xml_file = self.config.FILENAME_TEMPLATE_XML.format(self.query_info)
        with open(xml_file, "wb") as f:
            # 使用二进制写入模式（"wb"）来保存结果文件，因为response.content返回的是一个字节字符串
            logger.debug(f'successfully query, write xml @ {xml_file}')
            f.write(response.content)

    # step 05 parse xml info to update url
    def parse_xml(self):
        """
        Moves through the xml document <xml_file> and returns information
        about matches to elements in <DESIRED_CATEGORIES> if
        <filter_categories> is True, or all files otherwise

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        # Parse xml file for content to download

        xml_file = self.config.FILENAME_TEMPLATE_XML.format(self.query_info)

        logger.debug(f'start to parse_xml using parameter filter_categories = {self.filter_files}')
        display_cats = ['filename', 'url', 'size', 'label', 'sizeInBytes', 'timestamp', 'md5']
        # Choose between different XML parsers
        # Will only be used if --filter_files flag
        # if filter_files, user wants only those files in <_desired_categories>
        logger.debug('_xml_hunt...')
        found = self._xml_hunt(xml_file, logger=logger)
        found = self._format_found(found, self.filter_files, logger=logger)

        if not list(found.values()):
            return None
        logger.debug(f'successfully parsed xml @ {xml_file}')
        logger.debug('start to update file_list')

        category_id = 0

        for category, sub_cat in sorted(found.items()):
            if category not in self._desired_categories:
                category_id += 1
                self._desired_categories[category] = defaultdict(dict)
                self._desired_categories[category]["catID"] = category_id
            uid = 1
            for parent, children in sorted(sub_cat.items()):
                self._desired_categories[category]["results"][parent] = defaultdict(dict)
                results = self._desired_categories[category]["results"][parent]
                unique_children = self._uniqueify(children)
                for child in sorted(unique_children, key=lambda x: x['filename']):
                    try:
                        results[uid]
                    except KeyError:
                        results[uid] = {}
                    for dc in display_cats:
                        try:
                            results[uid][dc] = child[dc]
                        except KeyError:
                            continue
                    uid += 1

        logger.debug(
            f'successfully update self._desired_categories = {str(self._desired_categories)[:400]}...(omit)...')

        # Check if file has any categories of interest
        if not any(v["results"] for v in list(self._desired_categories.values())):
            logger.error(
                "no results found for '{}' in any of the following "
                "categories:\n---\n{}\n---"
                .format(self.query_info, "\n".join(self.config.info['categories']))
            )
            self._clean_exit(logger=logger)

    # step 06 download from url
    def download(self, logger=None):
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        self._decision_tree(logger=logger)

        self._urls_to_get = sorted(self._urls_to_get)
        filenames = [u.split("/")[-1] for u in self._urls_to_get]

        # logger.debug(f'self._desired_categories = {self._desired_categories}')
        file_sizes = self._get_sizes(self._desired_categories, sizes_by_url={})
        # logger.debug(f'file_sizes = {file_sizes}')
        total_size = sum(filter(None, [file_sizes[url] for url in self._urls_to_get]))
        size_string = self._byte_convert(total_size)
        num_files = len(self._urls_to_get)
        logger.info("Total download size for {} files: {}".format(num_files, size_string))

        if self.interactive:
            select = input("Continue? (y/n/b/p for yes/no/back/preview files): ").lower()
            if select == "p":
                while select == "p":
                    print('\n'.join(filenames))
                    select = input("Continue? (y/n/b/p for yes/no/back/preview files): ").lower()
            if select == "n":
                self._clean_exit("ABORTING DOWNLOAD", logger=logger)
            elif select == "b":
                logger.info('back to select files...')
                self.interactive = True
                logger.debug('!!! re-call self.download position1')
                self.download(logger=logger)
            elif select == 'y' or select == '':
                self._download_list(self._urls_to_get, logger=logger)
            else:
                logger.info('illegal selection, back to select files...')
                self.interactive = True
                logger.debug('!!! re-call self.download position2')
                self.download(logger=logger)
        else:  # non interactive
            if self.regex or self.all_get:
                self._download_list(self._urls_to_get, logger=logger)

        logger.info("Finished downloading {} files.".format(len(self._downloaded_files)))

        if self._failed_urls and self.interactive:
            n_broken = len(self._failed_urls)
            nretry_broken = input(
                "{} files failed to download; nretry them? (y/n): ".format(n_broken))
            if nretry_broken.lower() in ("yes", "y"):
                self._download_list(self._failed_urls, logger=logger)

        # Kindly offer to unpack files, if files remain after error check
        if self._downloaded_files and self.interactive:
            decompress = input(("Decompress all downloaded files? "
                                "(y/n/k, k=decompress and keep original): "))
            if decompress != "n":
                if decompress == "k":
                    keep_original = True
                else:
                    keep_original = False
                self._decompress_files(self._downloaded_files, keep_original, logger=logger)
                logger.info("Finished decompressing all files.")

        if self._failed_urls:
            """
            Write failed URLs to a local log file.

            """
            fail_log = self.FILENAME_TEMPLATE_LOG_FAIL.format(self.query_info)
            logger.info(
                "{} failed downloads logged to {}".format(len(self._failed_urls), fail_log))
            # write failed URLs to local file
            with open(fail_log, 'w') as f:
                f.write('\n'.join(self._failed_urls))
            failed_happen = True
        else:
            failed_happen = False

        # Clean up and exit
        # "cookies" file is always created
        exit_message = None
        remove_temp = True
        if self.interactive:
            keep_temp = input(f"Keep temporary files ('{self.config.FILENAME_TEMPLATE_XML.format(self.query_info)}' "
                              f"and '{self.config.FILENAME_TEMPLATE_COOKIE.format(self.query_info)}')? (y/n): ")
            if keep_temp.lower() in "y, yes":
                remove_temp = False
        elif failed_happen:  # failed files in non-interactive mode
            exit_message = (
                'Some files failed downloading')
            remove_temp = False

        exit_code = 1 if failed_happen else 0

        self._clean_exit(exit_message=exit_message, exit_code=exit_code, remove_temp=remove_temp, logger=logger)

    def _byte_convert(self, byte_size):
        """
        Converts a number of bytes to a more human-readable
        format.

        """
        # Calculate and display total size of selected data
        adjusted = byte_size / (1024 * 1024)  # bytes to MB
        if adjusted < 1:
            adjusted = byte_size / 1024
            unit = "KB"
        elif adjusted < 1024:
            unit = "MB"
        else:
            adjusted /= 1024
            unit = "GB"
        size_string = "{:.2f} {}".format(adjusted, unit)
        return size_string

    def _clean_exit(self, exit_message=None, exit_code=0, remove_temp=True, logger=None):
        """
        Perform a sys.exit() while removing temporary files and
        informing the user.

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        to_remove = [self.config.FILENAME_TEMPLATE_COOKIE.format(self.query_info)]

        # don't delete xml file if supplied by user
        if not self.xml and remove_temp is True:
            try:
                to_remove.append(self.config.FILENAME_TEMPLATE_XML.format(self.query_info))
            except NameError:
                pass
        for f in to_remove:
            try:
                os.remove(f)
            except OSError:
                continue
        if remove_temp is True:
            base_message = "Removing temp files and exiting"
        else:
            base_message = "Keeping temp files and exiting"
        if exit_message:
            logger.info(exit_message)
        logger.info(base_message)
        sys.exit(exit_code)

    def _check_for_xml(self, filename, logger=None):
        """
        Uses hex code at the beginning of a file to try to determine if it's an
        XML file or not. This seems to be occasionally necessary; if pulling
        files from JGI tape archives, the server may error out and provide an
        XML error document instead of the intended file. This function should
        return False on all downloaded files, although false positives have not
        been thoroughly investigated.

        Adapted from http://stackoverflow.com/a/13044946/3076552

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        xml_hex = "\x3c"  # hex code at beginning of XML file
        read_length = len(xml_hex)
        with open(filename) as f:
            try:
                file_start = f.read(read_length)

                if file_start.startswith(xml_hex):  # XML file
                    logger.debug('cheker for xml: XML file')
                    return True
                else:  # hopefully all other file types
                    logger.debug('cheker for xml: hopefully all other file types')
                    return False
            except UnicodeDecodeError:  # compressed file
                logger.debug('cheker for xml: compressed file')
                return False

    def _check_md5(self, filename, md5_hash, logger=None):
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        if not md5_hash:
            logger.warning(f"No MD5 hash listed for {filename}; skipping check")
            ret_val = True
        else:
            file_md5 = self._get_md5(filename)
            if file_md5 == md5_hash:
                logger.debug(f"MD5 hashes match for {filename} ({md5_hash})")
                ret_val = True
            else:
                logger.error(f"MD5 hash mismatch for {filename} (local: {file_md5}, remote: {md5_hash})")
                ret_val = False
        return ret_val

    def _check_sizeInBytes(self, filename, sizeInBytes, logger=None):
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        if not sizeInBytes:
            logger.info(f'No sizeInBytes listed for {filename}; skipping check')
            ret_val = True
        else:
            file_size_in_bytes = self._get_sizeInBytes(filename)
            if file_size_in_bytes == sizeInBytes:
                logger.debug(f'sizeInBytes match for {filename} ({sizeInBytes})')
                ret_val = True
            else:
                logger.error(
                    f'sizeInBytes mismatch for {filename} (local: {file_size_in_bytes}, remote: {sizeInBytes})')
                ret_val = False
        return ret_val

    def _decision_tree(self, logger=None):
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        # Decision tree depending on if non-interactive options given
        regex_filter = None
        user_choice = None
        self.interactive = True

        if self.all_get:  # 参数指定get all
            # non-interactive
            user_choice = "a"  # "all_get mode"
            self.interactive = False
        elif self.regex:  # 参数指定regex
            # non-interactive
            user_choice = "r"  # "regex mode"
            regex_filter = self.regex
            self.interactive = False
            # non-interactive
        elif self.log_fails is not None:  # 参数指定
            user_choice = "l"  # "re-download log_fails mode"
            self.interactive = False

        logger.debug(f'regex_filter = {regex_filter}')
        logger.debug(f'user_choice = {user_choice}')
        logger.debug(f'interactive_and_display_info = {self.interactive}')
        logger.debug(f"QUERY RESULTS FOR '{str(self._desired_categories)[:400]}...(omit)...'\n")

        for query_cat, v in sorted(iter(self._desired_categories.items()), key=lambda k_v: k_v[1]["catID"]):
            print_list = []

            if not v["results"]:
                continue

            catID = v["catID"]

            self._dict_to_get[catID] = {}
            print_list.append(f" {catID}: {query_cat} ".center(80, "="))
            results = v["results"]

            for sub_cat, items in sorted(iter(results.items()),
                                         key=lambda sub_cat_items: (sub_cat_items[0], sub_cat_items[1])):
                print_list.append("{}:".format(sub_cat))
                for index, i in sorted(items.items()):
                    url = i["url"]
                    self._dict_to_get[catID][index] = url
                    if "md5" in i:
                        self._url_to_validate[url]["md5"] = i["md5"]
                    # the following elif takes care of MD5 > sizeInBytes rank-order
                    # in downstream processing
                    if "sizeInBytes" in i:
                        self._url_to_validate[url]["sizeInBytes"] = int(i["sizeInBytes"])
                    print_index = " {}:[{}] ".format(str(catID), str(index))
                    date = self._format_timestamp(i["timestamp"])
                    date_string = "{:02d}/{}".format(date.tm_mon, date.tm_year)
                    size_date = "[{}|{}]".format(i["size"], date_string)
                    filename = i["filename"]
                    margin = 80 - (len(size_date) + len(print_index))
                    file_info = filename.ljust(margin, "-")
                    print_list.append("".join([print_index, file_info, size_date]))
            if self.interactive:
                print('\n'.join(print_list))
                print()  # padding

        if not user_choice:
            # Ask user which files to download from xml
            user_choice = self._get_user_choice(logger=logger)

        if user_choice == "r":  # regex-based filename matching
            regex_filter = self._get_regex()
        # special case for downloading all available files
        # or filtering with a regular expression
        if user_choice in (  # non interactive
                "a",  # all_get mode
                "r",  # regex mode
                "l",  # re-download log_fails mode
        ):
            for k, v in sorted(self._dict_to_get.items()):
                for u in v.values():
                    if regex_filter:
                        fn = re.search(r".+/([^\/]+$)", u).group(1)
                        match = regex_filter.search(fn)
                        if not match:
                            # 匹配失败
                            self._urls_to_get = set(self._urls_to_get)
                            self._urls_to_get.discard(u)
                            continue  # 进入下一个文件进行判断
                    elif user_choice == "l" and u not in self.log_fails:
                        continue  # 无需添加此文件，进入下一个文件进行判断
                    # final add to download
                    # logger.debug(f'self._urls_to_get = {self._urls_to_get}')
                    self._urls_to_get = set(self._urls_to_get)
                    self._urls_to_get.add(u)  # self._urls_to_get is a set()
        else:
            # interactive
            # 选2:3-5这种
            # Retrieve user-selected file urls from dict
            self._parse_selection(user_choice, logger=logger)  # update self._selections
            for k, v in sorted(self._selections.items()):
                for i in v:
                    self._urls_to_get.add(self._dict_to_get[k][i])

        # logger.debug('update self._dict_to_get, self._url_to_validate')
        # logger.debug(f'self._urls_to_get = {str(self._urls_to_get)[:1000]}...(omit)...')
        # logger.debug(f'self._dict_to_get = {str(self._dict_to_get)[:1000]}...(omit)...')

    def _decompress_files(self, local_file_list, keep_original=False, logger=None):
        """
        Decompresses list of files, and deletes compressed
        copies unless <keep_original> is True.

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        logger.debug('decompress_files...')
        for f in local_file_list:
            self._extract_file(f, keep_original, logger=logger)

    def _download_from_url(self, url, logger=None):
        """
        Attempts to download a file from JGI servers using cURL.

        Returns a tuple of (filename, cURL command used, success boolean)

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        url_no_prefix = url.replace("&amp;", "&")  # for query local xml info

        # get md5 value for this file from xml
        # logger.debug(f'self._url_to_validate = {self._url_to_validate}')
        md5_hash = self._url_to_validate[url_no_prefix].get("md5", None)  # local xml
        # get sizeInBytes (local xml) and file_size (internet) for size checking
        sizeInBytes = self._url_to_validate[url_no_prefix].get("sizeInBytes", None)  # local xml
        # counter
        error_counter = 0
        # download status
        loop = True
        success = False
        # filename to write in, parse from last of url
        filename = re.search('.+/(.+$)', url).group(1)

        while True:  # pass when loop == False
            if not loop:
                break

            if error_counter >= self.nretry:
                logger.critical(f'error_counter is reaching to max size -> self.nretry ={self.nretry}')
                # FUTURE, 支持了断点续传再改回来
                # logger.warning(
                #     'you can rerun your command and the unfinished file is needed for continuing download')
                success = False
                loop = False
                continue

            # temp_stemp_sizeize = 0  # 获取已写入的字节 temp_size
            if os.path.exists(filename):  # filename exists
                # check md5 for filename
                if not self._is_broken(filename, md5_hash=md5_hash, sizeInBytes=sizeInBytes, logger=logger):
                    logger.debug(f'File {filename} exists and passed md5 checking. Go on next task.')
                    success = True
                    loop = False
                    continue
                else:  # failed
                    if os.path.exists(filename + '.tmp'):  # filename.tmp exists at the sametime!
                        if not self._is_broken(filename + '.tmp', md5_hash=md5_hash, sizeInBytes=sizeInBytes,
                                               logger=logger):
                            logger.debug(f'File {filename}.tmp passed md5 checking! Renaming')
                            os.remove(filename)
                            os.rename(filename + '.tmp', filename)
                            logger.debug(f'File {filename} passed md5 checking. Go on next task.')
                            success = True
                            loop = False
                            continue
                        else:
                            # FUTURE: 由于JGI不支持断点续传，这个代码以后再用
                            # logger.error(f'File {filename} exists but is broken, file {filename}.tmp exists but is '
                            #              f'broken too. File {filename} is removed now.')
                            # os.remove(filename)
                            # temp_size = os.path.getsize(filename + '.tmp')
                            # logger.info(f'Try continuing download using {filename}.tmp. temp_size = {temp_size}')
                            # /FUTURE: 由于JGI不支持断点续传，这个代码以后再用
                            logger.warning(
                                f'File {filename} exists but is broken, file {filename}.tmp exists but is '
                                f'broken too.')
                            logger.info(
                                'Because JGI does not support continuing download now, the two files will be removed')
                            os.remove(filename)
                            logger.debug(f'File {filename} is removed now.')
                            os.remove(filename + '.tmp')
                            logger.debug(f'File {filename}.tmp is removed now.')
                            success = False
                            loop = True
                            continue
                    else:
                        logger.warning(f'File {filename} exists but is broken, file {filename}.tmp does not exists.'
                                       f'File {filename} is removed now.')
                        os.remove(filename)
                        logger.debug(f'Try to start a new download task.')
                        success = False
                        loop = True
                        continue
            else:  # filename don't exist
                if os.path.exists(filename + '.tmp'):  # filename.tmp exists
                    # 如果.tmp文件已存在，则打开文件并获取已写入的字节 temp_size
                    # print(f'found tmp file for {filename}, attempt to continuing download')
                    # temp_size = os.path.getsize(filename + '.tmp')

                    logger.warning(
                        'Because JGI does not support continuing download now, the tmp file will be removed')
                    os.remove(filename + '.tmp')
                    logger.debug(f'File {filename}.tmp is removed now.')
            # set requests
            # https://dabing1022.github.io/2016/12/24/%E8%81%8A%E4%B8%80%E8%81%8AHTTP%E7%9A%84Range,%20Content-Range/
            if not url.startswith('http'):
                url = f'{self.config.URL_JGI_MAIN}{url}'
            logger.debug(f'download aim file from url = {url}')
            cookie_file = self.config.FILENAME_TEMPLATE_COOKIE.format(self.query_info)
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/72.0.3626.109 Safari/537.36',
            }
            with open(cookie_file, 'rt') as f:
                cookies_dict = json.loads(f.read())
            cookies = requests.utils.cookiejar_from_dict(cookies_dict)

            # 想在Python中在不下载文件的情况下获取文件大小，可以使用requests库发送HEAD请求
            # HEAD请求不会下载文件，而是仅获取关于文件的一些元数据，如文件大小
            pre_response = requests.head(url, cookies=cookies, stream=True, headers=headers, timeout=self.timeout)
            # check response status
            try:
                pre_response.raise_for_status()  # 如果响应的状态码不是200，将引发HTTPError异常
            except HTTPError:
                error_counter += 1
                logger.warning(
                    'encounter HTTPError;'
                    f'error_counter = {error_counter}; pre_response.status_code = {pre_response.status_code}; '
                    'could not connect with server. Retry after 10s...'
                )
                time.sleep(10)
                logger.info('start next loop')
                loop = True
                pre_response.close()
                continue  # next try
            except Exception as e:
                error_counter += 1
                logger.warning(
                    f'encounter {e} (error);'
                    f'error_counter = {error_counter}; pre_response.status_code = {pre_response.status_code}; '
                    'unexpected error. Retry after 10s...'
                )
                time.sleep(10)
                logger.info('start next loop')
                success = False
                loop = True
                pre_response.close()
                continue  # next try

            # get file size, if pre_response.header do not have it, replace it from sizeInBytes in xml file
            remote_file_size = int(pre_response.headers.get("Content-Length", 0))
            remote_file_size = remote_file_size if remote_file_size != 0 else sizeInBytes
            pre_response.close()
            # FUTURE: 由于JGI不支持断点续传，这个代码以后再用
            # if temp_size > remote_file_size:
            #     # error
            #     logger.error('local size = {temp_size} > remote_file_size = {remote_file_size}, nretry...')
            #     error_counter += 1
            #     os.remove(filename + '.tmp')
            #     success = False
            #     loop = True
            #     pre_response.close()
            #     continue

            # core part for download
            # re-request use fixed headers
            # https://cizixs.com/2015/10/06/http-resume-download/
            # headers["Range"] = f"bytes={temp_size}-"
            # /FUTURE: 由于JGI不支持断点续传，这个代码以后再用
            file_response = requests.get(
                url,
                cookies=cookies,
                stream=True,  # 如果要下载大文件的话，就将steam设置为True，慢慢下载，而不是等整个文件下载完才返回。
                headers=headers,
                timeout=self.timeout,
            )  # 分段下载
            # FUTURE: 由于JGI不支持断点续传，代码以后再测试
            # print(file_response.headers)
            # print(file_response.status_code)
            # ################### core ########################
            with(open(filename + '.tmp', 'ab')) as f, tqdm(
                    desc=filename,
                    total=remote_file_size,
                    # initial=temp_size,  # FUTURE: 由于JGI不支持断点续传，这个代码以后再用
                    unit='iB',
                    unit_scale=True,
                    unit_divisor=1024,
            ) as pbar:
                try:
                    for data in file_response.iter_content(chunk_size=1024):
                        data_size = f.write(data)
                        pbar.update(data_size)
                    file_response.close()
                except InvalidChunkLength or ProtocolError or ChunkedEncodingError as e:
                    logger.error(f"Invalid chunk encoding {e}")
                    success = False
                    loop = True
                    file_response.close()
                    continue
            # ################### /core ########################
            # finish download
            # start check
            logger.debug(f'File {filename}.tmp seems done, checking md5 value')
            # 如果filename md5通过则删除tmp
            if not self._is_broken(filename + '.tmp', md5_hash=md5_hash, sizeInBytes=remote_file_size, logger=logger):
                logger.debug(f'File {filename}.tmp passed md5 checking! Renaming')
                os.rename(filename + '.tmp', filename)
                logger.debug(f'File {filename} passed md5 checking. Go on next task.')
                success = True
                loop = False
                file_response.close()
                continue  # next task
            else:
                error_counter += 1
                logger.info(f'File {filename}.tmp is broken! Removing and nretry')
                os.remove(filename + '.tmp')
                success = False
                loop = True
                file_response.close()
                continue  # next try
        return filename, success

    def _download_list(self, url_list, logger=None):
        """
        Attempts download command on a list of partial file
        URLs (completed by download_from_url()).

        Returns a list of successfully-downloaded files and a
        list of unsuccessful URLs

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        # Run curl commands to retrieve selected files
        # Make sure the URL formats conforms to the Genome Portal format

        broken_urls = []
        broken_files = []
        start_time = time.time()
        for url in url_list:
            current_time = time.time()
            # refresh the session cookie every 5 minutes
            if current_time - start_time > 300:
                logger.info('refresh the session cookie every 5 minutes')
                logger.info('run self.login.')
                self.login()
                logger.info('login succeed.')
                start_time = time.time()

            fn, success = self._download_from_url(url, logger=logger)

            if not success:
                logger.warning(f'File {fn} failed to download, appending to the file '
                               f'@ {self.config.FILENAME_TEMPLATE_LOG_FAIL.format(self.query_info)}. '
                               f'You can use parameter --log_fails to re-download these failed files'
                               f'after this run finish.')
                broken_urls.append(url)
                broken_files.append(fn)
            else:
                self._downloaded_files.append(fn)
        broken_urls = [
            u for u, f in zip(broken_urls, broken_files)
            if f not in self._downloaded_files
        ]

        return broken_urls

    def _extract_file(self, file_path, keep_compressed=False, logger=None):
        """
        Native Python file decompression for tar.gz and .gz files.

        TODO: implement .zip decompression

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        tar_pattern = "tar.gz$"  # matches tar.gz
        gz_pattern = "(?<!tar)\.gz$"  # excludes tar.gz
        endings_map = {"tar": (tarfile, "r:gz", ".tar.gz"),
                       "gz": (gzip, "rb", ".gz")
                       }
        relative_name = os.path.basename(file_path)
        if re.search(tar_pattern, file_path):
            logger.info(f'tar.gz file decompression for {file_path}')
            opener, mode, ext = endings_map["tar"]
            with opener.open(file_path) as f:
                file_count = len(f.getmembers())
                if file_count > 1:  # make sub-directory to unpack into
                    dir_name = relative_name.rstrip(ext)
                    try:
                        os.mkdir(dir_name)
                    except FileExistsError:
                        pass
                    destination = dir_name
                else:  # single file, extract into working directory
                    destination = "."

                def is_within_directory(directory, target):

                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)

                    prefix = os.path.commonprefix([abs_directory, abs_target])

                    return prefix == abs_directory

                def safe_extract(tar, path=".", members=None, *, numeric_owner=False):

                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not is_within_directory(path, member_path):
                            raise Exception("Attempted Path Traversal in Tar File")

                    tar.extractall(path, members, numeric_owner=numeric_owner)

                safe_extract(f, destination)
        elif re.search(gz_pattern, file_path):
            logger.info(f'gz file decompression for {file_path}')
            opener, mode, ext = endings_map["gz"]
            # out_name = file_path.rstrip(ext)
            out_name = relative_name.rstrip(ext)
            with opener.open(file_path) as f, open(out_name, "wb") as out:
                for l in f:
                    out.write(l)
        else:
            logger.info(f"Skipped decompression for '{file_path}'")
            return
        if not keep_compressed:
            logger.info('remove compressed files')
            os.remove(file_path)

    def _format_found(self, d, filter_found=False, logger=None):
        """
        Reformats the output from xml_hunt()

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        logger.debug('get categories from config (including possible user additions)')
        output = {}
        for p, c in sorted(d.items()):
            layers = [e for e in p.split(":") if e]
            if filter_found:
                if not any(cat in layers for cat in self.config.info["categories"]):
                    continue
            if len(layers) == 1:
                top = parent = layers[0]
            else:
                top = layers[-2]  # either -2 or -1 works well, != parent
                parent = layers[-1]  # either -2 or -1 works well, != top
            if top not in output:
                output[top] = defaultdict(dict)
            if parent not in output[top]:
                output[top][parent] = c
            else:
                output[top][parent].extend(c)
        return output

    def _format_timestamp(self, time_string):
        """
        Parses the timestamp string from an XML document
        of the form "Thu Feb 27 16:38:54 PST 2014"
        and returns a string of the form "2014".

        """
        # Remove platform-dependent timezone substring
        # of the general form "xxT"
        tz_pattern = re.compile("\s[A-Z]{3}\s")
        time_string = tz_pattern.sub(" ", time_string)
        # get the desired time info
        time_info = time.strptime(time_string, "%a %b %d %H:%M:%S %Y")
        # year = str(time_info.tm_year)
        return time_info

    def _get_sizes(self, d, sizes_by_url=None):
        """
        Builds a dictionary of url:sizes from
        output of get_file_list()

        """
        for k, v in d.items():
            if isinstance(v, dict):
                if "url" in v:
                    address = v["url"]
                    try:
                        size = int(v["sizeInBytes"])
                    except:
                        size = None
                    sizes_by_url[address] = size
                else:
                    self._get_sizes(v, sizes_by_url)
        return sizes_by_url

    def _get_md5(self, *fns, buffer_size=65536):
        hash = md5()
        for fn in fns:
            with open(fn, "rb") as f:
                while True:
                    data = f.read(buffer_size)
                    if not data:
                        break
                    hash.update(data)

        return hash.hexdigest()

    def _get_sizeInBytes(self, filename):
        try:
            file_size_in_bytes = os.path.getsize(filename)
        except:
            file_size_in_bytes = 0

        return file_size_in_bytes

    def _get_regex(self):
        """
        Get regex pattern from user, compile and return.

        """
        # TODO make this exit gracefully if user can't
        # manage to get a working regex
        compile_success = False
        pattern = ''
        while compile_success is False:
            if self.interactive and not self.regex:
                pattern = input("Regex pattern: ")
            else:
                pattern = self.regex
            try:
                pattern = re.compile(pattern)
                compile_success = True
                print(f'Your pattern = {pattern}')
            except:
                print("[!] ERROR: Regex pattern failed to compile.")

        return re.compile(pattern)

    def _get_user_choice(self, logger=None):
        """
        Get user file selection choice(s)

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        choice = input(
            "Enter file selection ('q' to quit, "
            "'u' to review syntax/usage, 'a' for all, "
            "'r' for regex-based filename matching):\n> ")
        logger.debug(f'choice = {choice}')
        if choice == "u":
            print()
            print(JGIDoc.select_blurb)
            print()
            return self._get_user_choice(logger=logger)
        elif choice == 'a':
            return choice
        elif choice == 'r':
            return choice
        elif choice == 'q':
            remove_temp = input("Remove index file? (y/n): ")
            remove_temp = remove_temp.lower() in ('y', 'yes', '')
            self._clean_exit(remove_temp=remove_temp, logger=logger)
        else:
            return choice

    def _is_broken(self, filename, min_size_bytes=20, md5_hash=None, sizeInBytes=None, logger=None):
        """
        Rudimentary check to see if a file appears to be broken.

        filename without ".tmp"
        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        if (
                not os.path.isfile(filename) or
                os.path.getsize(filename) < min_size_bytes or
                (self._check_for_xml(filename, logger=logger) and not filename.lower().endswith("xml")) or
                (not self._check_md5(filename, md5_hash, logger=logger) or
                 not self._check_sizeInBytes(filename, sizeInBytes, logger=logger))
        ):
            logger.debug('File is broken.')
            return True
        else:
            logger.debug('File is intact.')
            return False

    def _parse_selection(self, user_input, logger=None):
        """
        Parses the user choice string and returns a dictionary
        of categories (keys) and choices within each category
        (values).

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        parts = user_input.split(";")
        for p in parts:
            if len(p.split(":")) > 2 or p.count(':') == 0:
                logger.error(f"Can't parse desired input\n?-->'{p}'")
                user_choice = self._get_user_choice()
                self._parse_selection(user_choice, logger=logger)
            category, indices = p.split(":")
            category = int(category)
            self._selections[category] = []
            indices = indices.split(",")
            for i in indices:
                try:
                    self._selections[category].append(int(i))  # if it's already an integer
                except ValueError:
                    try:
                        start, stop = list(map(int, i.split("-")))
                    except:
                        logger.critical(f"can't parse desired input\n?-->'{i}'")
                        self._clean_exit(logger=logger)

                    add_range = list(range(start, stop + 1))

                    for e in add_range:
                        self._selections[category].append(e)

    def _xml_hunt(self, xml_file, logger=None):
        """
        Gets list of all XML entries with "filename" attribute,
        and returns a dictionary of the file attributes keyed
        by a ":"-joined string of parent names.

        """
        if not logger:
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
        logger.debug('xml_hunt...')
        root = ET.iterparse(xml_file, events=("start", "end"))
        parents = []
        matches = {}
        try:
            for event, element in root:
                if element.tag not in ["folder", "file"]:  # skip topmost categories
                    continue
                if element.tag == "folder":
                    if event == "start":  # add to parents
                        parents.append(element.attrib["name"])
                    elif event == "end":  # strip from parents
                        del parents[-1]
                    continue
                if event == "start" and element.tag == "file":
                    parent_string = ":".join(parents)
                    try:
                        matches[parent_string].append(element.attrib)
                    except KeyError:
                        matches[parent_string] = [element.attrib]
        except ET.ParseError as e:
            logger.critical('Parsing xml file failed! Maybe the username / password is wrong? '
                            'Try to use --overwrite_conf parameter to re-login')
            sys.exit(f'Exit with error {e}')
        return matches

    def _uniqueify(self, children):
        """
        Takes a list of child XML elements (dicts of attribs) as
        returns a filtered list of only unique filenames for a given
        month/year timestamp (e.g. duplicates are allowed if month/year
        is different).

        """
        unique = {}
        for child in children:
            try:
                fn = child['filename']
                date = self._format_timestamp(child['timestamp'])
                date_string = (date.tm_mon, date.tm_year)
                uid = (fn, date_string)
            except KeyError:
                continue
            if fn not in unique:
                unique[uid] = child
            else:
                existing = unique[uid].get('fileType', None)
                if existing == 'Unknown':
                    existing = None
                current = child.get('fileType', None)
                if current == 'Unknown':
                    current = None
                if current is not None and existing is None:
                    unique[uid] = child

        return unique.values()
