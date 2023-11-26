import tarfile
import gzip
import time
import os
import sys
import re
import json
import subprocess
import textwrap
import requests
from collections import defaultdict
from hashlib import md5
from bioat.lib.libpath import HOME
from bioat.logger import get_logger
import xml.etree.ElementTree as ET
from requests import HTTPError
from requests.utils import cookiejar_from_dict

# https://blog.csdn.net/m0_49960764/article/details/119460979

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
        This script will retrieve files from JGI using the cURL api. It will
        return a list of possible files for downloading.
    
        * This script depends upon cURL - it can be downloaded here:
        http://curl.haxx.se/
    
        # USAGE ///////////////////////////////////////////////////////////////////////
        $ jgi-query.py [<jgi_address>, <jgi_abbreviation>] [[-xml [<your_xml>]], -f]
    
        To get <jgi_address>, go to: http://genome.jgi.doe.gov/ and search for your
        species of interest. Click through until you are at the "Info" page. For
        \x1B[3mNematostella vectensis\x1B[23m, the appropriate page is
        "http://genome.jgi.doe.gov/Nemve1/Nemve1.info.html".
    
        To query using only the name simply requires the specific JGI organism
        abbreviation, as referenced in the full url.
    
        For the above example, the proper input syntax for this script would be:
    
        $ jgi-query.py http://genome.jgi.doe.gov/Nemve1/Nemve1.info.html
    
                                 -or-
    
        $ jgi-query.py Nemve1
    
        If you already have the XML file for the query in the directory, you may use
        the --xml flag to avoid redownloading it (particularly useful if querying
        large, top-level groups with many sub-species, such as "fungi"):
    
        $ jgi-query.py --xml <your_xml_index>
    
        If the XML filename is omitted when using the --xml flag, it is assumed that
        the XML file is named '<jgi_abbreviation>_jgi_index.xml'. In such cases, the
        organism name is required.
        # /USAGE //////////////////////////////////////////////////////////////////////
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
    # url home for jgi-img database
    URL_JGI_MAIN = "https://genome.jgi.doe.gov"
    # url login
    URL_JGI_LOGIN = 'https://signon.jgi.doe.gov/signon/create'
    # url fetch xml
    # https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=GuaBasSediGuay16
    URL_JGI_FETCH_XML = "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory"
    # remote XML file should be
    URL_TEMPLATE_XML_FILE = "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism={}"
    # see https://genome.jgi.doe.gov/portal/help/download.jsf#/api for CMD_TEMPLATE_LOGIN_JGI
    CMD_TEMPLATE_LOGIN_JGI = (
        "curl 'https://signon.jgi.doe.gov/signon/create' "
        "--data-urlencode 'login={}' --data-urlencode 'password={}'-c {} > /dev/null"
    )
    FILENAME_TEMPLATE_COOKIE = "jgi-xml-query.cookie_{}.cookie"
    FILENAME_TEMPLATE_XML = "jgi-xml-query.result_{}.xml"
    FILENAME_TEMPLATE_LOG_FAIL = "jgi-xml-query.failed_{}.log"
    FILENAME_CONFIG_PATH = os.path.join(HOME, '.bioat', 'JGI', "account.conf")

    def __init__(self, overwrite_conf: bool = False, log_level='INFO'):
        logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        self.log_level = log_level
        self.info: dict = {"user": None, "password": None, "categories": None}

        logger.debug('start to set / load user info...')
        if overwrite_conf:
            # 只要指定重写参数就重写，无需判断其他
            self.input_user_info()  # 从终端手动输入JGI用户信息
            self.save_config()  # 保存用户信息
            # 然后直接退出
            logger.info("Configuration complete.\n"
                        "Script may now be used to query JGI.\n"
                        "You need re-run your command without parameter `--overwrite_conf`\n")
            # self.load_config()  # 从文件加载用户信息
        else:
            if os.path.isfile(self.FILENAME_CONFIG_PATH) and os.path.getsize(self.FILENAME_CONFIG_PATH) > 0:
                # 直接加载配置信息
                self.load_config()
            else:
                # 配置信息保存文件异常
                logger.warning(f'errors occur when loading file: {self.FILENAME_CONFIG_PATH}!\n'
                               f'try to get user info by manual inputting...')
                self.input_user_info()  # 从终端手动输入JGI用户信息
                self.save_config()  # 保存用户信息
                self.load_config()  # 从文件加载用户信息
        logger.debug('set / load user info success!')

    def load_config(self):
        """
        Reads "user", "password" and "categories" entries
        from config file.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        logger.debug('loading JGI account info...')

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
            logger.critical(f"Config file present ({self.FILENAME_CONFIG_PATH}), but user and/or password not found.")
            sys.exit(1)

        logger.debug('loading JGI account info done')

    def input_user_info(self):
        """
        Dialog with user to gather user information for
        use with the curl query. Returns a dict.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

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
            logger.warning('user and password can not be None!')
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
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        logger.debug('saving JGI account info...')
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

        logger.debug('saving JGI account info done')


class JGIOperator:
    def __init__(
            self,
            # pick one from three
            query_info: str | None = None,
            xml: str | None = None,
            failed_log: str | None = None,
            # runtime params
            retry: int = 5,
            timeout: int = -1,
            regex: str | None = None,
            get_all: bool = False,
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
        self.failed_log = failed_log
        self.retry = retry
        self.timeout = timeout
        self.regex = regex
        self.get_all = get_all
        self.overwrite_conf = overwrite_conf
        self.filter_files = filter_files  # TODO 搞清楚用法
        self.syntax_help = syntax_help
        self.usage = usage
        self.log_level = log_level
        # self
        self._dict_to_get = {}
        self._urls_to_get = set()
        self._url_to_validate = defaultdict(dict)
        self._selections = {}
        self._downloaded_files = []
        self._failed_urls = []
        self._desired_categories = dict()
        # From other obj #
        # load configs; auto check if you need overwrite user info or not
        self.config = JGIConfig(overwrite_conf=overwrite_conf, log_level=log_level)
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

    # step 02 update self.query_info and self.failed_log and followed by self.login
    def parse_input(self):
        """Checker for input.

        A checker for this rule: ONLY ONE of the parameters ('query_info', 'xml' and 'failed_log') can be specified;
        Parse query_info/xml/failed_log to update self.query_info and self.failed_log.
        After this method, you should call self.login method and then self.query method
        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        logger.debug("update JGIOperator obj.query_info and obj.failed_log using parameters:"
                     " 'query_info', 'xml' and 'failed_log'")

        if sum(map(bool, [self.query_info, self.xml, self.failed_log])) != 1:
            # 这三个参数只可以指定一个！
            logger.error("ONLY ONE of the parameters ('query_info', 'xml' and 'failed_log') can be specified!")
            sys.exit('Done. exit.')
        else:
            logger.debug("pass param checker: 'query_info', 'xml' and 'failed_log'")

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
        elif self.failed_log:
            # load failed info from log file if provided
            # filename see: class JGIConfig.FILENAME_TEMPLATE_LOG_FAIL
            self.query_info = os.path.basename(self.failed_log).split('.')[1].split('_')[0]

            with open(self.failed_log, 'rt') as f:
                self.failed_log = f.read().splitlines()
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
                # DEBUG
                print(f'_org_line = {_org_line}')
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
        response = s.post(url, headers=headers, data=data)

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
            logger.info(f'successfully login, write cookie @ {cookie_file}')
        else:
            logger.critical("Couldn't connect with server. Please check Internet connection and retry.")
            self._clean_exit()

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
        )
        try:
            response.raise_for_status()  # 如果响应的状态码不是200，将引发HTTPError异常
        except HTTPError:
            logger.critical(f'response status: {response.status_code}')
            logger.critical("Couldn't connect with server. Please check Internet connection and retry.")

        xml_file = self.config.FILENAME_TEMPLATE_XML.format(self.query_info)
        with open(xml_file, "wb") as f:
            # 使用二进制写入模式（"wb"）来保存结果文件，因为response.content返回的是一个字节字符串
            logger.info(f'successfully query, write xml @ {xml_file}')
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
        found = self._xml_hunt(xml_file)
        found = self._format_found(found, self.filter_files)

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

        logger.debug(f'successfully update self._desired_categories = {self._desired_categories}')

        # Check if file has any categories of interest
        if not any(v["results"] for v in list(self._desired_categories.values())):
            logger.error(
                "no results found for '{}' in any of the following "
                "categories:\n---\n{}\n---"
                .format(self.query_info, "\n".join(self.config.info['categories']))
            )
            self._clean_exit()

        # Decision tree depending on if non-interactive options given
        regex_filter = None
        user_choice = None
        interactive_and_display_info = True

        if self.get_all:
            # non-interactive
            user_choice = "get_all mode"
            interactive_and_display_info = False
        elif self.regex:
            # non-interactive
            user_choice = "regex mode"
            regex_filter = self.regex
            interactive_and_display_info = False
            # non-interactive
        elif self.failed_log is not None:
            user_choice = "re-download failed_log mode"
            interactive_and_display_info = False
        # """
        # Prints info from dictionary data in a specific format.
        # Returns a dict with url information for every file
        # in desired categories, as well as a dict with md5 information for
        # each file (keyed by file URL).
        #
        # """
        logger.debug(f'regex_filter = {regex_filter}')
        logger.debug(f'user_choice = {user_choice}')
        logger.debug(f'interactive_and_display_info = {interactive_and_display_info}')
        logger.debug(f"\nQUERY RESULTS FOR '{self._desired_categories}'\n")

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
                    integrity_tag = ""
                    url = i["url"]
                    self._dict_to_get[catID][index] = url
                    if "md5" in i:
                        self._url_to_validate[url]["md5"] = i["md5"]
                    # the following elif takes care of MD5 > sizeInBytes rank-order
                    # in downstream processing
                    elif "sizeInBytes" in i:
                        self._url_to_validate[url]["sizeInBytes"] = int(i["sizeInBytes"])
                    print_index = " {}:[{}] ".format(str(catID), str(index))
                    date = self._format_timestamp(i["timestamp"])
                    date_string = "{:02d}/{}".format(date.tm_mon, date.tm_year)
                    size_date = "[{}|{}]".format(i["size"], date_string)
                    filename = i["filename"]
                    margin = 80 - (len(size_date) + len(print_index))
                    file_info = filename.ljust(margin, "-")
                    print_list.append("".join([print_index, file_info, size_date]))
            if interactive_and_display_info:
                print('\n'.join(print_list))
                print()  # padding

        if not user_choice:
            # Ask user which files to download from xml
            user_choice = self._get_user_choice()
            if user_choice == "regex mode":
                regex_filter = self._get_regex()

        # special case for downloading all available files
        # or filtering with a regular expression
        if user_choice in ("get_all mode", "regex mode", "re-download failed_log mode"):
            for k, v in sorted(self._dict_to_get.items()):
                for u in v.values():
                    if regex_filter:
                        fn = re.search(r".+/([^\/]+$)", u).group(1)
                        match = regex_filter.search(fn)
                        if not match:
                            continue
                    elif user_choice == "l" and u not in self.failed_log:
                        continue
                    self._urls_to_get.add(u)
        else:
            # Retrieve user-selected file urls from dict
            self._parse_selection(user_choice)
            for k, v in sorted(self._selections.items()):
                for i in v:
                    self._urls_to_get.add(self._dict_to_get[k][i])

        logger.debug('update self._dict_to_get, self._url_to_validate')
        logger.debug(f'self._urls_to_get = {self._urls_to_get}')
        logger.debug(f'self._dict_to_get = {self._dict_to_get}')

    # step 06 download from url
    def download(self):
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        self._urls_to_get = sorted(self._urls_to_get)
        filenames = [u.split("/")[-1] for u in self._urls_to_get]
        file_sizes = self._get_sizes(self._desired_categories, sizes_by_url={})
        total_size = sum(filter(None, [file_sizes[url] for url in self._urls_to_get]))
        size_string = self._byte_convert(total_size)
        num_files = len(self._urls_to_get)
        print(("Total download size for {} files: {}".format(num_files, size_string)))
        # if interactive:
        download = input("Continue? (y/n/[p]review files): ").lower()
        if download == "p":
            while download == "p":
                print('\n'.join(filenames))
                download = input("Continue with download? (y/n/[p]review files): ").lower()
        if download != "y":
            self._clean_exit("ABORTING DOWNLOAD")

        self._download_list(self._urls_to_get)

        print("Finished downloading {} files.".format(len(self._downloaded_files)))

        if self._failed_urls and self.interactive:
            n_broken = len(self._failed_urls)
            retry_broken = input(
                "{} files failed to download; retry them? (y/n): ".format(n_broken))
            if retry_broken.lower() in ("yes", "y"):
                self._download_list(self._failed_urls)

        # Kindly offer to unpack files, if files remain after error check
        if self._downloaded_files and self.interactive:
            decompress = input(("Decompress all downloaded files? "
                                "(y/n/k=decompress and keep original): "))
            if decompress != "n":
                if decompress == "k":
                    keep_original = True
                else:
                    keep_original = False
                self._decompress_files(self._downloaded_files, keep_original)
                print("Finished decompressing all files.")

        # TODO either offer to delete or append ".error" to local broken files
        if self._failed_urls:
            """
            Write failed URLs to a local log file.

            """
            fail_log = self.FILENAME_TEMPLATE_LOG_FAIL.format(self.query_info)
            print(
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

        self._clean_exit(exit_message=exit_message, exit_code=exit_code, remove_temp=remove_temp)

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

    def _clean_exit(self, exit_message=None, exit_code=0, remove_temp=True):
        """
        Perform a sys.exit() while removing temporary files and
        informing the user.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
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

    def _check_md5(self, filename, md5_hash, print_message=True):
        if not md5_hash:
            message = "INFO: No MD5 hash listed for {}; skipping check".format(filename)
            ret_val = True
        else:
            file_md5 = self._get_md5(filename)
            if file_md5 == md5_hash:
                message = (
                    "SUCCESS: MD5 hashes match for {} ({})".format(filename, md5_hash))
                ret_val = True
            else:
                message = ("ERROR: MD5 hash mismatch for {} (local: {}, remote: {})"
                           .format(filename, file_md5, md5_hash))
                ret_val = False

        if print_message is True:
            print(message)

        return ret_val

    def _check_sizeInBytes(self, filename, sizeInBytes, print_message=True):
        if not sizeInBytes:
            message = "INFO: No sizeInBytes listed for {}; skipping check".format(filename)
            ret_val = True
        else:
            file_sizeInBytes = self._get_sizeInBytes(filename)
            if file_sizeInBytes == sizeInBytes:
                message = (
                    "SUCCESS: sizeInBytes match for {} ({})".format(filename, sizeInBytes))
                ret_val = True
            else:
                message = ("ERROR: sizeInBytes mismatch for {} (local: {}, remote: {})"
                           .format(filename, file_sizeInBytes, sizeInBytes))
                ret_val = False

        if print_message is True:
            print(message)

        return ret_val

    def _decompress_files(self, local_file_list, keep_original=False):
        """
        Decompresses list of files, and deletes compressed
        copies unless <keep_original> is True.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        for f in local_file_list:
            self._extract_file(f, keep_original)

    def _download_from_url(self, url, timeout=None, min_file_bytes=20):
        """
        Attempts to download a file from JGI servers using cURL.

        Returns a tuple of (filename, cURL command used, success boolean)

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        success = True
        md5_hash = self._url_to_validate[url].get("md5", None)
        sizeInBytes = self._url_to_validate[url].get("sizeInBytes", None)
        url = url.replace("&amp;", "&")

        filename = re.search('.+/(.+$)', url).group(1)
        timeout = '' if timeout is None else f'-m {timeout}'

        cmd_download = (f"curl {timeout} '{self.config.URL_JGI_MAIN}{url}' "
                        f"-b {self.config.FILENAME_TEMPLATE_COOKIE.format(self.query_info)} "
                        f"> {filename}")

        if not self._is_broken(filename, md5_hash=md5_hash, sizeInBytes=sizeInBytes):
            success = True
            print("Skipping existing file {}".format(filename))
        else:
            print("Downloading '{}' using command:\n{}"
                  .format(filename, cmd_download))
            # The next line doesn't appear to be needed to refresh the cookies.
            #    subprocess.call(login, shell=True)
            logger.info(cmd_download)
            status = subprocess.run(cmd_download, shell=True).returncode

            logger.debug('!!!!!!cmd_download = {cmd_download}')

            if status != 0 or self._is_broken(
                    filename, min_file_bytes, md5_hash=md5_hash, sizeInBytes=sizeInBytes
            ):
                success = False
                if self.retry > 0:
                    # success = False
                    # this may be needed if initial download fails
                    alt_cmd = cmd_download.replace(
                        "blocking=true", "blocking=false")
                    current_retry = 1
                    while current_retry <= self.retry:
                        if current_retry % 2 == 1:
                            retry_cmd = alt_cmd
                        else:
                            retry_cmd = cmd_download
                        print(
                            "Trying '{}' again due to download error ({}/{}):\n{}"
                            .format(filename, current_retry, self.retry, retry_cmd)
                        )
                        logger.info(retry_cmd)
                        status = subprocess.run(retry_cmd, shell=True).returncode
                        if status == 0 and not self._is_broken(
                                filename, min_file_bytes, md5_hash=md5_hash, sizeInBytes=sizeInBytes
                        ):
                            success = True
                            break
                        current_retry += 1
                        time.sleep(10)

        return filename, cmd_download, success

    def _download_list(self, url_list, timeout=None):
        """
        Attempts download command on a list of partial file
        URLs (completed by download_from_url()).

        Returns a list of successfully-downloaded files and a
        list of unsuccessful URLs

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        # Run curl commands to retrieve selected files
        # Make sure the URL formats conforms to the Genome Portal format

        broken_urls = []
        broken_files = []

        cookie_file = self.config.FILENAME_TEMPLATE_COOKIE.format(self.query_info)
        cmd_login_jgi = JGIConfig.CMD_TEMPLATE_LOGIN_JGI.format(
            self.config.info["user"], self.config.info["password"], cookie_file
        )
        subprocess.run(cmd_login_jgi, shell=True)
        start_time = time.time()
        for url in url_list:
            current_time = time.time()
            # refresh the session cookie every 5 minutes
            if current_time - start_time > 300:
                subprocess.run(cmd_login_jgi, shell=True)
                start_time = time.time()
            fn, cmd, success = self._download_from_url(url, timeout=timeout)
            if not success:
                broken_urls.append(url)
                broken_files.append(fn)
            else:
                self._downloaded_files.append(fn)
        # in cases where multiple files with same name are present and any of them
        # succeed, we can remove corresponding URLs from the list of broken URLs
        # (otherwise, they would just overwrite one another).
        # TODO we could also rename any files with identical names, although then
        # we would need to differentiate between files with different content and
        # files that are just broken versions of the same file...
        broken_urls = [
            u for u, f in zip(broken_urls, broken_files)
            if f not in self._downloaded_files
        ]

        return broken_urls

    def _extract_file(self, file_path, keep_compressed=False):
        """
        Native Python file decompression for tar.gz and .gz files.

        TODO: implement .zip decompression

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        tar_pattern = "tar.gz$"  # matches tar.gz
        gz_pattern = "(?<!tar)\.gz$"  # excludes tar.gz
        endings_map = {"tar": (tarfile, "r:gz", ".tar.gz"),
                       "gz": (gzip, "rb", ".gz")
                       }
        relative_name = os.path.basename(file_path)
        if re.search(tar_pattern, file_path):
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
            opener, mode, ext = endings_map["gz"]
            # out_name = file_path.rstrip(ext)
            out_name = relative_name.rstrip(ext)
            with opener.open(file_path) as f, open(out_name, "wb") as out:
                for l in f:
                    out.write(l)
        else:
            print("Skipped decompression for '{}'"
                  .format(file_path))
            return
        if not keep_compressed:
            os.remove(file_path)

    def _format_found(self, d, filter_found=False):
        """
        Reformats the output from xml_hunt()

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
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
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        # Remove platform-dependent timezone substring
        # of the general form "xxT"
        tz_pattern = re.compile("\s[A-Z]{3}\s")
        time_string = tz_pattern.sub(" ", time_string)
        logger.debug('get the desired time info')
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
            file_sizeInBytes = os.path.getsize(filename)
        except:
            file_sizeInBytes = 0

        return file_sizeInBytes

    def _get_regex(self):
        """
        Get regex pattern from user, compile and return.

        """
        # TODO make this exit gracefully if user can't
        # manage to get a working regex
        compile_success = False
        while compile_success is False:
            pattern = input("Regex pattern: ")
            try:
                pattern = re.compile(pattern)
                compile_success = True
            except:
                print("[!] ERROR: Regex pattern failed to compile.")

        return re.compile(pattern)

    def _get_user_choice(self):
        """
        Get user file selection choice(s)

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        choice = input(
            "Enter file selection ('q' to quit, "
            "'usage' to review syntax, 'a' for all, "
            "'r' for regex-based filename matching):\n> ")
        if choice == "usage":
            print()
            print(JGIDoc.select_blurb)
            print()
            return self._get_user_choice()
        elif choice.lower() in ("q", "quit", "exit"):
            remove_temp = input("Remove index file? (y/n): ")
            remove_temp = remove_temp.lower() in ('y', 'yes', '')
            self._clean_exit(remove_temp=remove_temp)
        else:
            return choice

    def _is_broken(self, filename, min_size_bytes=20, md5_hash=None, sizeInBytes=None):
        """
        Rudimentary check to see if a file appears to be broken.

        """

        def _check_for_xml(filename):
            """
            Uses hex code at the beginning of a file to try to determine if it's an
            XML file or not. This seems to be occasionally necessary; if pulling
            files from JGI tape archives, the server may error out and provide an
            XML error document instead of the intended file. This function should
            return False on all downloaded files, although false positives have not
            been thoroughly investigated.

            Adapted from http://stackoverflow.com/a/13044946/3076552

            """
            logger = get_logger(level=self.log_level, module_name=__module_name__,
                                func_name=sys._getframe().f_code.co_name)
            logger.debug('cheker for xml')
            xml_hex = "\x3c"  # hex code at beginning of XML file
            read_length = len(xml_hex)

            with open(filename) as f:
                try:
                    file_start = f.read(read_length)

                    if file_start.startswith(xml_hex):  # XML file
                        logger.debug('cheker for xml: XML file')
                        return True
                    else:  # hopefully all other file types
                        logger.warning('cheker for xml: hopefully all other file types')
                        return False
                except UnicodeDecodeError:  # compressed file
                    logger.critical('cheker for xml: compressed file')
                    return False

        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        if (
                not os.path.isfile(filename) or
                os.path.getsize(filename) < min_size_bytes or
                (_check_for_xml(filename) and not filename.lower().endswith("xml")) or
                (not self._check_md5(filename, md5_hash) or
                 not self._check_sizeInBytes(filename, sizeInBytes))
        ):
            logger.warning('is broken!')
            return True
        else:
            logger.error('is intact!')
            return False

    def _parse_selection(self, user_input):
        """
        Parses the user choice string and returns a dictionary
        of categories (keys) and choices within each category
        (values).

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        parts = user_input.split(";")
        for p in parts:
            if len(p.split(":")) > 2:
                self._clean_exit("FATAL ERROR: can't parse desired input\n?-->'{}'".format(p))
            category, indices = p.split(":")
            category = int(category)
            self._selections[category] = []
            cat_list = self._selections[category]
            indices = indices.split(",")
            for i in indices:
                try:
                    cat_list.append(int(i))  # if it's already an integer
                except ValueError:
                    try:
                        start, stop = list(map(int, i.split("-")))
                    except:
                        logger.critical(f"can't parse desired input\n?-->'{i}'")
                        self._clean_exit()

                    add_range = list(range(start, stop + 1))

                    for e in add_range:
                        cat_list.append(e)

    def _xml_hunt(self, xml_file):
        """
        Gets list of all XML entries with "filename" attribute,
        and returns a dictionary of the file attributes keyed
        by a ":"-joined string of parent names.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
        logger.debug('xml_hunt...')
        root = ET.iterparse(xml_file, events=("start", "end"))
        parents = []
        matches = {}
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
