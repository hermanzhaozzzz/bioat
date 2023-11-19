import sys
import os
import re
import subprocess

import xml.etree.ElementTree as ET

from bioat import get_logger
from bioat.lib.libjgi import JGIDoc
from bioat.lib.libpath import HOME

__module_name__ = 'bioat.metatools'


class MetaTools:
    """Metagenome toolbox."""

    def JGI_query(
            self,
            organism_abbreviation: str | None = None,
            xml: str | None = None,
            configure: bool = False,
            syntax_help: bool = False,
            filter_files: bool = False,
            usage: bool = False,
            retry: int = 0,
            load_failed: str | None = None,
            regex: str | None = None,
            get_all: bool = False,
            log_level: str = 'INFO'
    ):
        """JGI_query, a tool for downloading files from JGI-IMG database.

        This script will list and retrieve files from JGI using the curl API.
        It will return a list of all files available for download for a given query organism.

        The source code is adapted from https://github.com/glarue/jgi-query

        :param organism_abbreviation: organism name formatted per JGI's abbreviation.
            For example, 'Nematostella vectensis' is abbreviated by JGI as 'Nemve1'.
            The appropriate abbreviation may be found by searching for the organism on JGI;
            the name used in the URL of the 'Info' page for that organism is the correct abbreviation.
            The full URL may also be used for this argument.
        :param xml: specify a local xml file for the query instead of retrieving a new copy from JGI
        :param configure: initiate configuration dialog to overwrite existing user/password configuration (interactive)
        :param syntax_help: syntax_help (no interactive)
        :param filter_files: filter organism results by config categories instead of reporting all files listed by JGI
            for the query (work in progress)
        :param usage: print verbose usage information and exit (no interactive)
        :param retry: number of times to retry downloading files with errors (0 to skip such files)
        :param load_failed: retry downloading from URLs listed in log file
        :param regex: (no interactive) Regex pattern to use to auto-select and download files
        :param get_all: (no interactive) Auto-select and download all files for query
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        run_docs = any([syntax_help, usage])
        # ------------------------------------------------------------------->>>>>>>>>>
        # doc mode
        # ------------------------------------------------------------------->>>>>>>>>>
        if run_docs:
            logger.info('user choose no interactive mode')
            # Check if user wants query help
            if syntax_help:
                logger.info(f'\n[syntax_help]:\n{JGIDoc.select_blurb}')
            if usage:
                logger.info(f'\n[usage]:\n{JGIDoc.usage_example_blurb}')
            sys.exit('Done. exit.')
        # finally exit
        # ------------------------------------------------------------------->>>>>>>>>>
        # no interactive mode
        # ------------------------------------------------------------------->>>>>>>>>>


        # ------------------------------------------------------------------->>>>>>>>>>
        # interactive mode
        # ------------------------------------------------------------------->>>>>>>>>>


        # CONFIG

        # Get script location info
        script_path = os.path.realpath(sys.argv[0])
        script_home = os.path.dirname(script_path)

        # Config should be in same directory as script
        config_filename = "jgi-query.config"
        config_filepath = script_home + "/{}".format(config_filename)

        # Does config file exist?
        if os.path.isfile(config_filepath) and not configure:  # use config file
            config_info = read_config(config_filepath)
        else:  # no config present or configure flag used; run config dialog
            config_info = get_user_info()
            config_info["categories"] = DEFAULT_CATEGORIES
            make_config(config_filepath, config_info)

        # Get user information for sign-on
        USER = config_info["user"]
        PASSWORD = config_info["password"]

        # Set curl login string using user and password as per https://goo.gl/oppZ2a

        # Old syntax
        # LOGIN_STRING = ("curl https://signon.jgi.doe.gov/signon/create --data-ascii "
        #                 "login={}\&password={} -b cookies -c cookies > "
        #                 "/dev/null".format(USER, PASSWORD))

        # New syntax
        LOGIN_STRING = (
            # "curl 'https://signon-old.jgi.doe.gov/signon/create' "
            "curl 'https://signon.jgi.doe.gov/signon/create' "
            "--data-urlencode 'login={}' "
            "--data-urlencode 'password={}' "
            "-s "  # suppress status output
            "-c cookies > /dev/null"
            .format(USER, PASSWORD)
        )

        LOCAL_XML = False

        # pull info from log file if provided
        if load_failed:
            logfile = load_failed
            org_input = os.path.basename(logfile).split('.')[0]
            with open(logfile) as f:
                load_failed = f.read().splitlines()
            # logfile = args.load_failed
            # print("Reading URLs from \'{}\'".format(logfile))
            # downloaded, failed = retry_from_failed(LOGIN_STRING, logfile)
            # clean_exit("All files in log attempted.")
        else:
            org_input = organism_abbreviation

        if not org_input:
            if configure:
                sys.exit("Configuration complete. Script may now be used to query JGI. "
                         "Exiting now.")
            elif xml and xml != 1:
                # Use org_input because is already checked further down
                # and avoids re-writing this whole block
                org_input = get_org_name(xml)
                if not org_input:
                    sys.exit("No organism specified. Exiting now.")
            else:
                sys.exit("No organism specified. Exiting now.")
        org_regex = re.compile(r'\.jgi.+\.(?:gov|org).*\/(.+)\/(?!\/)')
        try:  # see if it's in address form
            # organism = re.search("\.jgi.+\.(?:gov|org)/(.+)/", org_input).group(1)
            organism = org_regex.search(org_input).group(1)
        except AttributeError:  # not in address form, assume string is organism name
            organism = org_input

        # URL where remote XML file should be, if it exists
        org_url = ("https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?"
                   "organism={}".format(organism))

        # Get xml index of files, using existing local file or curl API
        if xml:
            LOCAL_XML = True  # global referenced by clean_exit()
            xml_arg = xml
            if xml_arg == 1:  # --xml flag used without argument
                xml_index_filename = "{}_jgi_index.xml".format(organism)
            else:
                xml_index_filename = xml_arg
            print(
                "Retrieving information from JGI for query "
                "'{}' using local file '{}'\n".format(organism, xml_index_filename))
        else:  # fetch XML file from JGI
            xml_index_filename = "{}_jgi_index.xml".format(organism)

            # Old syntax
            # xml_address = ("curl {} -b cookies -c cookies > {}"
            #                .format(org_url, xml_index_filename))

            # New syntax
            xml_address = ("curl '{}' -L -b cookies > {}"
                           .format(org_url, xml_index_filename))
            try:  # fails if unable to contact server
                subprocess.check_output(LOGIN_STRING, shell=True)
            except subprocess.CalledProcessError as error:
                clean_exit("Couldn't connect with server. Please check Internet "
                           "connection and retry.")
            print(
                "Retrieving information from JGI for query '{}' using command "
                "'{}'\n".format(organism, xml_address))
            subprocess.run(xml_address, shell=True)
            print()  # padding

        # Parse xml file for content to download
        xml_root = None
        if os.path.getsize(xml_index_filename) == 0:  # happens if user and/or pw wrong
            clean_exit("Invalid username/password combination (or other issue).\n"
                       "Restart script with flag '-c' to reconfigure credentials.")
        try:
            xml_in = ET.ElementTree(file=xml_index_filename)
            xml_root = xml_in.getroot()
        except ET.ParseError:  # organism not found/xml file contains errors
            clean_exit("Cannot parse XML file or no organism match found.\n"
                       "Ensure remote file exists and has content at the "
                       "following address:\n{}".format(org_url))

        # Get categories from config (including possible user additions)
        # Will only be used if --filter_files flag
        DESIRED_CATEGORIES = config_info["categories"]

        # Choose between different XML parsers
        # if args.filter_files, user wants only those files in <desired_categories>
        file_list = get_file_list(xml_index_filename, filter_categories=filter_files)

        # Check if file has any categories of interest
        if not any(v["results"] for v in list(file_list.values())):
            print(("ERROR: no results found for '{}' in any of the following "
                   "categories:\n---\n{}\n---"
                   .format(organism, "\n".join(DESIRED_CATEGORIES))))
            clean_exit()

        # Decision tree depending on if non-interactive options given
        regex_filter = None
        user_choice = None
        display_info = True
        if get_all:
            user_choice = "a"
            display_info = False
        elif regex:
            user_choice = "r"
            regex_filter = regex
            display_info = False
        elif load_failed is not None:
            user_choice = "l"
            display_info = False

        url_dict, url_to_validate = print_data(file_list, organism, display=display_info)

        if not user_choice:
            # Ask user which files to download from xml
            user_choice = get_user_choice()
            if user_choice == 'r':
                regex_filter = get_regex()

        urls_to_get = set()

        # special case for downloading all available files
        # or filtering with a regular expression
        if user_choice in ("a", "r", "l"):
            for k, v in sorted(url_dict.items()):
                for u in v.values():
                    if regex_filter:
                        fn = re.search(r".+/([^\/]+$)", u).group(1)
                        match = regex_filter.search(fn)
                        if not match:
                            continue
                    elif user_choice == "l" and u not in load_failed:
                        continue
                    urls_to_get.add(u)
        else:
            # Retrieve user-selected file urls from dict
            ids_dict = parse_selection(user_choice)
            for k, v in sorted(ids_dict.items()):
                for i in v:
                    urls_to_get.add(url_dict[k][i])

        # Calculate and display total size of selected data
        urls_to_get = sorted(urls_to_get)
        filenames = [u.split("/")[-1] for u in urls_to_get]
        file_sizes = get_sizes(file_list, sizes_by_url={})
        total_size = sum(filter(None, [file_sizes[url] for url in urls_to_get]))
        size_string = byte_convert(total_size)
        num_files = len(urls_to_get)
        print(("Total download size for {} files: {}".format(num_files, size_string)))
        if interactive:
            download = input("Continue? (y/n/[p]review files): ").lower()
            if download == "p":
                while download == "p":
                    print('\n'.join(filenames))
                    download = input("Continue with download? (y/n/[p]review files): ").lower()
            if download != "y":
                clean_exit("ABORTING DOWNLOAD")

        downloaded_files, failed_urls = download_list(
            urls_to_get, url_to_validate=url_to_validate, retries=retry)

        print("Finished downloading {} files.".format(len(downloaded_files)))

        if failed_urls and interactive:
            n_broken = len(failed_urls)
            retry_broken = input(
                "{} files failed to download; retry them? (y/n): ".format(n_broken))
            if retry_broken.lower() in ("yes", "y"):
                downloaded_files, failed_urls = download_list(
                    failed_urls, url_to_validate=url_to_validate, retries=1)

        # Kindly offer to unpack files, if files remain after error check
        if downloaded_files and interactive:
            decompress = input(("Decompress all downloaded files? "
                                "(y/n/k=decompress and keep original): "))
            if decompress != "n":
                if decompress == "k":
                    keep_original = True
                else:
                    keep_original = False
                decompress_files(downloaded_files, keep_original)
                print("Finished decompressing all files.")

        # TODO either offer to delete or append ".error" to local broken files
        if failed_urls:
            log_failed(organism, failed_urls)
            SOME_FAILED = True
        else:
            SOME_FAILED = False

        # Clean up and exit
        # "cookies" file is always created
        exit_message = None
        remove_temp = True
        if interactive:
            keep_temp = input("Keep temporary files ('{}' and 'cookies')? (y/n): "
                              .format(xml_index_filename))
            if keep_temp.lower() in "y, yes":
                remove_temp = False
        elif SOME_FAILED:  # failed files in non-interactive mode
            exit_message = (
                'Some files failed downloading')
            remove_temp = False

        exit_code = 1 if SOME_FAILED else 0

        clean_exit(exit_message=exit_message, exit_code=exit_code, remove_temp=remove_temp)
