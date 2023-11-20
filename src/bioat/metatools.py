import sys
import os
import re

from bioat import get_logger
from bioat.lib.libjgi import JGIDoc, JGIConfig, JGIOperator

__module_name__ = 'bioat.metatools'


class MetaTools:
    """Metagenome toolbox."""

    def JGI_query(
            self,
            query_info: str | None = None,
            xml: str | None = None,
            overwrite_conf: bool = False,
            syntax_help: bool = False,
            filter_files: bool = False,
            usage: bool = False,
            retry: int = 0,
            load_failed_log: str | None = None,
            regex: str | None = None,
            get_all: bool = False,
            log_level: str = 'WARNING'
    ):
        """JGI_query, a tool for downloading files from JGI-IMG database.

        This script will list and retrieve files from JGI using the curl API.
        It will return a list of all files available for download for a given query organism.

        The source code is adapted from https://github.com/glarue/jgi-query

        :param query_info: organism name formatted per JGI's abbreviation.
            For example, 'Nematostella vectensis' is abbreviated by JGI as 'Nemve1'.
            The appropriate abbreviation may be found by searching for the organism on JGI;
            the name used in the URL of the 'Info' page for that organism is the correct abbreviation.
            The full URL may also be used for this argument.
        :param xml: specify a local xml file for the query instead of retrieving a new copy from JGI
        :param overwrite_conf: initiate configuration dialog to overwrite existing user/password configuration (interactive)
        :param syntax_help: (doc mode) syntax_help
        :param filter_files: filter organism results by config categories instead of reporting all files listed by JGI
            for the query (work in progress)
        :param usage: (doc mode) print verbose usage information and exit
        :param retry: number of times to retry downloading files with errors (0 to skip such files)
        :param load_failed_log: retry downloading from URLs listed in log file
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
            logger.info('doc mode is selected')
            # Check if user wants query help
            if syntax_help:
                print(f'\n[syntax_help]:\n{JGIDoc.select_blurb}')
            if usage:
                print(f'\n[usage]:\n{JGIDoc.usage_example_blurb}')
            sys.exit('Done. exit.')
        # finally exit
        # ------------------------------------------------------------------->>>>>>>>>>
        # load or create JGI account info
        # ------------------------------------------------------------------->>>>>>>>>>
        # check if it needs exit
        if overwrite_conf and not any([query_info, xml, load_failed_log]):
            # exit
            logger.info("Configuration complete. Script may now be used to query JGI. ")
            sys.exit('Done. exit.')
        else:
            # go on
            logger.info("Select to run interactive mode")

        # param checker
        if sum(map(bool, [query_info, xml, load_failed_log])) != 1:
            logger.error("ONLY ONE of the parameters ('query_info', 'xml' and 'load_failed_log') can be specified!")
            sys.exit('Done. exit.')
        else:
            logger.debug("pass param checker: 'query_info', 'xml' and 'load_failed_log'")
        # ------------------------------------------------------------------->>>>>>>>>>
        # interactive mode
        # ------------------------------------------------------------------->>>>>>>>>>
        interactive = True
        # start to run main pipeline
        failed_logs = None

        if load_failed_log:
            # pull failed info from log file if provided
            query_info = os.path.basename(load_failed_log).split('.')[0]  # TODO 错误路径，记得修正

            with open(load_failed_log) as f:
                failed_logs = f.read().splitlines()
        elif xml:  # Get xml index of files, using existing local file or curl API
            pass  # TODO ????
        elif query_info:
            try:
                # if query_info is a URL
                logger.debug('attempt to parse query_info as a url')
                query_regex = re.compile(r'\.jgi.+\.(?:gov|org).*\/(.+)\/(?!\/)')
                query_info = query_regex.search(query_info).group(1)
            except AttributeError:
                # if query_info is an organism name
                logger.debug('not a url, try to define query_info as an organism name abbreviation')
                pass
        else:
            logger.error("one of the parameters 'query_info' and 'xml' should be specified")
            sys.exit('Done. exit.')

        # get JGI operator obj
        operator = JGIOperator(
            query_info=query_info,
            overwrite_conf=overwrite_conf,
            xml=xml,
            get_all=get_all,
            regex=regex,
            failed_logs=failed_logs,
            retry=retry,
            interactive=interactive,
        )
        operator.login()  # run login
        operator.query()  # run query
        # Moves through the xml document <xml_file> and returns information
        #         about matches to elements in <DESIRED_CATEGORIES> if
        #         <filter_categories> is True, or all files otherwise
        desired_categories = operator.get_file_list(filter_categories=filter_files)
        # Check if file has any categories of interest
        operator.category_picker(desired_categories)
        # start to download
        # Calculate and display total size of selected data
        operator.download(desired_categories)

