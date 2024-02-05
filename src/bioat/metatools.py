import sys
from bioat import get_logger
from bioat.lib.libjgi import JGIOperator

__module_name__ = "bioat.metatools"


class MetaTools:
    """Metagenome toolbox."""

    def JGI_query(
        self,
        # pick one from three
        query_info: str | None = None,
        xml: str | None = None,
        log_fails: str | None = None,
        # runtime params
        nretry: int = 4,
        timeout: int = 60,
        regex: str | None = None,
        all_get: bool = False,
        overwrite_conf: bool = False,
        filter_files: bool = False,
        proxy_pool: str | None = None,
        just_query_xml: bool = False,
        # doc helper
        syntax_help: bool = False,
        usage: bool = False,
        # log
        log_level: str = "INFO",
    ):
        """JGI_query, a tool for downloading files from JGI-IMG database.

        This script will list and retrieve files from JGI using the curl API.
        It will return a list of all files available for download for a given query organism.

        The source code is adapted from https://github.com/glarue/jgi-query

        :param query_info: (input) organism name formatted per JGI's abbreviation.
            For example, 'Nematostella vectensis' is abbreviated by JGI as 'Nemve1'.
            The appropriate abbreviation may be found by searching for the organism on JGI;
            the name used in the URL of the 'Info' page for that organism is the correct abbreviation.
            The full URL may also be used for this argument.
        :param xml: (input) specify a local xml file for the query instead of retrieving a new copy from JGI
        :param log_fails: (input) nretry downloading from URLs listed in log file
        :param nretry: number of times to nretry downloading files with errors (0 to skip such files)
        :param timeout: timeout (seconds) for downloading, set -1 to disable this
        :param regex: (no interactive) Regex pattern to use to auto-select and download files
        :param all_get: (no interactive) Auto-select and download all files for query
        :param overwrite_conf: (interactive) initiate configuration dialog to overwrite existing user/password
            configuration
        :param filter_files: (work in progress) filter organism results by config categories instead of reporting all
            files listed by JGI for the query
        :param proxy_pool: str | None = None, http://abc.com:port, see https://github.com/hermanzhaozzzz/proxy_pool
        :param just_query_xml: bool, set True if you just want save XML file.
        :param syntax_help: (doc mode) syntax_help
        :param usage: (doc mode) print verbose usage information and exit
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        # load or create JGI account info &
        # auto check if you need overwrite user info
        operator = JGIOperator(
            query_info=query_info,
            xml=xml,
            log_fails=log_fails,
            nretry=nretry,
            timeout=timeout,
            regex=regex,
            all_get=all_get,
            overwrite_conf=overwrite_conf,
            filter_files=filter_files,
            proxy_pool=proxy_pool,
            just_query_xml=just_query_xml,
            syntax_help=syntax_help,
            usage=usage,
            log_level=log_level,
        )
        logger = get_logger(
            level=log_level,
            module_name=__module_name__,
            func_name=sys._getframe().f_code.co_name,
        )
        logger.debug("run query")
        operator.query()
        logger.debug("parse xml to json")
        operator.parse_xml()
        logger.debug(
            "start to download; calculate and display total size of selected data"
        )
        operator.download()

    def mafft(self):
        pass
