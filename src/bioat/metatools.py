
from bioat.lib.libjgi import JGIOperator
from bioat.logger import get_logger

__module_name__ = "bioat.metatools"


class MetaTools:
    """Metagenome toolbox."""

    def JGI_query(
        self,
        query_info: str | None = None,
        xml: str | None = None,
        log_fails: str | None = None,
        nretry: int = 4,
        timeout: int = 60,
        regex: str | None = None,
        all_get: bool = False,
        overwrite_conf: bool = False,
        filter_files: bool = False,
        proxy_pool: str | None = None,
        just_query_xml: bool = False,
        syntax_help: bool = False,
        usage: bool = False,
        log_level: str = "INFO",
    ):
        """JGI_query: Tool for downloading files from the JGI-IMG database.

        This function lists and retrieves files from JGI using the curl API and
        returns a list of all files available for download for a given query organism.

        The source code is adapted from https://github.com/glarue/jgi-query.

        Args:
            query_info (str | None):
                Organism name formatted per JGI's abbreviation.
                Example: 'Nematostella vectensis' is abbreviated by JGI as 'Nemve1'.
                The correct abbreviation can be found by searching for the organism on JGI;
                the name used in the URL of the 'Info' page for that organism is
                the correct abbreviation. The full URL may also be used for this argument.
            xml (str | None):
                Specify a local XML file for the query instead of retrieving
                a new copy from JGI.
            log_fails (str | None):
                Log file containing URLs to retry downloading from in case of failure.
            nretry (int):
                Number of times to retry downloading files with errors.
                Use 0 to skip such files.
            timeout (int):
                Timeout (in seconds) for downloading. Set to -1 to disable.
            regex (str | None):
                Regex pattern to use for auto-selecting and downloading files
                without interaction.
            all_get (bool):
                If True, auto-select and download all files for the query
                without interaction.
            overwrite_conf (bool):
                If True, initiate configuration dialog to overwrite
                existing user/password configuration.
            filter_files (bool):
                Under development. Filter organism results by config categories
                instead of reporting all files listed by JGI for the query.
            proxy_pool (str | None):
                URL for the proxy pool, e.g., http://abc.com:port. See
                https://github.com/hermanzhaozzzz/proxy_pool.
            just_query_xml (bool):
                Set True if you just want to save the XML file.
            syntax_help (bool):
                If True, provide syntax help in doc mode.
            usage (bool):
                If True, print verbose usage information and exit.
            log_level (str):
                Set the logging level. Options include 'CRITICAL', 'ERROR',
                'WARNING', 'INFO', 'DEBUG', 'NOTSET'.
        """
        # Load or create JGI account info and
        # auto check if user info needs to be overwritten.

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
            func_name="JGI_query",
        )
        logger.debug("run query")
        operator.query()
        logger.debug("parse xml to json")
        operator.parse_xml()
        logger.debug(
            "start to download; calculate and display total size of selected data"
        )
        operator.download()
