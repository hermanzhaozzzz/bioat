import os
import sys
from datetime import datetime
from time import sleep

import matplotlib.pyplot as plt
import pandas as pd
import requests
from bs4 import BeautifulSoup

from bioat import BioatInvalidOptionError
from bioat.lib.libpatentseq import query_patent
from bioat.logger import get_logger

__module_name__ = "bioat.searchtools"


class SearchTools:
    """Search toolbox."""

    def __init__(self):
        pass

    def google_scholar(
            self,
            keyword: str = "CRISPR",
            sort_by: str = "CitePerYear",
            n_results: int = 100,
            output: str | None = None,
            save_table: bool = True,
            plot: bool = False,
            start_year: int = None,
            end_year: int = datetime.now().year,
            log_level: str = "WARNING",
    ):
        """Return a table with a list of publications from google scholar, sort_by cit/year.

        This code creates a database with a list of publications data from Google
        Scholar.
        The data acquired from GS is Title, Citations, Links and Rank.
        It is useful for finding relevant papers by sorting by the number of citations.

        As output this program will plot the number of citations in the Y axis and the
        rank of the result in the X axis. It also, optionally, export the database to
        a .csv file.

        :param keyword: Keyword to be searched. Use double quote followed by simple quote to search for an exact
                keyword. Example: "'exact keyword'"
        :param sort_by: Column to be sorted by.
                i.e., it will be sorted by the number of citations.
                If you want to sort by citations per year, use --sort_by "CitePerYear"
                Or --sort_by "Citations" to sort by citations totally.
        :param n_results: Number of articles to search on Google Scholar.
                Default is 100. (be careful with robot checking if value is too high)
        :param output: Path of the exported table file. format can be 'csv', 'tsv', 'xls' or 'xlsx'(default).
        :param save_table: By default, results are going to be exported to a csv file.
                Select this option to just print results but not store them
        :param plot: Use this flag in order to plot the results with the original rank in the x-axis
                and the number of citaions in the y-axis.
                Default is False
        :param start_year: Start year when searching. Default is None
        :param end_year: End year when searching. Default is current year
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        logger = get_logger(
            level=log_level,
            module_name=__module_name__,
            func_name="google_scholar",
        )

        def get_citations(content):
            out = 0
            for char in range(0, len(content)):
                if content[char: char + 9] == "Cited by ":
                    init = char + 9
                    for end in range(init + 1, init + 6):
                        if content[end] == "<":
                            break
                    out = content[init:end]
            return int(out)

        def get_year(content):
            for char in range(0, len(content)):
                if content[char] == "-":
                    out = content[char - 5: char - 1]
            if not out.isdigit():
                out = 0
            return int(out)

        def setup_driver():
            try:
                from selenium import webdriver
                from selenium.common.exceptions import StaleElementReferenceException
                from selenium.webdriver.chrome.options import Options
            except ImportError as e:
                logger.error(e)
                logger.error(
                    "Please install Selenium using `pip install selenium`"
                )
                sys.exit(0)

            logger.info("Loading...")
            options = Options()
            options.add_argument("disable-infobars")
            driver = webdriver.Chrome(options=options)
            return driver

        def get_author(content):
            for char in range(0, len(content)):
                if content[char] == "-":
                    out = content[2: char - 1]
                    break
            return out

        def get_element(driver, xpath, attempts=5, _count=0):
            """Safe get_element method with multiple attempts"""
            try:
                element = driver.find_element_by_xpath(xpath)
                return element
            except Exception:
                if _count < attempts:
                    sleep(1)
                    get_element(driver, xpath, attempts=attempts, _count=_count + 1)
                else:
                    logger.warning("Element not found")

        def get_content_with_selenium(url):
            if "driver" not in globals():
                global driver
                driver = setup_driver()
            driver.get(url)

            # Get element from page
            el = get_element(driver, "/html/body")
            c = el.get_attribute("innerHTML")

            if any(kw in el.text for kw in ROBOT_KW):
                logger.info(
                    "Solve captcha manually and press enter here to continue..."
                )
                el = get_element(driver, "/html/body")
                c = el.get_attribute("innerHTML")

            return c.encode("utf-8")

        # Websession Parameters
        GSCHOLAR_URL = (
            "https://scholar.google.com/scholar?start={}&q={}&hl=en&as_sdt=0,5"
        )
        YEAR_RANGE = ""  # &as_ylo={start_year}&as_yhi={end_year}'
        # GSCHOLAR_URL_YEAR = GSCHOLAR_URL + YEAR_RANGE
        STARTYEAR_URL = "&as_ylo={}"
        ENDYEAR_URL = "&as_yhi={}"
        ROBOT_KW = ["unusual traffic from your computer network", "not a robot"]

        # Create main URL based on command line arguments
        GSCHOLAR_MAIN_URL = (
            GSCHOLAR_URL + STARTYEAR_URL.format(start_year)
            if start_year
            else GSCHOLAR_URL
        )
        GSCHOLAR_MAIN_URL = (
            GSCHOLAR_MAIN_URL + ENDYEAR_URL.format(end_year)
            if end_year
            else GSCHOLAR_MAIN_URL
        )

        if log_level == "DEBUG":
            GSCHOLAR_MAIN_URL = (
                    "https://web.archive.org/web/20210314203256/" + GSCHOLAR_URL
            )

        # Start new session
        session = requests.Session()
        # headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36'}

        # Variables
        links = []
        title = []
        citations = []
        year = []
        author = []
        venue = []
        publisher = []
        rank = [0]

        # Get content from number_of_results URLs
        for n in range(0, n_results, 10):
            # if start_year is None:
            url = GSCHOLAR_MAIN_URL.format(str(n), keyword.replace(" ", "+"))
            logger = get_logger(
                level=log_level,
                module_name=__module_name__,
                func_name="google_scholar",
            )

            logger.debug("Opening URL:", url)
            # else:
            #    url=GSCHOLAR_URL_YEAR.format(str(n), keyword.replace(' ','+'), start_year=start_year, end_year=end_year)

            logger.info("Loading next {} results".format(n + 10))
            page = session.get(url)  # , headers=headers)
            c = page.content
            if any(kw in c.decode("ISO-8859-1") for kw in ROBOT_KW):
                logger.info(
                    "Robot checking detected, handling with selenium (if installed)"
                )
                try:
                    c = get_content_with_selenium(url)
                except Exception as e:
                    logger.error("No success. The following error was raised:")
                    logger.error(e)

            # Create parser
            soup = BeautifulSoup(c, "html.parser", from_encoding="utf-8")

            # Get stuff
            mydivs = soup.findAll("div", {"class": "gs_or"})

            for div in mydivs:
                try:
                    links.append(div.find("h3").find("a").get("href"))
                except:  # catch *all* exceptions
                    links.append("Look manually at: " + url)

                try:
                    title.append(div.find("h3").find("a").text)
                except:
                    title.append("Could not catch title")

                try:
                    citations.append(get_citations(str(div.format_string)))
                except:
                    logger.warning(
                        "Number of citations not found for {}. Appending 0".format(
                            title[-1]
                        )
                    )
                    citations.append(0)

                try:
                    year.append(get_year(div.find("div", {"class": "gs_a"}).text))
                except:
                    logger.warning(
                        "Year not found for {}, appending 0".format(title[-1])
                    )
                    year.append(0)

                try:
                    author.append(get_author(div.find("div", {"class": "gs_a"}).text))
                except:
                    author.append("Author not found")

                try:
                    publisher.append(
                        div.find("div", {"class": "gs_a"}).text.split("-")[-1]
                    )
                except:
                    publisher.append("Publisher not found")

                try:
                    venue.append(
                        " ".join(
                            div.find("div", {"class": "gs_a"})
                            .text.split("-")[-2]
                            .split(",")[:-1]
                        )
                    )
                except:
                    venue.append("Venue not fount")

                rank.append(rank[-1] + 1)

            # Delay
            sleep(0.5)

        # Create a dataset and sort by the number of citations
        data = pd.DataFrame(
            list(zip(author, title, citations, year, publisher, venue, links)),
            index=rank[1:],
            columns=[
                "Author",
                "Title",
                "Citations",
                "Year",
                "Publisher",
                "Venue",
                "Source",
            ],
        )
        data.index.name = "Rank"

        # Add columns with number of citations per year
        data["CitePerYear"] = data["Citations"] / (end_year + 1 - data["Year"])
        data["CitePerYear"] = data["CitePerYear"].round(0).astype(int)
        data = data[['Author', 'Citations', "CitePerYear", 'Year', 'Venue', 'Title', 'Publisher']]

        # Sort by the selected columns, if exists
        try:
            data_ranked = data.sort_values(by=sort_by, ascending=False)
        except Exception as e:
            print(
                "Column name to be sorted not found. Sorting by the number of citations..."
            )
            data_ranked = data.sort_values(by="Citations", ascending=False)
            print(e)
        # fix index
        data_ranked.reset_index(drop=True, inplace=True)
        # Print data
        print(data_ranked)

        # Plot by citation number
        if plot:
            plt.plot(rank[1:], citations, "*")
            plt.ylabel("Number of Citations")
            plt.xlabel("Rank of the keyword on Google Scholar")
            plt.title("Keyword: " + keyword)
            plt.show()

        # Save results
        if save_table:
            if not output:
                fpath_csv = os.path.join('.', keyword.replace(" ", "_") + ".xlsx")
                data_ranked.to_excel(fpath_csv, index=False)
            elif output.endswith(".csv"):
                data_ranked.to_csv(output, encoding="utf-8")
            elif output.endswith(".tsv"):
                data_ranked.to_csv(output, sep="\t", index=False)
            elif output.endswith(".xlsx"):
                data_ranked.to_excel(output, index=False)
            elif output.endswith(".xls"):
                data_ranked.to_excel(output, index=False)
            else:
                raise BioatInvalidOptionError("Output format not supported")

    def query_patent(
            self,
            seq: str,
            query_name: str | None = None,
            username: str | None = None,
            password: str | None = None,
            via_proxy: str | None = None,  # proxy_server = "socks5://127.0.0.1:8235",
            output: str | None = None,
            nobrowser: bool = True,
            retry: int = 3,
            local_browser: str | None = None,
            rm_fail_cookie: bool = False,
            log_level: str = "INFO",
    ):
        """Return a table with a list of patent blast hit from lens.org

        :param seq: protein sequence, e.g. MCRISQQKK
        :param query_name: queryName in output table
        :param username: ORCID username(usually a mail)
        :param password: ORCID password
        :param via_proxy: like http://127.0.0.1:8234 socks5://127.0.0.1:8235
        :param output: output table.csv/csv.gz
        :param nobrowser: wether or not to open browser for DEBUG
        :param retry: max retry times
        :param local_browser: local firefox browser executable file path
        :param rm_fail_cookie: remove cookies from local if query fail, default is False
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        logger = get_logger(
            level=log_level,
            module_name=__module_name__,
            func_name="query_patent",
        )
        logger.info("Run query patent sequence from lens.org...")
        query_patent(
            seq=seq,
            seq_header=query_name,
            username=username,
            password=password,
            proxy_server=via_proxy,
            output=output,
            headless=nobrowser,
            nretry=retry,
            local_browser=local_browser,
            rm_fail_cookie=rm_fail_cookie,
            log_level=log_level,
        )
