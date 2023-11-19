import os
import sys
from time import sleep
from datetime import datetime
from bioat import get_logger

import requests
from bs4 import BeautifulSoup
import matplotlib.pyplot as plt
import pandas as pd

__module_name__ = 'bioat.searchtools'


class SearchTools:
    """Search toolbox."""

    def __init__(self):
        pass

    def google_scholar(
            self,
            keyword: str = 'CRISPR',
            sort_by: str = 'cit/year',
            n_results: int = 100,
            csv_path: str = '.',
            save_table: bool = True,
            plot: bool = False,
            start_year: int = None,
            end_year: int = datetime.now().year,
            log_level: str = 'WARNING'
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
        :param sort_by: Column to be sorted by. Default is by the columns "Citations",
                i.e., it will be sorted by the number of citations.
                If you want to sort by citations per year, use --sort_by "cit/year"
                Or --sort_by "Citations" to sort by citations totally
        :param n_results: Number of articles to search on Google Scholar.
                Default is 100. (be careful with robot checking if value is too high)
        :param csv_path: Path to save the exported csv file. By default it is the current folder
        :param save_table: By default results are going to be exported to a csv file.
                Select this option to just print results but not store them
        :param plot: Use this flag in order to plot the results with the original rank in the x-axis
                and the number of citaions in the y-axis.
                Default is False
        :param start_year: Start year when searching. Default is None
        :param end_year: End year when searching. Default is current year
        :param log_level: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'
        """
        logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        def get_citations(content):
            out = 0
            for char in range(0, len(content)):
                if content[char:char + 9] == 'Cited by ':
                    init = char + 9
                    for end in range(init + 1, init + 6):
                        if content[end] == '<':
                            break
                    out = content[init:end]
            return int(out)

        def get_year(content):
            for char in range(0, len(content)):
                if content[char] == '-':
                    out = content[char - 5:char - 1]
            if not out.isdigit():
                out = 0
            return int(out)

        def setup_driver():
            try:
                from selenium import webdriver
                from selenium.webdriver.chrome.options import Options
                from selenium.common.exceptions import StaleElementReferenceException
            except Exception as e:
                logger.error(e)
                logger.error("Please install Selenium and chrome webdriver for manual checking of captchas")

            logger.info('Loading...')
            chrome_options = Options()
            chrome_options.add_argument("disable-infobars")
            driver = webdriver.Chrome(chrome_options=chrome_options)
            return driver

        def get_author(content):
            for char in range(0, len(content)):
                if content[char] == '-':
                    out = content[2:char - 1]
                    break
            return out

        def get_element(driver, xpath, attempts=5, _count=0):
            '''Safe get_element method with multiple attempts'''
            try:
                element = driver.find_element_by_xpath(xpath)
                return element
            except Exception as e:
                if _count < attempts:
                    sleep(1)
                    get_element(driver, xpath, attempts=attempts, _count=_count + 1)
                else:
                    logger.warning("Element not found")

        def get_content_with_selenium(url):
            if 'driver' not in globals():
                global driver
                driver = setup_driver()
            driver.get(url)

            # Get element from page
            el = get_element(driver, "/html/body")
            c = el.get_attribute('innerHTML')

            if any(kw in el.text for kw in ROBOT_KW):
                logger.info("Solve captcha manually and press enter here to continue...")
                el = get_element(driver, "/html/body")
                c = el.get_attribute('innerHTML')

            return c.encode('utf-8')

        MAX_CSV_FNAME = 255
        # Websession Parameters
        GSCHOLAR_URL = 'https://scholar.google.com/scholar?start={}&q={}&hl=en&as_sdt=0,5'
        YEAR_RANGE = ''  # &as_ylo={start_year}&as_yhi={end_year}'
        # GSCHOLAR_URL_YEAR = GSCHOLAR_URL + YEAR_RANGE
        STARTYEAR_URL = '&as_ylo={}'
        ENDYEAR_URL = '&as_yhi={}'
        ROBOT_KW = ['unusual traffic from your computer network', 'not a robot']

        # Create main URL based on command line arguments
        GSCHOLAR_MAIN_URL = GSCHOLAR_URL + STARTYEAR_URL.format(start_year) if start_year else GSCHOLAR_URL
        GSCHOLAR_MAIN_URL = GSCHOLAR_MAIN_URL + ENDYEAR_URL.format(end_year) if end_year else GSCHOLAR_MAIN_URL

        if log_level == 'DEBUG':
            GSCHOLAR_MAIN_URL = 'https://web.archive.org/web/20210314203256/' + GSCHOLAR_URL

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
            url = GSCHOLAR_MAIN_URL.format(str(n), keyword.replace(' ', '+'))
            logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

            logger.debug("Opening URL:", url)
            # else:
            #    url=GSCHOLAR_URL_YEAR.format(str(n), keyword.replace(' ','+'), start_year=start_year, end_year=end_year)

            logger.info("Loading next {} results".format(n + 10))
            page = session.get(url)  # , headers=headers)
            c = page.content
            if any(kw in c.decode('ISO-8859-1') for kw in ROBOT_KW):
                logger.info("Robot checking detected, handling with selenium (if installed)")
                try:
                    c = get_content_with_selenium(url)
                except Exception as e:
                    logger.error("No success. The following error was raised:")
                    logger.error(e)

            # Create parser
            soup = BeautifulSoup(c, 'html.parser', from_encoding='utf-8')

            # Get stuff
            mydivs = soup.findAll("div", {"class": "gs_or"})

            for div in mydivs:
                try:
                    links.append(div.find('h3').find('a').get('href'))
                except:  # catch *all* exceptions
                    links.append('Look manually at: ' + url)

                try:
                    title.append(div.find('h3').find('a').text)
                except:
                    title.append('Could not catch title')

                try:
                    citations.append(get_citations(str(div.format_string)))
                except:
                    logger.warning("Number of citations not found for {}. Appending 0".format(title[-1]))
                    citations.append(0)

                try:
                    year.append(get_year(div.find('div', {'class': 'gs_a'}).text))
                except:
                    logger.warning("Year not found for {}, appending 0".format(title[-1]))
                    year.append(0)

                try:
                    author.append(get_author(div.find('div', {'class': 'gs_a'}).text))
                except:
                    author.append("Author not found")

                try:
                    publisher.append(div.find('div', {'class': 'gs_a'}).text.split("-")[-1])
                except:
                    publisher.append("Publisher not found")

                try:
                    venue.append(" ".join(div.find('div', {'class': 'gs_a'}).text.split("-")[-2].split(",")[:-1]))
                except:
                    venue.append("Venue not fount")

                rank.append(rank[-1] + 1)

            # Delay
            sleep(0.5)

        # Create a dataset and sort by the number of citations
        data = pd.DataFrame(list(zip(author, title, citations, year, publisher, venue, links)), index=rank[1:],
                            columns=['Author', 'Title', 'Citations', 'Year', 'Publisher', 'Venue', 'Source'])
        data.index.name = 'Rank'

        # Add columns with number of citations per year
        data['cit/year'] = data['Citations'] / (end_year + 1 - data['Year'])
        data['cit/year'] = data['cit/year'].round(0).astype(int)

        # Sort by the selected columns, if exists
        try:
            data_ranked = data.sort_values(by=sort_by, ascending=False)
        except Exception as e:
            print('Column name to be sorted not found. Sorting by the number of citations...')
            data_ranked = data.sort_values(by='Citations', ascending=False)
            print(e)
        # fix index
        data_ranked.reset_index(drop=True, inplace=True)
        # Print data
        print(data_ranked)

        # Plot by citation number
        if plot:
            plt.plot(rank[1:], citations, '*')
            plt.ylabel('Number of Citations')
            plt.xlabel('Rank of the keyword on Google Scholar')
            plt.title('Keyword: ' + keyword)
            plt.show()

        # Save results
        if save_table:
            fpath_csv = os.path.join(csv_path, keyword.replace(' ', '_') + '.csv')
            fpath_csv = fpath_csv[:MAX_CSV_FNAME]
            data_ranked.to_csv(fpath_csv, encoding='utf-8')
