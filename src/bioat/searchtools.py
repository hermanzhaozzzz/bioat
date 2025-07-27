"""Some cmdline tools for searching, powered by playwright."""

import importlib
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from bioat.exceptions import BioatInvalidOptionError
from bioat.lib.libpatentseq import query_patent
from bioat.logger import LoggerManager

lm = LoggerManager(mod_name="bioat.searchtools")


class SearchTools:
    """Search toolbox."""

    lm.set_names(cls_name="SearchTools")

    def __init__(self):
        pass

    def google_scholar(
        self,
        keyword: str,
        output: str | Path | None = None,
        sort_by: str = "CitePerYear",
        n_results: int = 100,
        plot: bool = False,
        start_year: int | None = None,
        end_year: int = datetime.now().year,
        log_level: str = "WARNING",
        **kwargs,
    ):
        """Search Google Scholar.

        This method creates a DataFrame of publication data from Google Scholar.
        Each result includes title, citations, year, authors, venue, publisher, and link.
        It is useful for finding relevant papers by citation metrics.

        Optionally, it can generate a plot of citation counts versus rank and save the
        table in various formats.

        Args:
            keyword (str): Keyword to search for. For exact matches, wrap in double and single quotes,
                e.g., "'exact keyword'".
            output (str, optional): Output file path. Supported formats: .csv, .tsv, .xls, .xlsx.
                Default is None, which means no output and only print the table in the console.
            sort_by (str): Column to sort the result by, such as "Citations" or "CitePerYear". Default is "CitePerYear".
            n_results (int): Number of search results to retrieve. Default is 100.

            plot (bool): Whether to plot citation count vs. rank. Default is False.
            start_year (int, optional): Optional start year for publication filtering.
            end_year (int): End year for publication filtering. Default is current year.
            log_level (str): Logging level. One of: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'NOTSET'.

        Returns:
            None. Prints the table and optionally saves or plots it.
        """
        lm.set_names(func_name="google_scholar")
        lm.set_level(log_level)
        lm.logger.info("Run google scholar search...")
        # show params
        lm.logger.info(f"Keyword: {keyword}")
        lm.logger.info(f"Sort by: {sort_by}")
        lm.logger.info(f"Number of results: {n_results}")
        lm.logger.info(f"Start year: {start_year}")
        lm.logger.info(f"End year: {end_year}")
        lm.logger.info(f"Output: {output}")
        lm.logger.info(f"Plot: {plot}")

        try:
            playwright_module = importlib.import_module("playwright.sync_api")
            sync_playwright = playwright_module.sync_playwright
        except ModuleNotFoundError:
            lm.logger.exception(
                "Missing dependency 'playwright'. Please install it by running:\n"
                "    pip install playwright && playwright install"
            )
            sys.exit(1)

        GSCHOLAR_URL = (
            "https://scholar.google.com/scholar?start={}&q={}&hl=en&as_sdt=0,5"
        )
        STARTYEAR_URL = "&as_ylo={}"
        ENDYEAR_URL = "&as_yhi={}"

        def get_year(text: str) -> int:
            for token in text.split():
                if token.strip().isdigit() and 1900 < int(token) <= datetime.now().year:
                    return int(token)
            return 0

        def get_citations(text: str) -> int:
            """Extract citation number from a string like 'Cited by 123'."""
            if "Cited by" in text:
                try:
                    return int(text.split("Cited by")[1].split()[0])
                except Exception:
                    return 0
            return 0

        def get_author(text: str) -> str:
            return text.split("-")[0].strip()

        def get_venue(text: str) -> str:
            parts = text.split("-")
            return parts[-2].strip() if len(parts) >= 2 else ""

        def get_publisher(text: str) -> str:
            parts = text.split("-")
            return parts[-1].strip() if len(parts) >= 2 else ""

        def launch_browser(p, headless=True):
            return p.chromium.launch(headless=headless)

        GSCHOLAR_MAIN_URL = GSCHOLAR_URL
        if start_year:
            GSCHOLAR_MAIN_URL += STARTYEAR_URL.format(start_year)
        if end_year:
            GSCHOLAR_MAIN_URL += ENDYEAR_URL.format(end_year)

        with sync_playwright() as p:
            browser = launch_browser(p, headless=kwargs.get("headless", True))
            context = browser.new_context(
                user_agent=(
                    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) "
                    "Chrome/120.0.0.0 Safari/537.36"
                ),
                viewport={"width": 1920, "height": 1080},
                locale="en-US",
                java_script_enabled=True,
            )
            page = context.new_page()

            links, title, citations, year, author, venue, publisher, rank = (
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                [0],
            )

            for n in range(0, n_results, 10):
                url = GSCHOLAR_MAIN_URL.format(str(n), keyword.replace(" ", "+"))
                page.goto(url, timeout=60000)
                page.wait_for_timeout(1000)
                # 检查是否触发验证码：页面空、或包含“unusual traffic”、“not a robot”等文字
                if (
                    not page.locator(".gs_r.gs_or.gs_scl").count()
                    or "unusual traffic" in page.content()
                ):
                    lm.logger.warning(
                        "\nGoogle Scholar 验证触发, 进入人工处理模式...\n"
                    )
                    browser.close()

                    # 重开非 headless 模式
                    browser = launch_browser(p, headless=False)
                    context = browser.new_context(
                        user_agent=(
                            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                            "AppleWebKit/537.36 (KHTML, like Gecko) "
                            "Chrome/120.0.0.0 Safari/537.36"
                        ),
                        viewport={"width": 1920, "height": 1080},
                        locale="en-US",
                        java_script_enabled=True,
                    )
                    page = context.new_page()
                    page.goto(url)
                    input("请在弹出的浏览器中手动完成验证, 完成后按回车继续...")

                items = page.locator(".gs_r.gs_or.gs_scl").all()
                for it in items:
                    try:
                        link = (
                            it.locator("h3 > a").get_attribute("href") or "Unavailable"
                        )
                        ttl = it.locator("h3 > a").inner_text()
                    except Exception:
                        link = "Unavailable"
                        ttl = "No Title"
                    meta = it.locator(".gs_a").inner_text()
                    title.append(ttl)
                    links.append(link)
                    # citations.append(get_citations(it.inner_html()))
                    try:
                        cite_texts = it.locator(".gs_fl").all_inner_texts()
                        citation_line = next(
                            (t for t in cite_texts if "Cited by" in t), ""
                        )
                        citations.append(get_citations(citation_line))
                    except Exception:
                        citations.append(0)
                    year.append(get_year(meta))
                    author.append(get_author(meta))
                    venue.append(get_venue(meta))
                    publisher.append(get_publisher(meta))
                    rank.append(rank[-1] + 1)

            browser.close()

        df = pd.DataFrame(
            zip(author, title, citations, year, publisher, venue, links, strict=False),
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
        df["CitePerYear"] = df["Citations"] / (
            end_year + 1 - df["Year"].replace(0, end_year + 1)
        )
        df["CitePerYear"] = df["CitePerYear"].round().astype(int)

        df = df[
            [
                "Author",
                "Citations",
                "CitePerYear",
                "Year",
                "Venue",
                "Title",
                "Publisher",
                "Source",
            ]
        ]

        try:
            df_sorted = df.sort_values(by=sort_by, ascending=False).reset_index(
                drop=True
            )
        except Exception:
            lm.logger.warning("Sort failed, using default sort by 'Citations'")
            df_sorted = df.sort_values(by="Citations", ascending=False).reset_index(
                drop=True
            )

        if plot:
            try:
                # 保证 citation 来自最终排序结果（而不是原始顺序）
                citation_values = df_sorted["Citations"].tolist()

                plt.figure(figsize=(8, 5))
                plt.plot(
                    range(1, len(citation_values) + 1),
                    citation_values,
                    marker="o",
                    linestyle="None",
                )
                plt.yscale("log")  # 常用: log scale 能更清晰看分布
                plt.grid(True, which="both", linestyle="--", linewidth=0.5)

                plt.ylabel("Number of Citations (log scale)")
                plt.xlabel("Rank")
                plt.title(f"Keyword: {keyword} ({len(citation_values)} results)")
                plt.tight_layout()
                plt.show()
            except Exception as e:
                lm.logger.warning(f"Failed to plot citation distribution: {e}")

        if output:
            lm.logger.info(f"Exporting results to {output}")
            if not output:
                output = f"{keyword.replace(' ', '_')}.xlsx"
            if output.endswith(".csv"):
                df_sorted.to_csv(output, index=False)
            elif output.endswith(".tsv"):
                df_sorted.to_csv(output, sep="\t", index=False)
            elif output.endswith(".xlsx") or output.endswith(".xls"):
                df_sorted.to_excel(output, index=False)
            else:
                raise BioatInvalidOptionError("Unsupported output format: " + output)
        else:
            print(df_sorted)
        if kwargs.get("_for_test", False):
            return df_sorted
        return None

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
        """Return a table with a list of patent blast hit from lens.org.

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
        lm.set_names(func_name="query_patent")
        lm.set_level(log_level)

        lm.logger.info("Run query patent sequence from lens.org...")
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


if __name__ == "__main__":
    tool = SearchTools()
    tool.google_scholar(keyword="prime editing", n_results=20, plot=False)
