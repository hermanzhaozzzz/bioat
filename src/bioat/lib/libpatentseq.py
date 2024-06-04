"""_summary_

author: Herman Huanan Zhao
email: hermanzhaozzzz@gmail.com
homepage: https://github.com/hermanzhaozzzz

_description_

example 1:
    bioat list
        <in shell>:
            $ bioat list
        <in python consolo>:
            >>> from bioat.cli import Cli
            >>> bioat = Cli()
            >>> bioat.list()
            >>> print(bioat.list())

example 2:
    _example_
"""

import json
import os
import re
import sys
import time

import pandas as pd

from bioat.lib.libpandas import set_option
from bioat.lib.libpath import HOME
from bioat.logger import get_logger

__module_name__ = "bioat.lib.libpatentseq"

try:
    from playwright._impl._errors import TargetClosedError, TimeoutError
    from playwright.sync_api import Playwright, sync_playwright
except ImportError:
    logger = get_logger(__module_name__)
    logger.error(
        "Unable to import playwright. please exec `python -m pip install playwright`, then try again."
    )
    sys.exit(1)

STATUS = "RUN"
SEQ_HEADER = None
COOKIE = f"{HOME}/.bioat/LENS.ORG/cookie.json"
CSS_PATH = os.path.join(os.path.dirname(__file__), "patentseq")
IP_ERRORS = (
    "NS_ERROR_UNKNOWN_HOST",
    "NS_ERROR_NET_INTERRUPT",
    "NS_ERROR_PROXY_FORBIDDEN",
)


def save_cookies(context, log_level):
    logger = get_logger(
        level=log_level, module_name=__module_name__, func_name="save_cookies"
    )
    # 保存 cookies
    storage_state = context.storage_state()
    logger.info(f"Saving cookie to {COOKIE}")
    os.makedirs(os.path.dirname(COOKIE), exist_ok=True)
    with open(COOKIE, "wt") as f:
        f.write(json.dumps(storage_state))
    logger.info(f"Saving cookie time to {COOKIE}.time")
    with open(f"{COOKIE}.time", "wt") as f:
        f.write(str(int(time.time())))


def load_cookies(browser, proxy_ip=None, log_level="DEBUG") -> None | object:
    # 加载 cookies
    logger = get_logger(
        level=log_level, module_name=__module_name__, func_name="load_cookies"
    )
    try:
        with open(f"{COOKIE}.time", "rt") as f:
            login_time = int(f.read().strip())
    except FileNotFoundError:
        login_time = None

    if login_time is None:
        # run login
        return None
    else:
        # run check time
        current_time = int(time.time())
        elapsed_time = current_time - login_time
        if elapsed_time > 60 * 60:  # 1 hour
            logger.info("Cookies are expired, performing re-login...")
            return None
        else:
            if (
                os.path.exists(COOKIE)
                and os.path.isfile(COOKIE)
                and os.path.getsize(f"{COOKIE}") > 0
            ):
                logger.info(
                    "Cookies are still valid, skip login and load cookies"
                )
                with open(COOKIE, "rt") as f:
                    storage_state = json.loads(f.read().strip())
                logger.debug(f"new_context")
                context = browser.new_context(
                    storage_state=storage_state,
                    proxy={"server": proxy_ip} if proxy_ip else None,
                )
                return context
            else:
                return None


def remove_cookie(log_level):
    logger = get_logger(
        level=log_level, module_name=__module_name__, func_name="remove_cookie"
    )
    logger.info("Removing cookie")
    try:
        os.remove(COOKIE)
        os.remove(f"{COOKIE}.time")
    except FileNotFoundError:
        logger.warning("Cookies are not found, skip remove.")


def run(
    playwright: Playwright,
    username,
    password,
    seq,
    seq_header,
    proxy_server,
    output,
    headless,
    nretry,
    local_browser,
    rm_fail_cookie,
    log_level,
) -> None:
    logger = get_logger(
        level=log_level, module_name=__module_name__, func_name="run"
    )

    account = {
        "username": username,
        "password": password,
    }
    logger.debug(f"Set account: {account}")
    # headers = {  # not change it, for bypassing cloudflare challenge with firefox
    #     'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)'
    # }

    # logger.debug(f'Set headers: {headers}')
    url_login = "https://www.lens.org/lens/bio/"
    url_query = "https://www.lens.org/lens/bio/patseqfinder"
    logger.debug(f"set url_login: {url_login}")
    logger.debug(f"set url_query: {url_query}")
    data = {
        "queryString": seq,
        "database": "patseq-aa",
        "sequenceType": "AA",
        "blastTool": "blastp",
        "numAlignments": "50",
        "expValue": "0.1",
        "matrix": "BLOSUM62",
        "seqLength": f"{len(seq)}",
    }
    logger.debug(f"Set data: {data}\n")
    logger.debug(f"Set headless: {headless}")

    proxy_ip = None
    proxy_pool = None

    if proxy_server:
        proxy_pool = ProxyPool(url=proxy_server)
        while True:
            proxy_ip = proxy_pool.get_proxy().get("proxy")
            logger.info(f"Set proxy_server: {proxy_server}")
            if proxy_ip:
                logger.debug("Get valid proxy, go on")
                break
            else:
                logger.debug(
                    "No proxy found in proxy pool, waiting 10 seconds..."
                )
                time.sleep(10)
                logger.debug("Next try to get proxy")
                continue
        logger.info(f"Get proxy_ip: {proxy_ip}")

    # start to test
    browser = playwright.firefox.launch(
        headless=headless,
        executable_path=local_browser if local_browser else None,
    )
    logger.debug("Try to load cookies")

    context = load_cookies(browser, proxy_ip, log_level)

    # context = None
    if context is None:  # need login
        logger.info("Cookies not exist or out of time. try to login")

        # login part
        counter = 0

        while True:
            if counter > nretry:
                logger.error(
                    f"Could not load url {url_login}, already retry maximum ({nretry})number of retries"
                )
                logger.debug("browser close")
                browser.close()
                logger.error("Login failed")
                sys.exit(1)
            try:
                logger.debug(f"new_context")
                context = browser.new_context(
                    viewport=(
                        {"width": 1920, "height": 1080}
                        if not headless
                        else None
                    ),
                    proxy={"server": proxy_ip} if proxy_ip else None,
                )
                page = context.new_page()
                # page.set_extra_http_headers(headers)
                logger.debug(f"Goto {url_login}")
                page.goto(url_login)
                logger.debug(f'Click "Thanks, Got It')
                page.get_by_text("Thanks, Got It").click()
                logger.debug(f'Click "Sign in')
                page.get_by_role("link", name="Sign in").click()
                # logger.debug(f'Try to fill account info: {account}')
                # page.get_by_label("Email address or Username").click()
                # page.get_by_label("Email address or Username").fill(account['username'])
                # page.get_by_label("Email address or Username").press("Tab")
                # page.get_by_label("Password").fill(account['password'])
                # add remember me
                # page.locator("label").filter(has_text="Keep me logged in").locator("i").click()
                # bypass cloudflare challenge
                # page.frame_locator("iframe[title=\"Widget containing a Cloudflare security challenge\"]").get_by_label(
                #     "Verify you are human").check()

                # must sleep for passing cloudflare
                # time.sleep(2)
                # logger.debug(f'Try to sign in account')
                # page.get_by_role("button", name="Sign in").click()
                with page.expect_popup() as page1_info:
                    logger.debug(f'Click "SignIn with ORCID"')
                    page.locator("a").filter(
                        has_text="SignIn with ORCID"
                    ).click()
                page1 = page1_info.value
                logger.debug(f'Click "Proceed"')
                page1.get_by_role("link", name="Proceed").click()
                time.sleep(2)
                logger.debug('Click "Accept All Cookies"')
                page1.get_by_role("button", name="Accept All Cookies").click()
                time.sleep(2)
                logger.debug('Click "Email or 16-digit ORCID iD"')
                page1.get_by_label("Email or 16-digit ORCID iD").click()
                logger.debug('Fill "Email or 16-digit ORCID iD"')
                page1.get_by_label("Email or 16-digit ORCID iD").fill(
                    account["username"]
                )
                logger.debug("Press Tab")
                page1.get_by_label("Email or 16-digit ORCID iD").press("Tab")
                logger.debug('Fill "Password"')
                page1.get_by_label("Password").fill(account["password"])
                time.sleep(2)
                logger.debug('Click "SIGN IN"')
                page1.get_by_role("button", name="SIGN IN", exact=True).click()
                logger.debug("Waiting...")
                time.sleep(20)
                page1.close()

                save_cookies(context, log_level)
                break  # to else for break
            except TimeoutError as e:
                logger.debug("context close")
                context.close()
                logger.warning(f"Could not load url {url_login}: {e}, retry...")
                time.sleep(3)
                counter += 1
                continue
            except Exception as e:
                if str(e) in IP_ERRORS:
                    logger.error(
                        f"proxy error: {e}, delete ip ({proxy_ip}) from remote server and retry..."
                    )
                    proxy_pool.delete_proxy(proxy_ip)
                    proxy_ip = proxy_pool.get_proxy().get("proxy")
                    logger.debug(f"Get new proxy_ip: {proxy_ip} and retry...")

                logger.debug("context close")
                context.close()
                logger.warning(f"Could not load url {url_login}: {e}, retry...")
                time.sleep(3)
                counter += 1
                continue
    # now login status is ok
    page = context.new_page()
    # page.set_extra_http_headers(headers)
    logger.debug(f"Goto {url_query}")
    # query
    counter = 0
    while True:
        if counter > nretry:
            logger.error(
                f"Could not load url {url_query}, already retry maximum ({nretry})number of retries"
            )
            logger.debug("context close")
            context.close()
            logger.debug("browser close")
            browser.close()
            logger.error("Query failed")
            sys.exit(1)
        try:
            logger.debug("wait until networkidle")
            page.goto(url_query, wait_until="networkidle")
            checker1 = page.locator("a.parent", has_text="Signed in as").count()
            checker2 = (
                page.frame_locator("iframe")
                .get_by_placeholder("Enter a query sequence.")
                .count()
            )

            if checker1:
                if checker2:
                    logger.debug("checker1 and checker2 pass!")
                    login_status = True
                    rm_cookie = False
                else:
                    login_status = False
                    rm_cookie = False
                    logger.error(
                        "checker2 failed, it might because that the proxy IP has been banned!"
                    )
            else:
                login_status = False
                rm_cookie = True
                logger.error(
                    "checker1 failed, it might because that the cookies out of time!"
                )

            if not login_status:
                logger.error(
                    f"Login failed! checker1: {bool(checker1)}, checker2: {bool(checker2)}"
                )
                logger.debug("context close")
                context.close()
                logger.debug("browser close")
                browser.close()
                logger.error(
                    f"Query failed! checker1: {bool(checker1)}, checker2: {bool(checker2)}"
                )
                if rm_cookie:
                    logger.debug(
                        "Consider to remove cookies because checker1 is failed!"
                    )
                    if rm_fail_cookie:
                        logger.debug(f"Param rm_fail_cookie = {rm_fail_cookie}")
                        remove_cookie(log_level)
                    else:
                        logger.debug(f"Param rm_fail_cookie = {rm_fail_cookie}")
                        logger.warning(f"User set to not remove cookies.")
                sys.exit(1)
            else:
                logger.debug("Passed cookies and success login!")
            time.sleep(3)
            page.frame_locator("iframe").get_by_placeholder(
                "Enter a query sequence."
            ).click()
            logger.debug(f"Try to fill seq: {seq}")
            page.frame_locator("iframe").get_by_placeholder(
                "Enter a query sequence."
            ).fill(seq)
            logger.debug(f'Click "Protein"')
            page.frame_locator("iframe").get_by_text(
                "Protein", exact=True
            ).click()
            logger.debug(f'Click "advanced options"')
            page.frame_locator("iframe").locator("div").filter(
                has_text=re.compile(r"^advanced options$")
            ).locator("div").click()
            logger.debug(f"Set maximum number of hits to 50")
            page.frame_locator("iframe").get_by_label(
                "Maximum Number of Hits to"
            ).select_option("50")
            logger.debug(f"Set maximum e-value to 1.0")
            page.frame_locator("iframe").get_by_label(
                "Expectation value threshold"
            ).select_option("1.0")
            logger.debug("Submit blastp query info")
            page.frame_locator("iframe").get_by_role(
                "button", name="Submit search"
            ).nth(1).click()
            time.sleep(3)
            logger.info("Task has been submitted, waiting for completion")
            break
        except TimeoutError as e:
            logger.warning(f"Could not load url {url_query}: {e}, retry...")
            time.sleep(3)
            counter += 1
            continue
        except Exception as e:
            if str(e) in IP_ERRORS:
                logger.error(
                    f"proxy error: {e}, delete ip ({proxy_ip}) from remote server and retry..."
                )
                proxy_pool.delete_proxy(proxy_ip)
                proxy_ip = proxy_pool.get_proxy().get("proxy")
                logger.debug(f"Get new proxy_ip: {proxy_ip} and retry...")
            logger.warning(f"Could not load url {url_query}: {e}, retry...")
            time.sleep(3)
            counter += 1
            continue

    # set global var to check program status
    global STATUS
    global SEQ_HEADER
    STATUS = "RUN"
    SEQ_HEADER = seq_header

    def handle_response(response):
        global STATUS
        global SEQ_HEADER

        if (
            "#/results" in full_url
            and response.request.method == "POST"
            and response.status == 200
        ):
            try:
                # 获取响应文本
                response_text = response.text()
                # 尝试解析 JSON，但先确保它不是空的
                if response_text:
                    json_data = json.loads(response_text)
                    df = pd.DataFrame.from_dict(json_data["hits"])
                    set_option(log_level=log_level)
                    logger.info(f"Result hit:\n{df}")

                    if SEQ_HEADER:
                        SEQ_HEADER = SEQ_HEADER.replace(">", "").split(" ")[0]
                        df.insert(0, "queryName", SEQ_HEADER)
                        if output:
                            logger.info(f"Export hit info to {output}")
                            df.to_csv(output, index=False)
                    STATUS = "FINISHED"
                else:
                    logger.error("Response text is empty")
                    STATUS = "ResponseEmptyError"
            except json.JSONDecodeError as e:
                logger.error(f"Error decoding JSON: {e}")
                STATUS = "JSONDecodeError"
            except Exception as e:
                logger.error(f"Error while handling response: {e}")
                STATUS = "ERROR"
            finally:
                logger.info(f"STATUS: {STATUS}")
                logger.info("Task complete.")

    start_time = time.time()
    while True:  # more than 150 seconds, failed
        if STATUS != "RUN":
            logger.info(
                f"STATUS={STATUS}, SPEND={int(time.time() - start_time)}s"
            )
            break
        if int(time.time() - start_time) >= 300:
            logger.error(
                f"TimeoutError (waiting for more than 300s), "
                f"this might be a problem with browser was terminated external"
            )
            logger.error("You can try again later")
            exit(1)
        try:
            full_url = page.evaluate("window.location.href")
            # if time.time() - start_time % 5 == 0:
            #     logger.info(f'STATUS={STATUS}, SPEND={int(time.time() - start_time)}s, waiting for completion')
            #     logger.debug(f"Current URL: {page.url}")
            #     logger.debug(f"Full URL including hash: {full_url}")
            page.on("response", handle_response)
        except TargetClosedError as e:
            logger.error(
                f"Meet error{e}, this might be a problem with browser was terminated external"
            )
            logger.error("you can try again later")
            sys.exit(1)
    # ---------------------
    logger.debug("context close")
    context.close()
    logger.debug("browser close")
    browser.close()
    logger.info("End. Exit.")
    sys.exit(0)


def query_patent(
    seq: str,
    seq_header: str | None = None,
    username: str | None = None,
    password: str | None = None,
    proxy_server: str | None = None,
    output: str | None = None,
    headless: bool = True,
    nretry: int = 3,
    local_browser: str | None = None,
    rm_fail_cookie: bool = False,
    log_level: str = "INFO",
):
    with sync_playwright() as playwright:
        run(
            playwright,
            username=username,
            password=password,
            seq=seq,
            seq_header=seq_header,
            proxy_server=proxy_server,
            output=output,
            headless=headless,
            nretry=nretry,
            local_browser=local_browser,
            rm_fail_cookie=rm_fail_cookie,
            log_level=log_level,
        )


if __name__ == "__main__":
    query_patent(
        seq="MEDDKKTTDSIRYELKDKHFWAAFLNLARHNVYITVNHINKILEEGEINRDGY",
        username="your_account",
        password="abc123456",
        proxy_server="http://127.0.0.1:8234",
    )
