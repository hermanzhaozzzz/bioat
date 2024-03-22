import time
import re
import os
import sys
import json
import pandas as pd
from bioat import get_logger
from bioat.lib.libpath import HOME
from bioat.lib.libpandas import set_option
from playwright.sync_api import Playwright, sync_playwright
from playwright._impl._errors import TimeoutError, TargetClosedError

__module_name__ = 'bioat.lib.libpatentseq'

STATUS = 'RUN'
SEQ_HEADER = None
COOKIE = f'{HOME}/.bioat/LENS.ORG/cookie.json'


def save_cookies(context, log_level):
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
    # 保存 cookies
    storage_state = context.storage_state()
    logger.info(f'Saving cookie to {COOKIE}')
    os.makedirs(os.path.dirname(COOKIE), exist_ok=True)
    with open(COOKIE, 'wt') as f:
        f.write(json.dumps(storage_state))
    logger.info(f'Saving cookie time to {COOKIE}.time')
    with open(f'{COOKIE}.time', 'wt') as f:
        f.write(str(int(time.time())))


def load_cookies(browser, log_level) -> None | object:
    # 加载 cookies
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
    try:
        with open(f'{COOKIE}.time', 'rt') as f:
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
        if elapsed_time > 15 * 60:  # 1/4 hours 15min
            logger.info('Cookies are expired, performing re-login...')
            return None
        else:
            if os.path.exists(COOKIE) and os.path.isfile(COOKIE) and os.path.getsize(f'{COOKIE}') > 0:
                logger.info('Cookies are still valid, skip login and load cookies')
                with open(COOKIE, 'rt') as f:
                    storage_state = json.loads(f.read().strip())
                logger.debug(f'new_context')
                context = browser.new_context(storage_state=storage_state)
                return context
            else:
                return None


def remove_cookie(log_level):
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)
    logger.info('Removing cookie')
    try:
        os.remove(COOKIE)
        os.remove(f'{COOKIE}.time')
    except FileNotFoundError:
        logger.warning('Cookies are not found, skip remove.')


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
        log_level
) -> None:
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

    account = {
        'username': username,
        'password': password,
    }
    logger.debug(f'Set account: {account}')
    # headers = {  # not change it, for bypassing cloudflare challenge with firefox
    #     'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)'
    # }

    # logger.debug(f'Set headers: {headers}')
    url_login = "https://www.lens.org/lens/bio/"
    url_query = "https://www.lens.org/lens/bio/patseqfinder"
    logger.debug(f'set url_login: {url_login}')
    logger.debug(f'set url_query: {url_query}')
    data = {
        'queryString': seq,
        'database': 'patseq-aa',
        'sequenceType': 'AA',
        'blastTool': 'blastp',
        'numAlignments': '50',
        'expValue': '0.1',
        'matrix': 'BLOSUM62',
        'seqLength': f'{len(seq)}',
    }
    logger.debug(f'Set data: {data}\n')
    logger.debug(f'Set headless: {headless}')
    proxy_server = {"server": proxy_server} if proxy_server else None
    if proxy_server:
        logger.info(f'Set proxy_server: {proxy_server}')

    # start to test
    # browser = playwright.firefox.launch(headless=headless)
    browser = playwright.webkit.launch(headless=headless)
    logger.debug('Try to load cookies')

    context = load_cookies(browser, log_level)

    if context is None:  # need login

        logger.info('Cookies not exist or out of time. try to login')

        # login part
        counter = 0

        while True:
            if counter > nretry:
                logger.error(
                    f'Could not load url {url_login}, already retry maximum ({nretry})number of retries')
                remove_cookie(log_level)
                logger.debug('browser close')
                browser.close()
                logger.error('Login failed')
                sys.exit(1)
            try:
                logger.debug(f'new_context')
                context = browser.new_context(
                    viewport={"width": 1920, "height": 1080} if not headless else None,
                    proxy=proxy_server
                )
                page = context.new_page()
                # page.set_extra_http_headers(headers)
                logger.debug(f'Goto {url_login}')
                page.goto(url_login)
                page.get_by_text("Thanks, Got It").click()
                page.get_by_role("link", name="Sign in").click()
                logger.debug(f'Try to fill account info: {account}')
                page.get_by_label("Email address or Username").click()
                page.get_by_label("Email address or Username").fill(account['username'])
                page.get_by_label("Email address or Username").press("Tab")
                page.get_by_label("Password").fill(account['password'])
                # add remember me
                page.locator("label").filter(has_text="Keep me logged in").locator("i").click()
                # bypass cloudflare challenge
                # page.frame_locator("iframe[title=\"Widget containing a Cloudflare security challenge\"]").get_by_label(
                #     "Verify you are human").check()
                # must sleep for passing cloudflare
                time.sleep(2)
                logger.debug(f'Try to sign in account')
                page.get_by_role("button", name="Sign in").click()
                save_cookies(context, log_level)
                break  # to else for break
            except TimeoutError as e:
                logger.debug('context close')
                context.close()
                logger.warning(f'Could not load url {url_login}: {e}, retry...')
                time.sleep(3)
                counter += 1
                continue
            except Exception as e:
                logger.debug('context close')
                context.close()
                logger.warning(f'Could not load url {url_login}: {e}, retry...')
                time.sleep(3)
                counter += 1
                continue
    # now login status is ok
    page = context.new_page()
    # page.set_extra_http_headers(headers)
    logger.debug(f'Goto {url_query}')
    # query
    counter = 0
    while True:
        if counter > nretry:
            logger.error(
                f'Could not load url {url_query}, already retry maximum ({nretry})number of retries')
            remove_cookie(log_level)
            logger.debug('context close')
            context.close()
            logger.debug('browser close')
            browser.close()
            logger.error('Query failed')
            sys.exit(1)
        try:
            page.goto(url_query)
            page.frame_locator("iframe").get_by_placeholder("Enter a query sequence.").click()
            time.sleep(3)
            logger.debug(f'Try to fill seq: {seq}')
            page.frame_locator("iframe").get_by_placeholder("Enter a query sequence.").fill(seq)
            page.frame_locator("iframe").get_by_text("Protein", exact=True).click()
            page.frame_locator("iframe").locator("div").filter(has_text=re.compile(r"^advanced options$")).locator(
                "div").click()
            logger.debug(f'Set maximum number of hits to 50')
            page.frame_locator("iframe").get_by_label("Maximum Number of Hits to").select_option("50")
            logger.debug(f'Set maximum e-value to 1.0')
            page.frame_locator("iframe").get_by_label("Expectation value threshold").select_option("1.0")
            logger.debug('Submit blastp query info')
            page.frame_locator("iframe").get_by_role("button", name="Submit search").nth(1).click()
            logger.info('Task has been submitted, waiting for completion')
            break
        except TimeoutError as e:
            logger.warning(f'Could not load url {url_query}: {e}, retry...')
            time.sleep(3)
            counter += 1
            continue
        except Exception as e:
            logger.warning(f'Could not load url {url_query}: {e}, retry...')
            time.sleep(3)
            counter += 1
            continue

    # set global var to check program status
    global STATUS
    global SEQ_HEADER
    STATUS = 'RUN'
    SEQ_HEADER = seq_header

    def handle_response(response):
        global STATUS
        global SEQ_HEADER

        if (
                '#/results' in full_url and
                response.request.method == 'POST' and
                response.status == 200
        ):
            try:
                # 获取响应文本
                response_text = response.text()
                # 尝试解析 JSON，但先确保它不是空的
                if response_text:
                    json_data = json.loads(response_text)
                    df = pd.DataFrame.from_dict(json_data['hits'])
                    set_option(log_level=log_level)
                    logger.info(f'Result hit:\n{df}')

                    if SEQ_HEADER:
                        SEQ_HEADER = SEQ_HEADER.replace('>', '').split(' ')[0]
                        df.insert(0, 'queryName', SEQ_HEADER)
                        if output:
                            logger.info(f'Export hit info to {output}')
                            df.to_csv(output, index=False)
                    STATUS = 'FINISHED'
                else:
                    logger.error("Response text is empty")
                    STATUS = 'ResponseEmptyError'
            except json.JSONDecodeError as e:
                logger.error(f"Error decoding JSON: {e}")
                STATUS = 'JSONDecodeError'
            except Exception as e:
                logger.error(f"Error while handling response: {e}")
                STATUS = 'ERROR'
            finally:
                logger.info(f'STATUS: {STATUS}')
                logger.info('Task complete.')

    start_time = time.time()
    while True:  # more than 150 seconds, failed
        if STATUS != 'RUN':
            logger.info(f'STATUS={STATUS}, SPEND={int(time.time() - start_time)}s')
            break
        if int(time.time() - start_time) >= 150:
            logger.error(
                f"TimeoutError (waiting for more than 150s), "
                f"this might be a problem with browser was terminated external")
            logger.error('You can try again later')
            exit(1)
        try:
            full_url = page.evaluate('window.location.href')
            # if time.time() - start_time % 5 == 0:
            #     logger.info(f'STATUS={STATUS}, SPEND={int(time.time() - start_time)}s, waiting for completion')
            #     logger.debug(f"Current URL: {page.url}")
            #     logger.debug(f"Full URL including hash: {full_url}")
            page.on('response', handle_response)
        except TargetClosedError as e:
            logger.error(f"Meet error{e}, this might be a problem with browser was terminated external")
            logger.error('you can try again later')
            sys.exit(1)
    # ---------------------
    logger.debug('context close')
    context.close()
    logger.debug('browser close')
    browser.close()
    logger.info('End. Exit.')
    sys.exit(0)


def query_patent(
        seq: str,
        seq_header: str | None = None,
        username: str | None = None,
        password: str | None = None,
        proxy_server: str | None = None,  # proxy_server = "socks5://127.0.0.1:8235",
        output: str | None = None,
        headless: bool = True,
        nretry: int = 3,
        log_level: str = 'INFO'
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
            log_level=log_level
        )


if __name__ == '__main__':
    query_patent(
        seq='MEDDKKTTDSIRYELKDKHFWAAFLNLARHNVYITVNHINKILEEGEINRDGYETTLKNTWNEIKDINKKDRLSKLIIKHFPFLEAATYRLNPTDTTKQKEEKQAEAQSLESLRKSFFVFIYKLRDLRNHYSHYKHSKSLERPKFEEGLLEKMYNIFNASIRLVKEDYQYNKDINPDEDFKHLDRTEEEFNYYFTKDNEGNITESGLLFFVSLFLEKKDAIWMQQKLRGFKDNRENKKKMTNEVFCRSRMLLPKLRLQSTQTQDWILLDMLNELIRCPKSLYERLREEDREKFRVPIEIADEDYDAEQEPFKNTLVRHQDRFPYFALRYFDYNEIFTNLRFQIDLGTYHFSIYKKQIGDYKESHHLTHKLYGFERIQEFTKQNRPDEWRKFVKTFNSFETSKEPYIPETTPHYHLENQKIGIRFRNDNDKIWPSLKTNSEKNEKSKYKLDKSFQAEAFLSVHELLPMMFYYLLLKTENTDNDNEIETKKKENKNDKQEKHKIEEIIENKITEIYALYDTFANGEIKSIDELEEYCKGKDIEIGHLPKQMIAILKDEHKVMATEAERKQEEMLVDVQKSLESLDNQINEEIENVERKNSSLKSGKIASWLVNDMMRFQPVQKDNEGKPLNNSKANSTEYQLLQRTLAFFGSEHERLAPYFKQTKLIESSNPHPFLKDTEWEKCNNILSFYRSYLEAKKNFLESLKPEDWEKNQYFLKLKEPKTKPKTLVQGWKNGFNLPRGIFTEPIRKWFMKHRENITVAELKRVGLVAKVIPLFFSEEYKDSVQPFYNYHFNVGNINKPDEKNFLNCEERRELLRKKKDEFKKMTDKEKEENPSYLEFKSWNKFERELRLVRNQDIVTWLLCMELFNKKKIKELNVEKIYLKNINTNTTKKEKNTEEKNGEEKNIKEKNNILNRIMPMRLPIKVYGRENFSKNKKKKIRRNTFFTVYIEEKGTKLLKQGNFKALERDRRLGGLFSFVKTPSKAESKSNTISKLRVEYELGEYQKARIEIIKDMLALEKTLIDKYNSLDTDNFNKMLTDWLELKGEPDKASFQNDVDLLIAVRNAFSHNQYPMRNRIAFANINPFSLSSANTSEEKGLGIANQLKDKTHKTIEKIIEIEKPIETKE',
        username='hermanzhaozzzz',
        password='zhn123,.',
        proxy_server="socks5://127.0.0.1:8235",
    )
