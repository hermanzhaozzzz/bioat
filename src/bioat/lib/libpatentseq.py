import time
import re
import os
import sys
import json
import pandas as pd
from selenium.common import NoSuchElementException

from bioat import get_logger
from bioat.lib.libpath import HOME
from bioat.lib.libpandas import set_option
from playwright.sync_api import Playwright, sync_playwright

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
        if elapsed_time > 4 * 60 * 60:  # 4 hours
            logger.info('Cookies are expired, performing re-login...')
            return None
        else:
            if os.path.exists(COOKIE) and os.path.isfile(COOKIE) and os.path.getsize(f'{COOKIE}') > 0:
                logger.info('Cookies are still valid, skip login and load cookies')
                with open(COOKIE, 'rt') as f:
                    storage_state = json.loads(f.read().strip())
                context = browser.new_context(storage_state=storage_state)
                return context
            else:
                return None


def run(
        playwright: Playwright,
        username,
        password,
        seq,
        seq_header,
        proxy_server,
        output,
        headless,
        log_level
) -> None:
    logger = get_logger(level=log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

    account = {
        'username': username,
        'password': password,
    }
    logger.debug(f'set account: {account}')
    headers = {  # not change it, for bypassing cloudflare challenge with firefox
        'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)'
    }

    logger.debug(f'set headers: {headers}')
    url_login = "https://www.lens.org/lens/bio"
    # url_query = 'https://www.lens.org/lens/bio/psf/api/searchformdata'
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
    logger.debug(f'set data: {data}\n')
    logger.debug(f'set headless: {headless}')
    proxy_server = {"server": proxy_server} if proxy_server else None
    if proxy_server:
        logger.info(f'set proxy_server: {proxy_server}')

    browser = playwright.firefox.launch(headless=headless)

    logger.debug('try to load cookies')
    context = load_cookies(browser, log_level)
    if context is None:
        # need login
        logger.info('cookies not exist or out of time. try to login')
        context = browser.new_context(
            viewport={"width": 1920, "height": 1080} if not headless else None,
            proxy=proxy_server
        )
        page = context.new_page()
        page.set_extra_http_headers(headers)
        logger.debug(f'goto {url_login}')
        counter = 0
        while counter <= 3:
            try:
                page.goto(url_login)
                break
            except Exception as e:
                logger.warning(f'Could not load url {url_login}: {e}, retry...')
                time.sleep(3)
                counter += 1
                continue
        page.get_by_text("Thanks, Got It").click()
        page.get_by_role("link", name="Sign in", exact=True).click()
        logger.debug(f'try to fill account info: {account}')
        page.get_by_label("Email address or Username").click()
        page.get_by_label("Email address or Username").fill(account['username'])
        page.get_by_label("Email address or Username").press("Tab")
        page.get_by_label("Password").fill(account['password'])
        # add remember me
        page.get_by_role("group").locator("label").filter(has_text="Keep me logged in").locator("i").click()
        # bypass cloudflare challenge
        page.frame_locator("iframe[title=\"Widget containing a Cloudflare security challenge\"]").get_by_label(
            "Verify you are human").check()
        # must sleep for passing cloudflare
        time.sleep(2)
        logger.debug(f'try to sign in account')
        page.get_by_role("button", name="Sign in").click()
        save_cookies(context, log_level)
    # now login status is ok
    page = context.new_page()
    page.set_extra_http_headers(headers)
    logger.debug(f'goto {url_query}')
    counter = 0
    while counter <= 3:
        try:
            page.goto(url_query)
            break
        except Exception as e:
            logger.warning(f'Could not load url {url_query}: {e}, retry...')
            time.sleep(3)
            counter += 1
            continue
    page.frame_locator("iframe").get_by_placeholder("Enter a query sequence.").click()
    time.sleep(3)
    logger.debug(f'try to fill seq: {seq}')
    page.frame_locator("iframe").get_by_placeholder("Enter a query sequence.").fill(seq)
    page.frame_locator("iframe").get_by_text("Protein", exact=True).click()
    page.frame_locator("iframe").locator("div").filter(has_text=re.compile(r"^advanced options$")).locator(
        "div").click()
    logger.debug(f'set maximum number of hits to 50')
    page.frame_locator("iframe").get_by_label("Maximum Number of Hits to").select_option("50")
    logger.debug(f'set maximum e-value to 1.0')
    page.frame_locator("iframe").get_by_label("Expectation value threshold").select_option("1.0")
    logger.debug('submit blastp query info')
    page.frame_locator("iframe").get_by_role("button", name="Submit search").nth(1).click()
    logger.info('task has been submitted, waiting for completion')

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
                    set_option()
                    logger.info(f'result hit:\n{df}')
                    if SEQ_HEADER:
                        SEQ_HEADER = SEQ_HEADER.replace('>', '').split(' ')[0]
                        df.insert(0, 'queryName', SEQ_HEADER)
                        logger.info(f'export hit info to {output}')
                        df.to_csv(output, index=False)
                    STATUS = 'FINISHED'
                else:
                    logger.error("Response text is empty")
                    STATUS = 'ERROR'
            except json.JSONDecodeError as e:
                logger.error(f"Error decoding JSON: {e}")
                STATUS = 'ERROR'
            except Exception as e:
                logger.error(f"Error while handling response: {e}")
                STATUS = 'ERROR'
            finally:
                logger.info(f'STATUS: {STATUS}')
                logger.info('task complete.')

    counter = 0
    while STATUS == 'RUN' and counter * 2 <= 200:  # more than 200 seconds, failed
        time.sleep(2)
        logger.info(f'STATUS={STATUS}, waiting for completion')
        # print(f"Current URL: {page.url}")
        full_url = page.evaluate('window.location.href')
        # print(f"Full URL including hash: {full_url}")
        page.on('response', handle_response)
    # ---------------------
    context.close()
    browser.close()


def query_patent(
        seq: str,
        seq_header: str | None = None,
        username: str | None = None,
        password: str | None = None,
        proxy_server: str | None = None,  # proxy_server = "socks5://127.0.0.1:8235",
        output: str | None = None,
        headless: bool = True,
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
            log_level=log_level
        )


if __name__ == '__main__':
    query_patent(
        seq='MEDDKKTTDSIRYELKDKHFWAAFLNLARHNVYITVNHINKILEEGEINRDGYETTLKNTWNEIKDINKKDRLSKLIIKHFPFLEAATYRLNPTDTTKQKEEKQAEAQSLESLRKSFFVFIYKLRDLRNHYSHYKHSKSLERPKFEEGLLEKMYNIFNASIRLVKEDYQYNKDINPDEDFKHLDRTEEEFNYYFTKDNEGNITESGLLFFVSLFLEKKDAIWMQQKLRGFKDNRENKKKMTNEVFCRSRMLLPKLRLQSTQTQDWILLDMLNELIRCPKSLYERLREEDREKFRVPIEIADEDYDAEQEPFKNTLVRHQDRFPYFALRYFDYNEIFTNLRFQIDLGTYHFSIYKKQIGDYKESHHLTHKLYGFERIQEFTKQNRPDEWRKFVKTFNSFETSKEPYIPETTPHYHLENQKIGIRFRNDNDKIWPSLKTNSEKNEKSKYKLDKSFQAEAFLSVHELLPMMFYYLLLKTENTDNDNEIETKKKENKNDKQEKHKIEEIIENKITEIYALYDTFANGEIKSIDELEEYCKGKDIEIGHLPKQMIAILKDEHKVMATEAERKQEEMLVDVQKSLESLDNQINEEIENVERKNSSLKSGKIASWLVNDMMRFQPVQKDNEGKPLNNSKANSTEYQLLQRTLAFFGSEHERLAPYFKQTKLIESSNPHPFLKDTEWEKCNNILSFYRSYLEAKKNFLESLKPEDWEKNQYFLKLKEPKTKPKTLVQGWKNGFNLPRGIFTEPIRKWFMKHRENITVAELKRVGLVAKVIPLFFSEEYKDSVQPFYNYHFNVGNINKPDEKNFLNCEERRELLRKKKDEFKKMTDKEKEENPSYLEFKSWNKFERELRLVRNQDIVTWLLCMELFNKKKIKELNVEKIYLKNINTNTTKKEKNTEEKNGEEKNIKEKNNILNRIMPMRLPIKVYGRENFSKNKKKKIRRNTFFTVYIEEKGTKLLKQGNFKALERDRRLGGLFSFVKTPSKAESKSNTISKLRVEYELGEYQKARIEIIKDMLALEKTLIDKYNSLDTDNFNKMLTDWLELKGEPDKASFQNDVDLLIAVRNAFSHNQYPMRNRIAFANINPFSLSSANTSEEKGLGIANQLKDKTHKTIEKIIEIEKPIETKE',
        username='hermanzhaozzzz',
        password='zhn123,.',
        proxy_server="socks5://127.0.0.1:8235",
    )
