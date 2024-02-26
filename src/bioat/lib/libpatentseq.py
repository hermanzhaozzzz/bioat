import time
import re
import sys
import json
import pandas as pd
from bioat import get_logger
from playwright.sync_api import Playwright, sync_playwright
from bioat.lib.libspider import get_random_user_agents, ProxyPool

__module_name__ = 'bioat.lib.libpatentseq'


STATUS = 'RUN'
SEQ_HEADER = None
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
    # headers = {
    #     'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)'
    # }
    headers = {
        'User-Agent': get_random_user_agents()
    }
    logger.debug(f'set headers: {headers}')
    url_login = "https://www.lens.org/lens/bio/patseqfinder"
    url_query = 'https://www.lens.org/lens/bio/psf/api/searchformdata'
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
    browser = playwright.firefox.launch(headless=headless)
    proxy_server = {"server": proxy_server} if proxy_server else None
    if proxy_server:
        logger.info(f'set proxy_server: {proxy_server}')

    context = browser.new_context(
        viewport={"width": 1920, "height": 1080} if not headless else None,
        proxy=proxy_server
    )
    page = context.new_page()
    page.set_extra_http_headers(headers)
    logger.debug(f'goto {url_login}')
    page.goto(url_login)
    page.get_by_text("Thanks, Got It").click()
    page.get_by_role("link", name="Sign in", exact=True).click()
    logger.debug(f'try to fill account info: {account}')
    page.get_by_label("Email address or Username").click()
    page.get_by_label("Email address or Username").fill(account['username'])
    page.get_by_label("Email address or Username").press("Tab")
    page.get_by_label("Password").fill(account['password'])
    page.get_by_label("Password").press("Tab")
    page.frame_locator("iframe[title=\"Widget containing a Cloudflare security challenge\"]").get_by_label(
        "Verify you are human").check()
    # must sleep for passing cloudflare
    time.sleep(2)
    logger.debug(f'try to sign in account')
    page.get_by_role("button", name="Sign in").click()

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
