import tarfile
import gzip
import time
import os
import sys
import textwrap
from collections import defaultdict
from hashlib import md5
from bioat.lib.libpath import HOME
from bioat.logger import get_logger

__module_name__ = 'bioat.lib.libjgi'


class JGIDoc:
    DEFAULT_CATEGORIES = ["ESTs",
                          "EST Clusters",
                          "Assembled scaffolds (unmasked)",
                          "Assembled scaffolds (masked)",
                          "Transcripts",
                          "Genes",
                          "CDS",
                          "Proteins",
                          "Additional Files"]
    usage_example_blurb = """\
        This script will retrieve files from JGI using the cURL api. It will
        return a list of possible files for downloading.
    
        * This script depends upon cURL - it can be downloaded here:
        http://curl.haxx.se/
    
        # USAGE ///////////////////////////////////////////////////////////////////////
        $ jgi-query.py [<jgi_address>, <jgi_abbreviation>] [[-xml [<your_xml>]], -f]
    
        To get <jgi_address>, go to: http://genome.jgi.doe.gov/ and search for your
        species of interest. Click through until you are at the "Info" page. For
        \x1B[3mNematostella vectensis\x1B[23m, the appropriate page is
        "http://genome.jgi.doe.gov/Nemve1/Nemve1.info.html".
    
        To query using only the name simply requires the specific JGI organism
        abbreviation, as referenced in the full url.
    
        For the above example, the proper input syntax for this script would be:
    
        $ jgi-query.py http://genome.jgi.doe.gov/Nemve1/Nemve1.info.html
    
                                 -or-
    
        $ jgi-query.py Nemve1
    
        If you already have the XML file for the query in the directory, you may use
        the --xml flag to avoid redownloading it (particularly useful if querying
        large, top-level groups with many sub-species, such as "fungi"):
    
        $ jgi-query.py --xml <your_xml_index>
    
        If the XML filename is omitted when using the --xml flag, it is assumed that
        the XML file is named '<jgi_abbreviation>_jgi_index.xml'. In such cases, the
        organism name is required.
        # /USAGE //////////////////////////////////////////////////////////////////////
        """
    select_blurb = """\
        # SYNTAX ///////////////////////////////////////////////////////////////////////
        Using the following format syntax to download selected file:
            <category number>:<i>[,<i>, <i>];<category number>:<i>-<i>;...
        
        Indices (<i>) may be a mixture of comma-separated values and hyphen-separated 
        ranges. 
        
        Example:
            '3:4,5; 7:1-10,13' will select elements 4 and 5 from category 3, and 1-10 
            plus 13 from category 7.
        # /SYNTAX ///////////////////////////////////////////////////////////////////////
        """


class JGIConfig:
    def __init__(self, overwrite_conf: bool = False):
        self.path = os.path.join(HOME, '.bioat', 'JGI', "account.conf")
        self.log_level = 'CRITICAL'
        self.info: dict = {"user": None, "password": None, "categories": None}

        # Check config exist or not
        if os.path.isfile(self.path) and not overwrite_conf:
            self.load()
        else:  # no config present or configure flag used; run config dialog
            self.get_user_info()
            self.save()
            self.load()

    def load(self):
        """
        Reads "user", "password" and "categories" entries
        from config file.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        logger.debug('loading JGI account info...')

        with open(self.path, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith("user"):
                    self.info["user"] = line.split("=")[1]
                if line.startswith("password"):
                    self.info["password"] = line.split("=")[1]
                if line.startswith("categories"):
                    cats = line.strip().split("=")[1]
                    self.info["categories"] = [e.strip() for e in cats.split(",")]

        if not all([self.info["user"], self.info["password"]]):
            logger.critical(f"Config file present ({self.path}), but user and/or password not found.")
            sys.exit(1)

        logger.debug('loading JGI account info done')

    def get_user_info(self):
        """
        Dialog with user to gather user information for
        use with the curl query. Returns a dict.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        blurb = """
        === USER SETUP ===

        JGI access configuration:

        Before continuing, you will need to provide your JGI login credentials.
        These are required by JGI's curl api, and will be stored in a config
        file for future use (unless you choose to delete them).

        If you need to sign up for a JGI account, use the registration link at
        https://contacts.jgi.doe.gov/registration/new

        === CREDENTIALS ===
        """
        # Print left-justified triple-quoted text blocks
        print(textwrap.dedent(blurb))
        user_query = "JGI account username/email (or 'q' to quit): "

        # write user and passwd from interactive terminal
        self.info["user"] = input(user_query)

        if self.info["user"] == "q":
            sys.exit("Exiting now.")

        self.info["password"] = input("JGI account password (or 'q' to quit): ")

        if self.info["password"] == "q":
            sys.exit("Exiting now.")

        input_blurb = (
            f"Proceed with USER={self.info['user']}, PASSWORD={self.info['password']} to configure "
            "script?\n([y]es, [n]o, [r]estart): "
        )

        self.info["categories"]: list[str] = JGIDoc.DEFAULT_CATEGORIES

        while True:  # catch invalid responses
            choice = input(input_blurb)

            if choice.lower() == "y":
                pass
            elif choice.lower() == "n":
                sys.exit("Exiting now.")
            elif choice.lower() == "r":
                self.get_user_info()

    def save(self):
        """
        Creates a config file <config_path> using
        credentials from dict <config_info>.

        """
        logger = get_logger(level=self.log_level, module_name=__module_name__, func_name=sys._getframe().f_code.co_name)

        logger.debug('saving JGI account info...')
        u = self.info["user"]
        p = self.info["password"]
        c = self.info["categories"]
        c = ",".join(c)
        header = ("# bioat meta JGI_query: JGI account info and essentials {}\n".format("#" * 24))
        info = f"user={u}\npassword={p}\ncategories={c}"

        with open(self.path, 'wt') as f:
            f.write(header)
            f.write(info)

        logger.debug('saving JGI account info done')

    def clean_exit(exit_message=None, exit_code=0, remove_temp=True):
        """
        Perform a sys.exit() while removing temporary files and
        informing the user.

        """
        to_remove = ["cookies"]
        # don't delete xml file if supplied by user
        if not LOCAL_XML and remove_temp is True:
            try:
                to_remove.append(xml_index_filename)
            except NameError:
                pass
        for f in to_remove:
            try:
                os.remove(f)
            except OSError:
                continue
        if remove_temp is True:
            base_message = "Removing temp files and exiting"
        else:
            base_message = "Keeping temp files and exiting"
        if exit_message:
            print(exit_message)
        print(base_message)

        sys.exit(exit_code)

    def get_user_choice():
        """
        Get user file selection choice(s)

        """
        choice = input(
            "Enter file selection ('q' to quit, "
            "'usage' to review syntax, 'a' for all, "
            "'r' for regex-based filename matching):\n> ")
        if choice == "usage":
            print()
            print(select_blurb)
            print()
            return get_user_choice()
        elif choice.lower() in ("q", "quit", "exit"):
            remove_temp = input("Remove index file? (y/n): ")
            remove_temp = remove_temp.lower() in ('y', 'yes', '')
            clean_exit(remove_temp=remove_temp)
        else:
            return choice

    def is_xml(filename):
        """
        Uses hex code at the beginning of a file to try to determine if it's an
        XML file or not. This seems to be occasionally necessary; if pulling
        files from JGI tape archives, the server may error out and provide an
        XML error document instead of the intended file. This function should
        return False on all downloaded files, although false positives have not
        been thoroughly investigated.

        Adapted from http://stackoverflow.com/a/13044946/3076552

        """
        xml_hex = "\x3c"  # hex code at beginning of XML files
        read_length = len(xml_hex)
        with open(filename) as f:
            try:
                file_start = f.read(read_length)
            except UnicodeDecodeError:  # compressed files
                return False
            if file_start.startswith(xml_hex):  # XML file
                return True
            else:  # hopefully all other file types
                return False

    def is_broken(filename, min_size_bytes=20, md5_hash=None, sizeInBytes=None):
        """
        Rudimentary check to see if a file appears to be broken.

        """
        if (
                not os.path.isfile(filename) or
                os.path.getsize(filename) < min_size_bytes or
                (is_xml(filename) and not filename.lower().endswith("xml")) or
                ((not check_md5(filename, md5_hash)) or
                 (not check_sizeInBytes(filename, sizeInBytes)))
        ):
            return True
        else:
            return False

    def check_md5(filename, md5_hash, print_message=True):
        if not md5_hash:
            message = "INFO: No MD5 hash listed for {}; skipping check".format(filename)
            ret_val = True
        else:
            file_md5 = get_md5(filename)
            if file_md5 == md5_hash:
                message = (
                    "SUCCESS: MD5 hashes match for {} ({})".format(filename, md5_hash))
                ret_val = True
            else:
                message = ("ERROR: MD5 hash mismatch for {} (local: {}, remote: {})"
                           .format(filename, file_md5, md5_hash))
                ret_val = False

        if print_message is True:
            print(message)

        return ret_val

    def check_sizeInBytes(filename, sizeInBytes, print_message=True):
        if not sizeInBytes:
            message = "INFO: No sizeInBytes listed for {}; skipping check".format(filename)
            ret_val = True
        else:
            file_sizeInBytes = get_sizeInBytes(filename)
            if file_sizeInBytes == sizeInBytes:
                message = (
                    "SUCCESS: sizeInBytes match for {} ({})".format(filename, sizeInBytes))
                ret_val = True
            else:
                message = ("ERROR: sizeInBytes mismatch for {} (local: {}, remote: {})"
                           .format(filename, file_sizeInBytes, sizeInBytes))
                ret_val = False

        if print_message is True:
            print(message)

        return ret_val

    def download_from_url(url, timeout=120, retry=0, min_file_bytes=20, url_to_validate={}):
        """
        Attempts to download a file from JGI servers using cURL.

        Returns a tuple of (filename, cURL command used, success boolean)

        """
        success = True
        md5_hash = url_to_validate[url].get("md5", None)
        sizeInBytes = url_to_validate[url].get("sizeInBytes", None)

        url = url.replace("&amp;", "&")

        filename = re.search('.+/(.+$)', url).group(1)
        url_prefix = "https://genome.jgi.doe.gov"
        download_command = (
            "curl -m {} '{}{}' -b cookies "
            "> {}".format(timeout, url_prefix, url, filename)
        )
        if not is_broken(filename, md5_hash=md5_hash, sizeInBytes=sizeInBytes):
            success = True
            print("Skipping existing file {}".format(filename))
        else:
            print("Downloading '{}' using command:\n{}"
                  .format(filename, download_command))
            # The next line doesn't appear to be needed to refresh the cookies.
            #    subprocess.call(login, shell=True)
            status = subprocess.run(download_command, shell=True).returncode
            if status != 0 or is_broken(
                    filename, min_file_bytes, md5_hash=md5_hash, sizeInBytes=sizeInBytes
            ):
                success = False
                if retry > 0:
                    # success = False
                    # this may be needed if initial download fails
                    alt_cmd = download_command.replace(
                        "blocking=true", "blocking=false")
                    current_retry = 1
                    while current_retry <= retry:
                        if current_retry % 2 == 1:
                            retry_cmd = alt_cmd
                        else:
                            retry_cmd = download_command
                        print(
                            "Trying '{}' again due to download error ({}/{}):\n{}"
                            .format(filename, current_retry, retry, retry_cmd)
                        )
                        status = subprocess.run(retry_cmd, shell=True).returncode
                        if status == 0 and not is_broken(
                                filename, min_file_bytes, md5_hash=md5_hash, sizeInBytes=sizeInBytes
                        ):
                            success = True
                            break
                        current_retry += 1
                        time.sleep(10)

        return filename, download_command, success

    def get_regex():
        """
        Get regex pattern from user, compile and return.

        """
        # TODO make this exit gracefully if user can't
        # manage to get a working regex
        compile_success = False
        while compile_success is False:
            pattern = input("Regex pattern: ")
            try:
                pattern = re.compile(pattern)
                compile_success = True
            except:
                print("[!] ERROR: Regex pattern failed to compile.")

        return re.compile(pattern)

    def log_failed(organism, failed_urls):
        """
        Write failed URLs to a local log file.

        """
        fail_log = "{}.failed.log".format(organism)
        print(
            "{} failed downloads logged to {}".format(len(failed_urls), fail_log))
        # write failed URLs to local file
        with open(fail_log, 'w') as f:
            f.write('\n'.join(failed_urls))

    def download_list(url_list, url_to_validate={}, timeout=120, retries=3):
        """
        Attempts download command on a list of partial file
        URLs (completed by download_from_url()).

        Returns a list of successfully-downloaded files and a
        list of unsuccessful URLs

        """
        # Run curl commands to retrieve selected files
        # Make sure the URL formats conforms to the Genome Portal format

        downloaded_files = []
        broken_urls = []
        broken_files = []
        subprocess.run(LOGIN_STRING, shell=True)
        start_time = time.time()
        for url in url_list:
            current_time = time.time()
            # refresh the session cookie every 5 minutes
            if current_time - start_time > 300:
                subprocess.run(LOGIN_STRING, shell=True)
                start_time = time.time()
            fn, cmd, success = download_from_url(
                url, timeout=timeout, retry=retries, url_to_validate=url_to_validate)
            if not success:
                broken_urls.append(url)
                broken_files.append(fn)
            else:
                downloaded_files.append(fn)
        # in cases where multiple files with same name are present and any of them
        # succeed, we can remove corresponding URLs from the list of broken URLs
        # (otherwise, they would just overwrite one another).
        # TODO we could also rename any files with identical names, although then
        # we would need to differentiate between files with different content and
        # files that are just broken versions of the same file...
        broken_urls = [
            u for u, f in zip(broken_urls, broken_files)
            if f not in downloaded_files
        ]

        return downloaded_files, broken_urls

    def get_md5(*fns, buffer_size=65536):
        hash = md5()
        for fn in fns:
            with open(fn, "rb") as f:
                while True:
                    data = f.read(buffer_size)
                    if not data:
                        break
                    hash.update(data)

        return hash.hexdigest()

    def get_sizeInBytes(filename):
        try:
            file_sizeInBytes = os.path.getsize(filename)
        except:
            file_sizeInBytes = 0

        return file_sizeInBytes

    def byte_convert(byte_size):
        """
        Converts a number of bytes to a more human-readable
        format.

        """
        # Calculate and display total size of selected data
        adjusted = byte_size / (1024 * 1024)  # bytes to MB
        if adjusted < 1:
            adjusted = byte_size / 1024
            unit = "KB"
        elif adjusted < 1024:
            unit = "MB"
        else:
            adjusted /= 1024
            unit = "GB"
        size_string = "{:.2f} {}".format(adjusted, unit)
        return size_string

    def get_org_name(xml_file):
        """
        Checks an XML file for organism name information,
        for cases where an XML file is used without organism
        information supplied by the user. Returns None if
        no organism name is found.

        XML entry format is: <organismDownloads name="org_name">

        """
        name_pattern = r"name=\"(.+)\""
        org_line = None
        with open(xml_file) as f:
            for l in f:
                if "organismDownloads" in l:  # standardized name indicator
                    org_line = l.strip()
                    break  # don't keep looking, already found
        try:
            org_name = re.search(name_pattern, org_line).group(1)
            return org_name
        except TypeError:  # org_line still None
            return None

    def parse_selection(user_input):
        """
        Parses the user choice string and returns a dictionary
        of categories (keys) and choices within each category
        (values).

        """
        selections = {}
        parts = user_input.split(";")
        for p in parts:
            if len(p.split(":")) > 2:
                clean_exit("FATAL ERROR: can't parse desired input\n?-->'{}'"
                           .format(p))
            category, indices = p.split(":")
            category = int(category)
            selections[category] = []
            cat_list = selections[category]
            indices = indices.split(",")
            for i in indices:
                try:
                    cat_list.append(int(i))  # if it's already an integer
                except ValueError:
                    try:
                        start, stop = list(map(int, i.split("-")))
                    except:
                        clean_exit("FATAL ERROR: can't parse desired "
                                   "input\n?-->'{}'".format(i))
                    add_range = list(range(start, stop + 1))
                    for e in add_range:
                        cat_list.append(e)
        return selections

    def print_data(data, org_name, display=True):
        """
        Prints info from dictionary data in a specific format.
        Returns a dict with url information for every file
        in desired categories, as well as a dict with md5 information for
        each file (keyed by file URL).

        """
        print("\nQUERY RESULTS FOR '{}'\n".format(org_name))
        dict_to_get = {}
        url_to_validate = defaultdict(dict)
        for query_cat, v in sorted(iter(data.items()),
                                   key=lambda k_v: k_v[1]["catID"]):
            print_list = []
            if not v["results"]:
                continue
            catID = v["catID"]
            dict_to_get[catID] = {}
            print_list.append(" {}: {} ".format(catID, query_cat).center(80, "="))
            results = v["results"]
            for sub_cat, items in sorted(iter(results.items()),
                                         key=lambda sub_cat_items:
                                         (sub_cat_items[0], sub_cat_items[1])):
                print_list.append("{}:".format(sub_cat))
                for index, i in sorted(items.items()):
                    integrity_tag = ""
                    url = i["url"]
                    dict_to_get[catID][index] = url
                    if "md5" in i:
                        url_to_validate[url]["md5"] = i["md5"]
                    # the following elif takes care of MD5 > sizeInBytes rank-order
                    # in downstream processing
                    elif "sizeInBytes" in i:
                        url_to_validate[url]["sizeInBytes"] = int(i["sizeInBytes"])
                    print_index = " {}:[{}] ".format(str(catID), str(index))
                    date = fmt_timestamp(i["timestamp"])
                    date_string = "{:02d}/{}".format(date.tm_mon, date.tm_year)
                    size_date = "[{}|{}]".format(i["size"], date_string)
                    filename = i["filename"]
                    margin = 80 - (len(size_date) + len(print_index))
                    file_info = filename.ljust(margin, "-")
                    print_list.append("".join([print_index, file_info, size_date]))
            if display is True:
                print('\n'.join(print_list))
                print()  # padding

        return dict_to_get, url_to_validate

    def xml_hunt(xml_file):
        """
        Gets list of all XML entries with "filename" attribute,
        and returns a dictionary of the file attributes keyed
        by a ":"-joined string of parent names.

        """
        root = ET.iterparse(xml_file, events=("start", "end"))
        parents = []
        matches = {}
        for event, element in root:
            if element.tag not in ["folder", "file"]:  # skip topmost categories
                continue
            if element.tag == "folder":
                if event == "start":  # add to parents
                    parents.append(element.attrib["name"])
                elif event == "end":  # strip from parents
                    del parents[-1]
                continue
            if event == "start" and element.tag == "file":
                parent_string = ":".join(parents)
                try:
                    matches[parent_string].append(element.attrib)
                except KeyError:
                    matches[parent_string] = [element.attrib]
        return matches

    def format_found(d, filter_found=False):
        """
        Reformats the output from xml_hunt()

        """
        output = {}
        for p, c in sorted(d.items()):
            layers = [e for e in p.split(":") if e]
            if filter_found:
                if not any(cat in layers for cat in DESIRED_CATEGORIES):
                    continue
            if len(layers) == 1:
                top = parent = layers[0]
            else:
                top = layers[-2]  # either -2 or -1 works well, != parent
                parent = layers[-1]  # either -2 or -1 works well, != top
            if top not in output:
                output[top] = defaultdict(dict)
            if parent not in output[top]:
                output[top][parent] = c
            else:
                output[top][parent].extend(c)
        return output

    def get_file_list(xml_file, filter_categories=False):
        """
        Moves through the xml document <xml_file> and returns information
        about matches to elements in <DESIRED_CATEGORIES> if
        <filter_categories> is True, or all files otherwise

        """
        descriptors = {}
        display_cats = ['filename', 'url', 'size',
                        'label', 'sizeInBytes', 'timestamp', 'md5']
        found = xml_hunt(xml_file)
        found = format_found(found, filter_categories)
        if not list(found.values()):
            return None
        category_id = 0
        for category, sub_cat in sorted(found.items()):
            c = category
            if c not in descriptors:
                category_id += 1
                descriptors[c] = defaultdict(dict)
                descriptors[c]["catID"] = category_id
            uid = 1
            for parent, children in sorted(sub_cat.items()):
                descriptors[c]["results"][parent] = defaultdict(dict)
                results = descriptors[c]["results"][parent]
                unique_children = uniqueify(children)
                for child in sorted(unique_children, key=lambda x: x['filename']):
                    try:
                        results[uid]
                    except KeyError:
                        results[uid] = {}
                    for dc in display_cats:
                        try:
                            results[uid][dc] = child[dc]
                        except KeyError:
                            continue
                    uid += 1

        return descriptors

    def uniqueify(children):
        """
        Takes a list of child XML elements (dicts of attribs) as
        returns a filtered list of only unique filenames for a given
        month/year timestamp (e.g. duplicates are allowed if month/year
        is different).

        """
        unique = {}
        for child in children:
            try:
                fn = child['filename']
                date = fmt_timestamp(child['timestamp'])
                date_string = (date.tm_mon, date.tm_year)
                uid = (fn, date_string)
            except KeyError:
                continue
            if fn not in unique:
                unique[uid] = child
            else:
                existing = unique[uid].get('fileType', None)
                if existing == 'Unknown':
                    existing = None
                current = child.get('fileType', None)
                if current == 'Unknown':
                    current = None
                if current is not None and existing is None:
                    unique[uid] = child

        return unique.values()

    def get_sizes(d, sizes_by_url=None):
        """
        Builds a dictionary of url:sizes from
        output of get_file_list()

        """
        for k, v in d.items():
            if isinstance(v, dict):
                if "url" in v:
                    address = v["url"]
                    try:
                        size = int(v["sizeInBytes"])
                    except:
                        size = None
                    sizes_by_url[address] = size
                else:
                    get_sizes(v, sizes_by_url)
        return sizes_by_url

    def extract_file(file_path, keep_compressed=False):
        """
        Native Python file decompression for tar.gz and .gz files.

        TODO: implement .zip decompression

        """
        tar_pattern = "tar.gz$"  # matches tar.gz
        gz_pattern = "(?<!tar)\.gz$"  # excludes tar.gz
        endings_map = {"tar": (tarfile, "r:gz", ".tar.gz"),
                       "gz": (gzip, "rb", ".gz")
                       }
        relative_name = os.path.basename(file_path)
        if re.search(tar_pattern, file_path):
            opener, mode, ext = endings_map["tar"]
            with opener.open(file_path) as f:
                file_count = len(f.getmembers())
                if file_count > 1:  # make sub-directory to unpack into
                    dir_name = relative_name.rstrip(ext)
                    try:
                        os.mkdir(dir_name)
                    except FileExistsError:
                        pass
                    destination = dir_name
                else:  # single file, extract into working directory
                    destination = "."

                def is_within_directory(directory, target):

                    abs_directory = os.path.abspath(directory)
                    abs_target = os.path.abspath(target)

                    prefix = os.path.commonprefix([abs_directory, abs_target])

                    return prefix == abs_directory

                def safe_extract(tar, path=".", members=None, *, numeric_owner=False):

                    for member in tar.getmembers():
                        member_path = os.path.join(path, member.name)
                        if not is_within_directory(path, member_path):
                            raise Exception("Attempted Path Traversal in Tar File")

                    tar.extractall(path, members, numeric_owner=numeric_owner)

                safe_extract(f, destination)
        elif re.search(gz_pattern, file_path):
            opener, mode, ext = endings_map["gz"]
            # out_name = file_path.rstrip(ext)
            out_name = relative_name.rstrip(ext)
            with opener.open(file_path) as f, open(out_name, "wb") as out:
                for l in f:
                    out.write(l)
        else:
            print("Skipped decompression for '{}'"
                  .format(file_path))
            return
        if not keep_compressed:
            os.remove(file_path)

    def decompress_files(local_file_list, keep_original=False):
        """
        Decompresses list of files, and deletes compressed
        copies unless <keep_original> is True.

        """
        for f in local_file_list:
            extract_file(f, keep_original)

    def fmt_timestamp(time_string):
        """
        Parses the timestamp string from an XML document
        of the form "Thu Feb 27 16:38:54 PST 2014"
        and returns a string of the form "2014".

        """
        # Remove platform-dependent timezone substring
        # of the general form "xxT"
        tz_pattern = re.compile("\s[A-Z]{3}\s")
        time_string = tz_pattern.sub(" ", time_string)

        # Get the desired time info
        time_info = time.strptime(time_string, "%a %b %d %H:%M:%S %Y")
        # year = str(time_info.tm_year)
        return time_info