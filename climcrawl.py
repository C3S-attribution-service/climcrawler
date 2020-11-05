#!/usr/bin/env python
import csv
import string
import os
import logging
import enum
import requests
import argparse

file_template = "${source}_${operator}_${variable}_${frequency}_${resolution}_${start}_${end}.nc"

info_fields = ["source", "operator", "variable", "frequency", "resolution"]


class DownloadMode(enum.Enum):
    replace = 1
    update = 2
    skip = 3


log = logging.getLogger(__name__)


def map_row_values(d):
    res = d["resolution"]
    d["resolution"] = res.replace('.', 'd')
    return d


def read_csv(fname):
    result = {}
    token_values = {s: set() for s in info_fields}
    with open(fname) as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        template = string.Template(file_template)
        for row in reader:
            try:
                tokens = map_row_values(row)
                for s in info_fields:
                    token_values[s].add(tokens[s])
                target_file = template.substitute(tokens)
            except KeyError:
                log.fatal("Badly formatted input csv file %s detected, exiting" % fname)
                return {}
            result[target_file] = row.get("url", None)
    return result, token_values


def download_file(url, local_filename):
    log.info("Downloading data from %s" % url)
    if os.environ.get("CLIMCRAWL_DRYRUN", "False").lower() == "true":
        return
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename


def download_files(files, target_dir, mode):
    flist = []
    for fname, source_url in files.items():
        target = os.path.join(target_dir, fname)
        if os.path.isfile(target):
            if mode == DownloadMode.replace:
                log.info("File %s already exists and will be replaced" % target)
                os.remove(target)
                download_file(source_url, target)
            elif mode == DownloadMode.skip:
                log.info("File %s already exists and will be skipped" % target)
            elif mode == DownloadMode.update:
                if target.split('_')[-1].split('.')[0] == "now":
                    log.info("File %s will be updated" % target)
                    download_file(source_url, target)
                else:
                    log.info("File %s already up to date and will be skipped" % target)
        else:
            log.info("File %s not found and will be downloaded" % target)
            download_file(source_url, target)
        flist.append(fname)
    return flist


def write_manifest(flist, base_url, filepath):
    with open(filepath, 'w') as f:
        f.writelines([base_url + f + '\n' for f in flist])


def main(args=None):
    if args is None:
        pass
    parser = argparse.ArgumentParser(description="Download climate explorer data, rename files and create manifest "
                                                 "for CDS")
    parser.add_argument("--dir", metavar="DIR", type=str, help="Download directory, will be created if not existing"
                        , default=".")
    parser.add_argument("--url", metavar="URL", type=str, help="Base url to be prepended in the manifest file",
                        default="https://attribution.climate.copernicus.eu/")
    parser.add_argument("--mode", metavar="MODE", type=str, default="skip", choices=["replace", "skip", "update"],
                        help="MODE:replace|skip|update whenever file already exists locally")
    parser.add_argument("file", metavar="FILE.csv", type=str, help="File (csv) containing the data descriptions",
                        default="./csv/obs_datasets.csv")
    args = parser.parse_args(args)

    input_file = args.file
    if not os.path.isfile(input_file):
        log.fatal("Could not find input file %s" % input_file)

    target_dir = os.path.abspath(args.dir)
    if not os.path.isdir(target_dir):
        os.mkdir(target_dir)

    mode = DownloadMode[args.mode]

    mapping, token_values = read_csv(input_file)
    flist = download_files(mapping, target_dir=target_dir, mode=mode)

    write_manifest(flist, args.url, 'MANIFEST.txt')

    for field, valueset in token_values.items():
        print("All possible values for %s: %s" % (field, ",".join(sorted(valueset))))


if __name__ == "__main__":
    main()