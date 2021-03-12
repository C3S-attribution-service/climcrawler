#!/usr/bin/env python
import csv
import string
import os
import logging
import enum
import requests
import argparse
import tempfile
import subprocess
import tqdm
import netCDF4 as nc
import datetime
import shutil
from dateutil.parser import parse
from dateutil.relativedelta import relativedelta

file_template = "${source}_${operator}_${variable}_${frequency}_${resolution}_${region}_${start}_${end}_${version}.nc"

info_fields = ["source", "operator", "variable", "frequency", "resolution", "version"]

logging.basicConfig(level=logging.INFO)


class DownloadMode(enum.Enum):
    replace = 1
    skip = 2

log = logging.getLogger(__name__)

read_start_token = "RDSTART"
read_end_token = "RDEND"


def map_row_values(d):
    res = d["resolution"]
    d["resolution"] = res.replace('.', 'd')
    d["start"] = read_start_token
    if d["updates"].lower() in ["yes", "y"]:
        d["end"] = "present"
    elif d["updates"].lower() in ["no", "n"]:
        d["end"] = read_end_token
    else:
        raise ValueError("Updates field {} not recognized: choose from yes,no,y,n".format(d["updates"]))
    return d


def read_filters(sources, vars, freqs):
    def remap(s):
        if s.lower() in ["tmp", "tas", "t"]:
            return ["temperature", "temperature-anomaly"]
        if s.lower() == "pr":
            return ["precipitation"]
        if s == "daily":
            return ["day"]
        if s == "monthly":
            return ["mon"]
        return [s]
    result = {}
    for key,lst in zip(["source", "variable", "frequency"], [sources, vars, freqs]):
        if lst:
            result[key] = set()
            for s in lst.split(','):
                result[key].update(remap(s))
    return {k: list(v) for k, v in result.items() if any(v)}


def read_csv(fname, filters):
    print(filters)
    result = []
    with open(fname) as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
        count = 0
        for row in reader:
            count += 1
            include = True
            for key, values in filters.items():
                if any(values) and row[key] not in values:
                    log.info("Skipping row {} with {} {}, which is outside filters...".format(count,key,row[key]))
                    include = False
                    break
            if include:
                result.append(row)
    return result


def download_file(url, local_filename):
    log.info("Downloading data from %s" % url)
    if os.environ.get("CLIMCRAWL_DRYRUN", "False").lower() == "true":
        return
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total_length = int(r.headers.get('content-length', 0))
            chunk_size = 1024
            progress_bar = tqdm.tqdm(total=total_length, unit='iB', unit_scale=True)
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size):
                    progress_bar.update(len(chunk))
                    f.write(chunk)
            progress_bar.close()
        return local_filename
    except requests.RequestException as e:
        log.error("URL request error encountered: {}".format(e.msg))
        return None



def fix_time_units(ncfile, row):
    ds = nc.Dataset(ncfile, 'r')
    run_options = []
    for varname, variable in ds.variables.items():
        if getattr(variable, "standard_name", "").lower() == "time" or varname.lower() == "time":
            units, freq = variable.units.strip().lower(), row["frequency"]
            substring_index = units.index("since") + 5
            refdate_string = units[substring_index:].strip()
            if units.startswith("years since"):
                log.info("Changing time units for {} {} data...".format(row["source"], row["frequency"]))
                refyear = int(variable[0])
                new_unit = "months" if freq == "mon" else "days"
                ncap2_option = "\'@units=\"{} since {}-01-01 00:00:00\";time=udunits(time,@units);time@units=@units\'".format(new_unit, refyear)
                run_options = ["ncap2", "-O", "-s", ncap2_option]
    ds.close()
    if any(run_options):
        proc = subprocess.run(' '.join(run_options + [ncfile, ncfile]), shell=True)
    return ncfile


def get_time_bounds(ncvar):
    def add_months(refdate, nmonths):
        floor = int(nmonths)
        frac = float(nmonths) - floor
        d1 = refdate + relativedelta(months=floor)
        if frac > 0.:
            d2 = refdate + relativedelta(months=floor+1)
            interval = datetime.timedelta(seconds=frac * (d2 - d1).seconds)
            return d1 + interval
        return d1
    units, calendar = ncvar.units.lower(), getattr(ncvar, "calendar", "proleptic_gregorian")
    if units.startswith("months since"):
        substring_index = units.index("since") + 5
        refdatestr = units[substring_index:].strip()
        refdate = parse(refdatestr)
        return [add_months(refdate, ncvar[0]), add_months(refdate, ncvar[-1])]
    else:
        return nc.num2date([ncvar[0],ncvar[-1]], ncvar.units, calendar=calendar)


def get_file_name(ncfile, row):
    ds = nc.Dataset(ncfile, 'r')
    result = os.path.basename(ncfile)
    for varname, variable in ds.variables.items():
        if getattr(variable, "standard_name", "").lower() == "time" or varname.lower() == "time":
            tempdict = row
            bnds = get_time_bounds(variable)
            tempdict["start"] = bnds[0].strftime("%Y%m") if row["frequency"] == "mon" else bnds[0].strftime("%Y%m%d")
            if row["updates"].lower() in ["yes", "y"]:
                tempdict["end"] = "present"
            else:
                tempdict["end"] = bnds[1].strftime("%Y%m") if row["frequency"] == "mon" else bnds[1].strftime("%Y%m%d")
            result = string.Template(file_template).substitute(tempdict)
    ds.close()
    return result

def make_long_name(row):
    timemap = {"day": "daily", "mon": "monthly", "year": "annual"}
    varmap = {"temperature": "near-surface temperature", "temperature-anomaly": "near-surface temperature anomaly"}
    return ' '.join([row["operator"], 
                    timemap.get(row["frequency"], row["frequency"]), 
                    varmap.get(row["variable"], row["variable"])])


def process_file(ncfile, row):
    log.info("Correcting metadata for file {}".format(ncfile))
    ds = nc.Dataset(ncfile, 'r+')
    rename_tuple = None
    for varname, variable in ds.variables.items():
        if getattr(variable, "standard_name", "").lower() == "time" or varname.lower() == "time":
            bnds = get_time_bounds(variable)
            bnd_strings = [b.strftime("%Y-%m-%d %H:%M:%S") for b in bnds]
            setattr(ds,"time_coverage_start",bnd_strings[0])
            setattr(ds,"time_coverage_end",bnd_strings[1])
        elif getattr(variable, "standard_name", "").lower() not in ["latitude", "longitude", "ensemble member"]:
            if row["variable"] == "temperature":
                setattr(variable, "standard_name", "air_temperature")
                newname = "tas"
                if row["operator"] in ["maximum", "minimum"]:
                    newname = newname + row["operator"][:3]
                rename_tuple = (varname, newname)
                setattr(variable, "long_name", make_long_name(row))
            elif row["variable"] == "temperature-anomaly":
                setattr(variable, "standard_name", "air_temperature_anomaly")
                newname = "tas"
                if row["operator"] in ["maximum", "minimum"]:
                    newname = newname + row["operator"][:3]
                rename_tuple = (varname, newname)
                setattr(variable, "long_name", make_long_name(row))
            elif row["variable"] == "precipitation":
                setattr(variable, "standard_name", "lwe_precipitation_rate")
                if row["frequency"] == "day":
                    setattr(variable, "units", "mm/day")
                elif row["frequency"] == "mon":
                    setattr(variable, "units", "mm/month")
                elif row["frequency"] == "year":
                    setattr(variable, "units", "mm/year")
                rename_tuple = (varname, "pr")
                setattr(variable, "long_name", make_long_name(row))
            else:
                setattr(variable, "long_name", make_long_name(row))
    if rename_tuple is not None:
        ds.renameVariable(*rename_tuple)
    hist = getattr(ds, "history", "")
    setattr(ds, "history", "{}: climcrawler v0.1\n".format(datetime.datetime.now()) + hist)
    ds.close()


def download_files(rows, target_dir, mode, tmpdir=None):
    for row in rows:
        tmpnc = tempfile.NamedTemporaryFile(suffix=".nc", delete=False, dir=tmpdir)
        if download_file(row["url"], tmpnc.name) is None:
            continue
        fix_time_units(tmpnc.name, row)
        fname = get_file_name(tmpnc.name, row)
        target = os.path.join(target_dir, fname)
        if os.path.isfile(target) and mode == DownloadMode.skip:
            log.info("File %s already up to date and will be skipped" % target)
            os.remove(tmpnc.name)
        else:
            shutil.move(tmpnc.name, target)
            process_file(target, row)


def main(args=None):
    if args is None:
        pass
    parser = argparse.ArgumentParser(description="Download climate explorer data, rename files and create manifest "
                                                 "for CDS")
    parser.add_argument("--dir", metavar="DIR", type=str, help="Download directory, will be created if not existing", default=".")
    parser.add_argument("--mode", metavar="MODE", type=str, default="replace", choices=["replace", "skip"],
                        help="MODE:replace|skip whenever file already exists locally")
    parser.add_argument("--source", metavar="SRC", type=str, help="Filter on the source providers")
    parser.add_argument("--var", metavar="tmp|pr|co2", choices=["tmp", "pr", "co2"], type=str, help="Filter on variable (temperature or precipitation)")
    parser.add_argument("--freq", metavar="mon|day", choices=["mon", "day"], type=str, help="Filter on the time series frequency (monthly or daily)")
    parser.add_argument("--tmpdir", metavar="DIR", type=str, default=None, help="Directory for temporary files (default: /tmp)")
    parser.add_argument("--file", metavar="FILE.csv", type=str, help="File (csv) containing the data descriptions",
                        default="./csv/obs_datasets.csv")
    args = parser.parse_args(args)

    input_file = args.file
    if not os.path.isfile(input_file):
        log.fatal("Could not find input file %s" % input_file)

    target_dir = os.path.abspath(args.dir)
    if not os.path.isdir(target_dir):
        os.mkdir(target_dir)

    mode = DownloadMode[args.mode]

    filters = read_filters(args.source, args.var, args.freq)
    download_list = read_csv(input_file, filters)
    download_files(download_list, target_dir=target_dir, mode=mode, tmpdir=args.tmpdir)

if __name__ == "__main__":
    main()
