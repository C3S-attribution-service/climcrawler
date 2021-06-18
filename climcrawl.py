#!/usr/bin/env python
import argparse
import datetime
import enum
import logging
import os
import shutil
import string
import tempfile

import cdo
import netCDF4 as nc
import netcdftime
import numpy as np
import requests
import tqdm
from dateutil.parser import parse
from dateutil.relativedelta import relativedelta

import csv

obs_file_template = "${source}_${operator}_${variable}_${frequency}_${resolution}_${region}_${start}_${end}_${" \
                    "version}.nc"
mod_file_template = "${source}_${operator}_${variable}_${frequency}_${resolution}_${region}_${experiment}_${" \
                    "member}_${start}_${end}_${version}.nc"

info_fields = ["source", "operator", "variable", "frequency", "resolution", "version", "member", "experiment"]

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
    def remap(cell):
        if cell.lower() in ["tmp", "tas", "t"]:
            return ["temperature", "temperature-anomaly"]
        if cell.lower() == "pr":
            return ["precipitation"]
        if cell == "daily":
            return ["day"]
        if cell == "monthly":
            return ["mon"]
        return [cell]

    result = {}
    for key, lst in zip(["source", "variable", "frequency"], [sources, vars, freqs]):
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
                    log.info("Skipping row {} with {} {}, which is outside filters...".format(count, key, row[key]))
                    include = False
                    break
            if include:
                result.append(row)
    return result


def download_file(url, local_filename):
    log.info("Downloading data from %s" % url)
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total_length = int(r.headers.get('content-length', 0))
            chunk_size = 1024
            if os.environ.get("CLIMCRAWL_DRYRUN", "False").lower() == "true":
                log.info("Checked url {0}: OK".format(url))
                return None
            progress_bar = tqdm.tqdm(total=total_length, unit='iB', unit_scale=True)
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size):
                    progress_bar.update(len(chunk))
                    f.write(chunk)
            progress_bar.close()
        return local_filename
    except requests.RequestException as e:
        log.error("URL request error encountered: {}".format(str(e)))
        return None


def fix_time_units(ncfile, row):
    ds = nc.Dataset(ncfile, 'r+')
    for varname, variable in ds.variables.items():
        if getattr(variable, "standard_name", "").lower() == "time" or varname.lower() == "time":
            units, freq = variable.units.strip().lower(), row["frequency"]
            #            substring_index = units.index("since") + 5
            if units.split(' ')[0] in ['years', 'months']:
                log.info("Changing time units for {} {} data...".format(row["source"], row["frequency"]))
                substring_index = units.index("since") + 5
                refyearstr = units[substring_index:].strip().split(' ')[0].split('-')[0]
                refyear = int(refyearstr)
                subtract_year = (refyear == 0)
                if subtract_year:
                    units = units[:substring_index] + " " + "0001" + units[substring_index + 1 + len(refyearstr):]
                    log.info("My new units are {0}".format(units))
                dates = num2date(variable[:], units, calendar=getattr(variable, "calendar", "Standard"))
                if subtract_year:
                    dates = [add_year(d, -1) for d in dates]
                refyear = dates[0].year
                new_unit = "days since {}-01-01 00:00:00".format(refyear)
                variable[:] = date2num(dates, new_unit, calendar=getattr(variable, "calendar", "Standard"))
                variable.units = new_unit
    ds.close()
    return ncfile


def get_time_bounds(ncvar):
    def add_months(rdate, nmonths):
        floor = int(nmonths)
        frac = float(nmonths) - floor
        d1 = rdate + relativedelta(months=floor)
        if frac > 0.:
            d2 = rdate + relativedelta(months=floor + 1)
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
        return nc.num2date([ncvar[0], ncvar[-1]], ncvar.units, calendar=calendar)


def get_file_name(ncfile, row, name_template, split_years):
    if split_years:
        short_template = '_'.join(name_template.split('_')[:-3])
        return string.Template(short_template).substitute(row)
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
            result = string.Template(name_template).substitute(tempdict)
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
            setattr(ds, "time_coverage_start", bnd_strings[0])
            setattr(ds, "time_coverage_end", bnd_strings[1])
        elif getattr(variable, "standard_name", "").lower() not in ["latitude", "longitude", "ensemble member",
                                                                    "height"]:
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
                setattr(variable, "units", "degrees Celsius")
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
    if rename_tuple is not None and rename_tuple[0] != rename_tuple[1]:
        ds.renameVariable(*rename_tuple)
    hist = getattr(ds, "history", "")
    setattr(ds, "history", "{}: climcrawler v0.1\n".format(datetime.datetime.now()) + hist)
    ds.close()


def download_files(rows, target_dir, mode, name_template, split_yrs=False, tmpdir=None):
    for row in rows:
        tmpnc = tempfile.NamedTemporaryFile(suffix=".nc", delete=False, dir=tmpdir)
        if download_file(row["url"], tmpnc.name) is None:
            continue
        fix_time_units(tmpnc.name, row)
        fname = get_file_name(tmpnc.name, row, name_template, split_yrs)
        target = os.path.join(target_dir, fname)
        if os.path.isfile(target) and mode == DownloadMode.skip:
            log.info("File %s already up to date and will be skipped" % target)
            os.remove(tmpnc.name)
        else:
            if split_yrs:
                process_file(tmpnc.name, row)
                filelist = split_years(tmpnc.name, target + '_')
                for filepath in filelist:
                    shutil.move(filepath, filepath[:-3] + "_" + row["version"] + ".nc")
            else:
                shutil.move(tmpnc.name, target)
                process_file(target, row)


def split_years(fname, target):
    cdoapp = cdo.Cdo()
    cdoapp.debug = True
    return cdoapp.splityear(input=fname, output=target)


def num2date(num_axis, units, calendar):
    """
    A wrapper from ``netCDF4.num2date`` able to handle "years since" and "months since" units.
    If time units are not "years since" or "months since", calls usual ``netcdftime.num2date``.

    :param numpy.array num_axis: The numerical time axis following units
    :param str units: The proper time units
    :param str calendar: The NetCDF calendar attribute
    :returns: The corresponding date axis
    :rtype: *array*

    """
    if not units.split(' ')[0] in ['years', 'months']:
        # If units are not 'years' or 'months since', call usual netcdftime.num2date:
        return nc.num2date(num_axis, units=units, calendar=calendar)
    else:
        # Return to time reference with 'days since'
        units_as_days = 'days ' + ' '.join(units.split(' ')[1:])
        # Convert the time reference 'units_as_days' as datetime object
        start_date = nc.num2date(0.0, units=units_as_days, calendar=calendar)
        # Control num_axis to always get an Numpy array (even with a scalar)
        num_axis_mod = np.atleast_1d(np.array(num_axis))
        if units.split(' ')[0] == 'years':
            # If units are 'years since'
            # Define the number of maximum and minimum years to build a date axis covering
            # the whole 'num_axis' period
            max_years = np.floor(np.max(num_axis_mod)) + 1
            min_years = np.ceil(np.min(num_axis_mod)) - 1
            # Create a date axis with one year that spans the entire period by year
            years_axis = np.array([add_year(start_date, years_to_add)
                                   for years_to_add in np.arange(min_years, max_years + 2)])
            # Convert rebuilt years axis as 'number of days since'
            cdftime = netcdftime.utime(units_as_days, calendar=calendar)
            years_axis_as_days = cdftime.date2num(years_axis)
            # Index of each years
            yind = np.vectorize(np.int)(np.floor(num_axis_mod))
            # Rebuilt num_axis as 'days since' adding the number of days since referenced time
            # with an half-increment (num_axis_mod - yind) = 0 or 0.5
            num_axis_mod_days = (years_axis_as_days[yind - int(min_years)] +
                                 (num_axis_mod - yind) *
                                 np.diff(years_axis_as_days)[yind - int(min_years)])
            # Convert result as date axis
            return nc.num2date(num_axis_mod_days, units=units_as_days, calendar=calendar)
        elif units.split(' ')[0] == 'months':
            # If units are 'months since'
            # Define the number of maximum and minimum months to build a date axis covering
            # the whole 'num_axis' period
            max_months = np.floor(np.max(num_axis_mod)) + 1
            min_months = np.ceil(np.min(num_axis_mod)) - 1
            # Create a date axis with one month that spans the entire period by month
            months_axis = np.array([add_month(start_date, months_to_add)
                                    for months_to_add in np.arange(min_months, max_months + 12)])
            # Convert rebuilt months axis as 'number of days since'
            cdftime = netcdftime.utime(units_as_days, calendar=calendar)
            months_axis_as_days = cdftime.date2num(months_axis)
            # Index of each months
            mind = np.vectorize(np.int)(np.floor(num_axis_mod))
            # Rebuilt num_axis as 'days since' adding the number of days since referenced time
            # with an half-increment (num_axis_mod - mind) = 0 or 0.5
            num_axis_mod_days = (months_axis_as_days[mind - int(min_months)] +
                                 (num_axis_mod - mind) *
                                 np.diff(months_axis_as_days)[mind - int(min_months)])
            # Convert result as date axis
            return nc.num2date(num_axis_mod_days, units=units_as_days, calendar=calendar)


# Copied from nctime source code
def date2num(date_axis, units, calendar):
    """
    A wrapper from ``netCDF4.date2num`` able to handle "years since" and "months since" units.
    If time units are not "years since" or "months since" calls usual ``netcdftime.date2num``.

    :param numpy.array date_axis: The date axis following units
    :param str units: The proper time units
    :param str calendar: The NetCDF calendar attribute
    :returns: The corresponding numerical time axis
    :rtype: *array*

    """
    # date_axis is the date time axis incremented following units (i.e., by years, months, etc).
    if not units.split(' ')[0] in ['years', 'months']:
        # If units are not 'years' or 'months since', call usual netcdftime.date2num:
        return nc.date2num(date_axis, units=units, calendar=calendar)
    else:
        # Return to time reference with 'days since'
        units_as_days = 'days ' + ' '.join(units.split(' ')[1:])
        # Convert date axis as number of days since time reference
        days_axis = nc.date2num(date_axis, units=units_as_days, calendar=calendar)
        # Convert the time reference 'units_as_days' as datetime object
        start_date = nc.num2date(0.0, units=units_as_days, calendar=calendar)
        # Create years axis from input date axis
        years = np.array([date.year for date in np.atleast_1d(np.array(date_axis))])
        if units.split(' ')[0] == 'years':
            # If units are 'years since'
            # Define the number of maximum and minimum years to build a date axis covering
            # the whole 'num_axis' period
            max_years = np.max(years - start_date.year + 1)
            min_years = np.min(years - start_date.year - 1)
            # Create a date axis with one year that spans the entire period by year
            years_axis = np.array([add_year(start_date, yid)
                                   for yid in np.arange(min_years, max_years + 2)])
            # Convert years axis as number of days since time reference
            cdftime = netcdftime.utime(units_as_days, calendar=calendar)
            years_axis_as_days = cdftime.date2num(years_axis)
            # Find closest index for years_axis_as_days in days_axis
            closest_index = np.searchsorted(years_axis_as_days, days_axis)
            # Compute the difference between closest value of year axis and start date, in number of days
            num = days_axis - years_axis_as_days[closest_index]
            # Number of days of the corresponding closest year
            den = np.diff(years_axis_as_days)[closest_index]
            return min_years + closest_index + num / den
        elif units.split(' ')[0] == 'months':
            # If units are 'months since'
            # Define the number of maximum and minimum months to build a date axis covering
            # the whole 'num_axis' period
            max_months = np.max(12 * (years - start_date.year + 12))
            min_months = np.min(12 * (years - start_date.year - 12))
            # Create a date axis with one month that spans the entire period by month
            months_axis = np.array([add_month(start_date, mid)
                                    for mid in np.arange(min_months, max_months)])
            # Convert months axis as number of days since time reference
            cdftime = netcdftime.utime(units_as_days, calendar=calendar)
            months_axis_as_days = cdftime.date2num(months_axis)
            # Find closest index for months_axis_as_days in days_axis
            closest_index = np.searchsorted(months_axis_as_days, days_axis)
            # Compute the difference between closest value of months axis and start date, in number of days
            num = days_axis - months_axis_as_days[closest_index]
            # Number of days of the corresponding closest month
            den = np.diff(months_axis_as_days)[closest_index]
            return min_months + closest_index + num / den


# Copied from nctime source code
def add_month(date, months_to_add):
    """
    Finds the next month from date.

    :param netcdftime.datetime date: Accepts datetime or phony datetime from ``netCDF4.num2date``.
    :param int months_to_add: The number of months to add to the date
    :returns: The final date
    :rtype: *netcdftime.datetime*

    """
    years_to_add = int((date.month + months_to_add - np.mod(date.month + months_to_add - 1, 12) - 1) / 12)
    new_month = int(np.mod(date.month + months_to_add - 1, 12)) + 1
    new_year = date.year + years_to_add
    date_next = netcdftime.datetime(year=new_year,
                                    month=new_month,
                                    day=date.day,
                                    hour=date.hour,
                                    minute=date.minute,
                                    second=date.second)

    return date_next


# Copied from nctime source code
def add_year(date, years_to_add):
    """
    Finds the next year from date.

    :param netcdftime.datetime date: Accepts datetime or phony datetime from ``netCDF4.num2date``.
    :param int years_to_add: The number of years to add to the date
    :returns: The final date
    :rtype: *netcdftime.datetime*

    """
    new_year = date.year + years_to_add
    date_next = netcdftime.datetime(year=new_year,
                                    month=date.month,
                                    day=date.day,
                                    hour=date.hour,
                                    minute=date.minute,
                                    second=date.second)
    return date_next


def main(args=None):
    if args is None:
        pass
    parser = argparse.ArgumentParser(description="Download climate explorer data, rename files and create manifest "
                                                 "for CDS")
    parser.add_argument("--dir", metavar="DIR", type=str, help="Download directory, will be created if not existing",
                        default=".")
    parser.add_argument("--mode", metavar="MODE", type=str, default="replace", choices=["replace", "skip"],
                        help="MODE:replace|skip whenever file already exists locally")
    parser.add_argument("--source", metavar="SRC", type=str, help="Filter on the source providers")
    parser.add_argument("--var", metavar="tmp|pr|co2", choices=["tmp", "pr", "co2"], type=str,
                        help="Filter on variable (temperature or precipitation)")
    parser.add_argument("--freq", metavar="mon|day", choices=["mon", "day"], type=str,
                        help="Filter on the time series frequency (monthly or daily)")
    parser.add_argument("--tmpdir", metavar="DIR", type=str, default=None,
                        help="Directory for temporary files (default: /tmp)")
    parser.add_argument("--splityrs", action="store_true", help="Split files into yearly chunks")
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

    if os.path.basename(input_file).startswith('mod'):
        name_template = mod_file_template
    else:
        name_template = obs_file_template

    download_files(download_list, target_dir=target_dir, mode=mode, name_template=name_template,
                   split_yrs=args.splityrs, tmpdir=args.tmpdir)


if __name__ == "__main__":
    main()
