#!/usr/bin/env python

import csv

vars = {"tas": "mean temperature",
        "tasmin": "minimum temperature",
        "tasmax": "minimum temperature",
        "pr": "total precipitation",
        "psl": "mean sea-level-pressure",
        "rsds": "mean surface-downward-shortwave-radiation",
        "ssr": "mean surface-shortwave-radiation",
        "evspsbl": "total evaporation",
        "sfcWind": "mean surface-wind-speed",
        "sfcWindmax": "maximum surface-wind-speed",
        "uas": "mean eastward-surface-wind-speed",
        "vas": "mean northward-surface-wind-speed"}

freqs = {   "day": "day",
            "mon": "Amon"}

header = [  "source",
            "version",
            "provider",
            "variable",
            "operator",
            "frequency",
            "resolution",
            "region",
            "updates",
            "experiment",
            "member",
            "url"]

def main(args=None):
    if args is None:
        pass
    with open("mod_ece23.csv", mode='w') as ofile:
            writer = csv.writer(ofile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(header)
            source = "ECEARTH23"
            version = "1.0"
            provider = "EC-Earth-Consortium"
            resolution = "T159L62"
            region="global"
            updates="no"
            experiment="rcp85"
            top_url = "https://climexp.knmi.nl/KNMI14Data/CMIP5/output/KNMI/{0}/{1}/".format(source, experiment)
            for freqname, freqtag in freqs.items():
                frequency = freqname
                for varname, vardescr in vars.items():
                    if varname == "sfcWindmax" and frequency == "mon":
                        continue
                    operator = vardescr.split(' ')[0]
                    variable = vardescr.split(' ')[1]
                    for mem in range(1,17):
                        member = "r{0}i1p1".format(mem)
                        bot_url = "{0}/atmos/{1}/{2}/v1/{3}/".format(freqname, freqtag, member, varname)
                        time_spans = ["186001-210012"]
                        if frequency == "day":
                            time_spans = [  "18600101-18991231",
                                            "19000101-19501231",
                                            "19510101-20001231",
                                            "20010101-20501231",
                                            "20510101-21001231"]
                        for time_span in time_spans:
                            fname = "{0}_{1}_{2}_{3}_{4}_{5}.nc".format(varname, freqtag, source, experiment, member, time_span)
                            row = [ source,
                                    version,
                                    provider,
                                    variable,
                                    operator,
                                    frequency,
                                    resolution,
                                    region,
                                    updates,
                                    experiment,
                                    member,
                                    top_url + bot_url + fname]
                            writer.writerow(row)

if __name__ == "__main__":
    main()
