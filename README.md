# climcrawler
Tool to setup the attribution datastore for CDS: download from climexp, rename and create manifest file

## Installation
First install CDO (Climate Data Operators) and netCDF4 manually, either from source or by using your OS package manager (e.g. apt-get on ubuntu). Then create a virtual python-3 environment and install the python dependencies
```shell
virtualenv climcrawl
source climcrawl/bin/activate
pip install -r requirements.txt
```

## Usage
Basic downloading of observations into folder `DIR` and splitting into yearly chunks is achieved by running
```
./climcrawl.py --dir DIR --file csv/obs_datasets.csv --splityrs
```
For all options, run
```shell
./climcrawl.py --help
```
