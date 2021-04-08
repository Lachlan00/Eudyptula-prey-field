"""
Acoustic Survey Data Sorting
2019-01-29

This script crawls through the acoustic survey data on
the master data drive and sorts all the raw data by survey
(based on datetime stamps) and then by filetypev (e.g. 
.xml, .raw and calibration files etc.)

Requires
    - input directory
    - output directory
    - dictionary of survey dates
        * produced from CTD metadata

Notes:
    - Raw acoustic filename timestamps are in UTC
    - best way to sort file dates is by nearest survey min/max timestamp 
      to account for anomalies
    - all .raw filenames contain 'DYYYYmmdd-THHMMSS' datetime strings
"""

import os
import re
import numpy as np
import pandas as pd
from datetime import datetime
from progressbar import ProgressBar
from pathlib import Path
from shutil import copyfile

# __Setup__

# directories
input_dir = '/Volumes/LP_MstrData/master-data/survey/acoustics/acoustics_raw'
output_dir = '/Volumes/LP_MstrData/master-data/survey/acoustics/acoustics_sorted'
# CTD meta to produce survey data
CTD_meta_dir = '../../data/CTD/processed/meta/cast_meta.csv'

# __Survey Dates__

# read in data
CTD_meta = pd.read_csv(CTD_meta_dir, parse_dates=['cast_time_UTC']).loc[:,['cast_time_UTC','survey_id']]
survey_min = CTD_meta.groupby('survey_id').min()
survey_max = CTD_meta.groupby('survey_id').max()
# make into new dataframe
survey_dates = survey_min.rename({'cast_time_UTC': 'dt_min'}, axis='columns')
survey_dates['dt_max'] = survey_max.cast_time_UTC
# make combined dataframe
survey_date_marks = pd.DataFrame({
    'survey':np.tile(survey_dates.index,2),
    'dt_UTC':pd.concat([survey_dates.dt_min, survey_dates.dt_max])})
survey_date_marks = survey_date_marks.reset_index(drop=True)

# __Create directories__

# survey level
for survey in survey_dates.index:
    Path(output_dir+'/'+survey+'/raw/transect').mkdir(parents=True, exist_ok=True)
    Path(output_dir+'/'+survey+'/raw/calibration').mkdir(parents=True, exist_ok=True)
    Path(output_dir+'/'+survey+'/raw/self_noise').mkdir(parents=True, exist_ok=True)
    Path(output_dir+'/'+survey+'/xml/transect').mkdir(parents=True, exist_ok=True)
    Path(output_dir+'/'+survey+'/xml/calibration').mkdir(parents=True, exist_ok=True)
    Path(output_dir+'/'+survey+'/xml/software').mkdir(parents=True, exist_ok=True)

# __File crawling__

# Obtain list of all files in directories
file_ls = [os.path.join(dp, f) for dp, dn, filenames in os.walk(input_dir) for f in filenames]
# remove all the OSX '._' and .DS_Store file turds
file_ls = [file for file in file_ls if os.path.basename(file)[0:2] != '._']
file_ls = [file for file in file_ls if os.path.basename(file) != '.DS_Store']

# .raw
# get .raw files
files_raw = pd.DataFrame({'filepath':[file for file in file_ls if file.endswith(".raw")]})
files_raw['basename'] = [os.path.basename(file) for file in files_raw.filepath]
# remove any duplicate base names (due to old data being copied into new directories)
files_raw = files_raw.loc[~files_raw.duplicated(['basename'])].reset_index(drop=True)

# get timstamps for all the .raw files
files_raw['dt_UTC'] = [re.search(r'D\d{8}-T\d{6}', file).group(0) for file in files_raw.basename]
files_raw['dt_UTC'] = [datetime.strptime(t, 'D%Y%m%d-T%H%M%S') for t in files_raw.dt_UTC]
# get what survey data falls into
print('\nFinding survey period for each ".raw" survey file..')
pbar = ProgressBar(max_value=len(files_raw))
pbar.start()
for i, file in files_raw.iterrows():
    # find nearest stamp
    nearest_dt = min(survey_date_marks.dt_UTC, key=lambda x: abs(x - file.dt_UTC))
    files_raw.loc[i,'survey'] = survey_date_marks.loc[survey_date_marks.dt_UTC == nearest_dt].survey.values
    pbar.update(i)
pbar.finish()
# is .raw file a calibration file
files_raw['calibration'] = False
calibration_test = [re.search(r'cal', file, re.IGNORECASE) for file in files_raw.filepath]
calibration_test = [idx for idx, match in enumerate(calibration_test) if match is not None]
files_raw.loc[calibration_test, 'calibration'] = True
# is .raw file a "self noise" file
files_raw['self_noise'] = False
self_noise_test = [re.search(r'PASSIVE', file) for file in files_raw.filepath]
self_noise_test = [idx for idx, match in enumerate(self_noise_test) if match is not None]
files_raw.loc[self_noise_test, 'self_noise'] = True

# .xml
# get .xml files
files_xml = pd.DataFrame({'filepath':[file for file in file_ls if file.endswith(".xml")]})
files_xml['basename'] = [os.path.basename(file) for file in files_xml.filepath]
# duplicates
xml_dups = files_xml.duplicated(['basename'])
# remove the first duplicate of xml0_env.xml as it actually isn't a duplicate (confirmed manually)
xml_dups[2] = False
files_xml = files_xml.loc[~xml_dups].reset_index(drop=True)
# get timstamps for all the .xml files
files_xml['dt_UTC'] = [re.search(r'D\d{8}-T\d{6}', file) for file in files_xml.basename]
files_xml['dt_UTC'] = [match.group(0) if not match is None else None for match in files_xml.dt_UTC]
files_xml['dt_UTC'] = [datetime.strptime(t, 'D%Y%m%d-T%H%M%S') if not t is None else None for t in files_xml.dt_UTC]
# get what survey data falls into
print('\nFinding survey period for each ".xml" survey file..')
pbar = ProgressBar(max_value=len(files_xml))
pbar.start()
for i, file in files_xml.iterrows():
    # find nearest stamp
    if file.dt_UTC is None or pd.isnull(file.dt_UTC):
        files_xml.loc[i,'survey'] = None
    else:
        nearest_dt = min(survey_date_marks.dt_UTC, key=lambda x: abs(x - file.dt_UTC))
        files_xml.loc[i,'survey'] = survey_date_marks.loc[survey_date_marks.dt_UTC == nearest_dt].survey.values
    pbar.update(i)
pbar.finish()
# find surveys for those with no datetime codes by looking for year strings in directory
dir_year = [re.search(r'\d{4}', file).group(0) for file in files_xml.loc[pd.isnull(files_xml.survey)].filepath]
dir_year = [year + '_S1' for year in dir_year]
files_xml.loc[pd.isnull(files_xml.survey), 'survey'] = dir_year
# is .xml file a calibration file
files_xml['calibration'] = False
calibration_test = [re.search(r'cal', file, re.IGNORECASE) for file in files_xml.filepath]
calibration_test = [idx for idx, match in enumerate(calibration_test) if match is not None]
files_xml.loc[calibration_test, 'calibration'] = True
# is .xml file a software file
files_xml['software'] = False
software_test = [re.search(r'xml0', file) for file in files_xml.filepath]
software_test = [idx for idx, match in enumerate(software_test) if match is not None]
files_xml.loc[software_test, 'software'] = True

# __Generate new file path destinations__

# .raw
files_raw['destination'] = ['calibration' if row.calibration else 'self_noise' 
            if row.self_noise else 'transect' for i, row in files_raw.iterrows()]
files_raw['newpath'] = output_dir+'/'+files_raw.survey+'/raw/'+files_raw.destination+'/'+files_raw.basename

# .xml
files_xml['destination'] = ['calibration' if row.calibration else 'software' 
            if row.software else 'transect' for i, row in files_xml.iterrows()]
files_xml['newpath'] = output_dir+'/'+files_xml.survey+'/xml/'+files_xml.destination+'/'+files_xml.basename

# __Copy acorss new files__

# .raw
print('\nCopying ".raw" data files to new directory..')
pbar = ProgressBar(max_value=len(files_raw))
pbar.start()
for i, row in files_raw.iterrows():
    copyfile(row.filepath, row.newpath)
    pbar.update(i)
pbar.finish()

# .xml
print('\nCopying ".raw" data files to new directory..')
pbar = ProgressBar(max_value=len(files_xml))
pbar.start()
for i, row in files_xml.iterrows():
    copyfile(row.filepath, row.newpath)
    pbar.update(i)
pbar.finish()

print('Done!')




