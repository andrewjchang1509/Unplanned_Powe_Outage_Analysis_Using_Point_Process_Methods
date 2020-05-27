import geopandas as gpd
import pandas as pd
import requests
import zipfile
import json
import io
import numpy as np
##########################################################################
######################## DATA COLLECTING FROM API ########################
##########################################################################

#https://gis.data.ca.gov/datasets/CalEMA::power-outage-incidents?geometry=-93.516%2C-88.438%2C93.516%2C88.438&selectedAttribute=StartDate
response = requests.get(
    "https://opendata.arcgis.com/datasets/439afad071eb4754903906aff1946719_0.geojson")
data = response.json()
new = gpd.GeoDataFrame.from_features(data["features"])
new.crs={'init': 'epsg:4326'}# add crs
new = new.to_crs(epsg=3857) # reproject the crs
new.IncidentId = new.IncidentId.astype(str) # Some incident id contains letters, probably by mistakes
new.FID = new.FID.astype(str)
########### Load in the latest data
 
with open('/Users/andrewchang/Desktop/222/data/latest_data.json', "r") as dataFile:
    latest = json.load(dataFile)
latest = gpd.GeoDataFrame.from_features(latest["features"]) 
#### Merge the newly obtained data and the current latest daa, then drop duplicates  
#### OBJECTID are unique to each incident           
latest = latest.merge(new, how = "outer")
latest.drop_duplicates(subset = "OBJECTID", inplace = True, keep = "last")
latest.reset_index(drop = True, inplace = True)


latest.to_file('/Users/andrewchang/Desktop/222/data/latest_data.json', driver = "GeoJSON")

###########################################################################
########## DOWNLOAD THE CALIFORNIA GIS SHAPE ZIP FILE  ####################
######### ONLY RUN ONCE! AFTER THAT, ONLY RUN THE DATA COLLECTING CHUNK ###
###########################################################################

url = 'https://data.ca.gov/dataset/e212e397-1277-4df3-8c22-40721b095f33/resource/b0007416-a325-4777-9295-368ea6b710e6/download/ca-county-boundaries.zip'
local_path = 'tmp/'
print('Downloading shapefile...')
r = requests.get(url)
z = zipfile.ZipFile(io.BytesIO(r.content))
print("Done")

z.extractall(path=local_path) # extract to folder

######## READ IN SHAPE FILE #########

filenames = [y for y in sorted(z.namelist()) 
             for ending in ['dbf', 'prj', 'shp', 'shx'] 
             if y.endswith(ending)] 
print(filenames)

dbf, prj, shp, shx = [filename for filename in filenames]
ca = gpd.read_file(local_path + shp)
