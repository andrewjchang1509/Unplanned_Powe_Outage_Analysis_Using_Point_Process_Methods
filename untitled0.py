import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import json
import contextily as ctx
import seaborn as sns
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from pyproj import Proj, transform
from scipy import stats
###########################################################################
####### IF THE SHAPE FILE ALREADY EXISTS THEN RUN THE FOLLOWING############
###########################################################################

caBaseMap = gpd.read_file("/Users/andrewchang/Desktop/222/tmp/CA_Counties/CA_Counties_TIGER2016.shp")
crs = "+proj=longlat +datum=WGS84 +no_defs"####Scale the crs
ca = caBaseMap.to_crs(crs=crs) # this is for plotting without basemape

###########################################################################
############################    DATA WRANGLING    #########################
###########################################################################
with open('/Users/andrewchang/Desktop/222/data/latest_data.json', "r") as latest:
    data = json.load(latest)
data = gpd.GeoDataFrame.from_features(data["features"]) 
data.crs={'init': 'epsg:3857'} # assign crs

latest = data[data.OutageType == "Not Planned"].copy()
latest = latest.drop(['EstimatedRestoreDate','Cause', 'OutageStatus', 'OutageType', 'GlobalID', 'OutageTypeColor','OutageStatusColor', 'IncidentId'],axis = 1)

scaler = MinMaxScaler() #default is min = 0, max =  1

# extract coordinates, then scale it to [0,1] range
latest['x'] = scaler.fit_transform(latest.geometry.centroid.x.values.reshape(-1,1))
latest['y'] = scaler.fit_transform(latest.geometry.centroid.y.values.reshape(-1,1))
latest['xp'] = latest.geometry.centroid.x.values.reshape(-1,1)
latest['yp'] = latest.geometry.centroid.y.values.reshape(-1,1)
# create columns for lonlat form
inProj = Proj(init = 'epsg:3857')
outProj = Proj(init = 'epsg:4326')
latest['lon'],latest['lat'] = transform(inProj,outProj,latest.geometry.centroid.x.values,latest.geometry.centroid.y.values)
# process datetime data
latest.StartDate = pd.to_datetime(latest.StartDate, format='%Y-%m-%dT%H:%M:%S.%fZ')
latest.rename(columns={"StartDate":"StartDateTime"}, inplace = True)
latest['StartDate'] = latest.StartDateTime.dt.date
latest['StartTime'] = latest.StartDateTime.dt.time
latest = latest.sort_values(by = "StartDate", ascending = True)

# save to csv file to use with R
latest.copy().drop('geometry',axis=1).to_csv('/Users/andrewchang/Desktop/222/data/check.csv')

# add a column to show the power outage density for each county
caBaseMap["NAME"]=caBaseMap.NAME.str.upper()
caBaseMap["Count"] = caBaseMap.NAME.map(latest.County.value_counts())
caBaseMap["Count"].fillna(0, inplace = True)