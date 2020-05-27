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
###########################################################################
###################    DATA VISUALIZATIONS/ANALYSIS   #####################
###########################################################################

fig,ax = plt.subplots(figsize = (12,12))
caBaseMap.plot(ax=ax, alpha=0.5, edgecolor='black', column = "Count", cmap='Blues')
latest.plot(ax = ax, color = "red", markersize = 0.2)
ctx.add_basemap(ax) # url=ctx.sources. can be added to change the background
ax.set_axis_off()
# Create colorbar as a legend
vmin, vmax = caBaseMap["Count"].min(),caBaseMap["Count"].max()
sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=vmin, vmax=vmax))
# empty array for the data range
sm._A = []
# add the colorbar to the figure
cbar = fig.colorbar(sm, alpha=0.5) # same alpha as in caBaseMap.plot(.)
# add a title
ax.set_title('Unplanned Outages in California from April to June, 2020', fontdict={'fontsize': '22', 'weight' : 'bold','fontname': 'Times New Roman'})
# create an annotation for the data source
ax.annotate('Data Source: California State Geoportal\nColor bar shows the density of unplanned power outage per county',xy=(0.1, .08),
            xycoords='figure fraction', horizontalalignment='left',
            verticalalignment='top', fontsize=12, color='#555555')
plt.show()
fig.savefig("/Users/andrewchang/Desktop/222/plots/general.png")
### Kernel Smoothing, serves for identifying clustering pattern  ###

fig_k,ax_k = plt.subplots(figsize = (12,12))
sns.kdeplot(latest.geometry.centroid.x.values,latest.geometry.centroid.y.values,\
            shade=True, cmap='Greys', ax=ax_k, bw = "silverman") #silverman bandwidth
caBaseMap.plot(ax=ax_k, edgecolor='black', color = None, alpha = 0.1)
ax_k.set_axis_off()
ax_k.set_title('Density Map of Unplanned Outages in California', fontdict={'fontsize': '22', 'weight' : 'bold','fontname': 'Times New Roman'})
plt.show()
fig_k.savefig("/Users/andrewchang/Desktop/222/plots/pykern.png")

### Kernel Smoothing, using scaled x,y.  ###
xx, yy = np.mgrid[0:1:.01, 0:1:.01]
kde = stats.gaussian_kde((latest[['x','y']]).T, bw_method = 'silverman')
density = kde(np.c_[xx.flat, yy.flat].T).reshape(xx.shape)
fig_s, ax_s = plt.subplots()
cset = ax_s.contourf(xx, yy, density, cmap="Greys")
fig_s.colorbar(cset)
ax_s.set_xlabel("Scaled Longtitude")
ax_s.set_ylabel("Scaled Latitude")
plt.show()
fig_s.savefig("/Users/andrewchang/Desktop/222/plots/pykernScaled.png")
