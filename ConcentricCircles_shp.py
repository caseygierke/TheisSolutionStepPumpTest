# ConcentricCircles_shp.py

# This program takes point inputs and creates 
# concentric circles around them as is specified 
# in a contours file.  The circles represent
# potentiometric surface contours derived from 
# a Theis solution analysis of pumping FCI's PTF
# wells.

# Scripted by Casey Gierke of LWA on
# 6/18/2018

# This code relied heavily on material found at 
# https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-a-new-shapefile-and-add-data

# With Notepad++, use F5 then copy this into box
# C:\Python27\python.exe -i "$(FULL_CURRENT_PATH)"

# ------------------------------------------------------
# IMPORTS
# ------------------------------------------------------

import os
import csv
import osgeo.ogr as ogr
import osgeo.osr as osr

# ------------------------------------------------------
# DEFINE FUNCTIONS
# ------------------------------------------------------

# Define last position finder
def find_last(s,t):
	last_pos = -1
	while True:
		pos = s.find(os.sep, last_pos +1)
		if pos == -1:
			return last_pos
		last_pos = pos

# Create finder function
def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

# ------------------------------------------------------
# INPUTS
# ------------------------------------------------------

# Define which well you are analyzing
# Well = 'R-01'
Well = 'R-03'
# Well = 'R-05'
# Well = 'R-07'

# Define path
path = os.path.abspath(os.path.dirname(__file__))

# Read in contour info
Contours = open(path+os.sep+Well+'- Outfile- Contours.txt','r')

# Create a well dictionary
Wells = {
		"R-01": [647905.525,746225.951], 
		"R-03": [648046.814,746084.337], 
		"R-05": [647905.594,745942.885], 
		"R-07": [647763.998,746084.498], 
		}

# Set up the shapefile driver
driver = ogr.GetDriverByName("ESRI Shapefile")

# Create the spatial reference, WGS84
srs = osr.SpatialReference()
srs.ImportFromEPSG(26749)

# Open a loop to move through contours file
i = 0
for line in Contours:
	x1 = find_nth(line, '\t', 1)
	x2 = find_nth(line, '\t', 2)
	DD = float(line[:x1])
	bufferDistance = float(line[x1+1:-1])
	
	# Delete files so that they can be written
	if os.path.exists(path+os.sep+'Buffer Files'+os.sep+Well+"_Buffer-"+str(DD)+"ft.shp"):
		os.remove(path+os.sep+'Buffer Files'+os.sep+Well+"_Buffer-"+str(DD)+"ft.dbf")
		os.remove(path+os.sep+'Buffer Files'+os.sep+Well+"_Buffer-"+str(DD)+"ft.prj")
		os.remove(path+os.sep+'Buffer Files'+os.sep+Well+"_Buffer-"+str(DD)+"ft.shp")
		os.remove(path+os.sep+'Buffer Files'+os.sep+Well+"_Buffer-"+str(DD)+"ft.shx")

		
	# create the data source
	outFile = driver.CreateDataSource(path+os.sep+'Buffer Files'+os.sep+Well+"_Buffer-"+str(DD)+"ft.shp")
	
	# create the layer
	layer = outFile.CreateLayer(str(DD)+"ft", srs, ogr.wkbPolygon)

	# Add the fields we're interested in
	layer.CreateField(ogr.FieldDefn("Drawdown", ogr.OFTReal))
	layer.CreateField(ogr.FieldDefn("Distance", ogr.OFTReal))
	
	# create the feature
	feature = ogr.Feature(layer.GetLayerDefn())
	# Set attribute names
	feature.SetField("Drawdown", DD)
	feature.SetField("Distance", bufferDistance)
	
	# Define the point that buffer should surround
	wkt = "POINT(%f %f)" % (Wells[Well][0], Wells[Well][1])
	
	# Make that point a geometry
	point = ogr.CreateGeometryFromWkt(wkt)

	# Then create the 1st buffer around it
	poly = point.Buffer(bufferDistance)
	poly.ExportToWkt()
	pkt = poly.ExportToWkt()
	# Make it a geometry
	POLYGON = ogr.CreateGeometryFromWkt(pkt)

	# Set the feature geometry using the polygon
	feature.SetGeometry(POLYGON)
	# Create the feature in the layer (shapefile)
	layer.CreateFeature(feature)
	# Dereference the feature
	feature = None

	# Save and close the data source
	outFile = None

