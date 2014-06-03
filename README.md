
GDAL based node/way density rasterizer for OSM data
===================================================

This tool is a merger between [gdal_rasterize](http://gdal.org/gdal_rasterize.html) 
and the Osmium [nodedensity](https://github.com/joto/osmium/blob/master/examples/nodedensity.cpp)
tool.  It allows generating density plots of OSM data in arbitrary projections and resolutions.

In addition to basic node density it provides a way density mode generating density data proportional
to the length of ways rather than to the number of nodes.

Building it requires additional source files from GDAL (just like gdal_rasterize), i.e. you need a 
GDAL source package.

