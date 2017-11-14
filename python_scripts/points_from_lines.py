"""
Create a series of points at a specified distance along a line using shapely
FJC
24/01/17
"""

def get_points_along_line(DataDirectory, baseline_shapefile, distance, output_shapefile):
    """
    Interpolate a series of points at equal distances along an input line shapefile. Arguments that need
    to be supplied are:
    * DataDirectory: the directory of the input/output shapefiles
    * baseline_shapefile: the name of the input line shapefile with extension
    * distance: the distance to place points at
    * output_shapefile: the name of the output points shapefile with extension
    """

    from fiona import collection
    from shapely.geometry import shape, Point, LineString, mapping

    lines = []
    points = []
    distances = []
    # read in the baseline shapefile
    c = collection(DataDirectory+baseline_shapefile, 'r')
    rec = c.next()
    line = LineString(shape(rec['geometry']))
    # get the coordinate system from the input shapefile
    crs = c.crs

    total_distance = line.length
    # handle exceptions
    if distance < 0.0 or distance >= total_distance:
        print "\tNot a valid distance, sorry pal!"

    # get the points at the specified distance along the line
    temp_distance = 0
    n_points = int(total_distance/distance)
    print "The total distance is", total_distance, ": returning ", n_points, "points"
    # have a point at the start of the line
    for j in range(n_points+1):
        point = line.interpolate(temp_distance)
        points.append(Point(point))
        distances.append(temp_distance)
        temp_distance+=distance

    #output schema
    schema={'geometry': 'Point', 'properties': {'distance': 'float'} }

    # write the points to a shapefile
    with collection(DataDirectory+output_shapefile, 'w', crs=crs, driver='ESRI Shapefile', schema=schema) as output:
        for i in range (n_points+1):
            #print point
            output.write({'properties':{'distance':distances[i]},'geometry': mapping(points[i])})

if __name__ == '__main__':

    DataDirectory = '/media/fionaclubb/terrace_lidar/DEMs_for_analysis/Upper_Miss_reach11/'
    get_points_along_line(DataDirectory,baseline_shapefile='Upper_Miss_baseline_reach11.shp',distance=5,output_shapefile='Upper_Miss_reach11_baseline_points.shp')
