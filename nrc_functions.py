################### FUNCTION DEFINITIONS #################################################
## BRIGHTEN PIXEL#################################
def brighten(band):
    alpha=0.13
    beta=0
    
    return np.clip(alpha*band+beta, 0,255)
###################################################

## NORMALISE ######################################
def normaliser(band):
    band_min =np.nanmin(band)
    band_max = np.nanmax(band)
    nband = (band-band_min)/(band_max - band_min)
       
    return nband
####################################################

##  WRITE IMAGE ####################################################
def writeImage (imageData, x, y, strFileName) :
    
    transform = rasterio.Affine.translation(x, y)* Affine.scale(10.0, -10.0)
    with rasterio.open(os.path.join('c:\Temp\dhk\Sentinel55', strFileName+'.tif'), 
                       'w', 
                       driver='GTiff',
                       height = imageData.shape[0], 
                       width = imageData.shape[1],
                       count=1, dtype=str(imageData.dtype),
                       #count=1,
                       #dtype='int16',
                       crs='epsg:2193',
                      transform=transform) as dst:
                      dst.write(imageData.astype(rasterio.float32),1)
#####################################################################

##  WRITE MASK IMAGE ####################################################
def writeMaskImage (imageData, x, y, strFileName) :
    
    transform = rasterio.Affine.translation(x, y)* Affine.scale(10.0, -10.0)
    with rasterio.open(os.path.join('c:\Temp\dhk\Sentinel', strFileName+'.tif'), 
                       'w', 
                       driver='GTiff',
                       height = imageData.shape[0], 
                       width = imageData.shape[1],
                       #count=1, dtype=str(imageData.dtype),
                       count=1,
                       dtype='int16',
                       crs='epsg:2193',
                      transform=transform) as dst:
                      dst.write(imageData.astype(rasterio.float32),1)
#####################################################################

## WRITE RGB IMAGE ##################################################

def writeRGB (red,green,blue, x, y, strFileName):    

    #PRE PROCESS
    red_b=brighten(red)
    blue_b=brighten(blue)
    green_b=brighten(green)
    
    red_bn = normaliser(red_b)
    green_bn = normaliser(green_b)
    blue_bn = normaliser(blue_b)

    rgb_composite_bn= np.dstack((red_bn, green_bn, blue_bn))
    rgb_composite_bn2 =  np.rollaxis(rgb_composite_bn, axis=2)  # Roll axis 2 to 0th position

    transform = rasterio.Affine.translation(x, y)* Affine.scale(10.0, -10.0)
    with rasterio.open(os.path.join('c:\Temp\dhk\Sentinel', strFileName+'.tif'), 
                       'w',
                        driver='GTiff',
                        height=rgb_composite_bn2.shape[1],
                        width=rgb_composite_bn2.shape[2],
                        count=3,
                        dtype=rgb_composite_bn2.dtype,
                        nodata=0,
                        crs='epsg:2193',
                        transform=transform) as dst:
                        dst.write(rgb_composite_bn2)


#####################################################################

########################################################################
def writeFeatures ( sLayer, sID, sDate, sClass, polygons ):
    
    targetFields = ['SHAPE@', 'Sentinel_ID','Date', 'SCL_Class']
    rowcount = 0

    print ("Layer : ",sLayer)
    print (sID, sDate, sClass)
    
    with arcpy.da.InsertCursor(sLayer,targetFields) as cursor:

        for polygon in polygons:
            
            rowcount = rowcount + 1

            poly_geom  = arcgis.geometry.Geometry.from_shapely(polygon)
            poly_geom.spatialReference = {'wkid': 2193}

            array = arcpy.Array([])

            for i in range(len(poly_geom['rings'][0])):

                x = poly_geom['rings'][0][i][0]
                y = poly_geom['rings'][0][i][1]

                array.append(arcpy.Point(x, y))

            poly = arcpy.Polygon(array)

            # Open an InsertCursor and insert the new geometry
            # Comment out for testing
            cursor.insertRow((poly,sID,sDate,sClass))
            
    return rowcount
 #################################################################### 
    
## MERGE FEATURES ###################################################
def DissolveGeoms(geoms):
    cnt = 0
    for geom in geoms:
        cnt += 1
        if cnt == 1:
            diss_geom = geom
        else:
            diss_geom = diss_geom.union(geom)
    return diss_geom
    

#####################################################################

## WRITE POLYGON MASK FEATURES ######################################

def writeMaskShapes ( inImage, inImageName, inLayer, intClass, strClass, imageDate, x, y):
    
    #Temporary
    inLayer = 'C:\Temp\Sentinel\Sentinel\Sentinel.gdb\Sentinel_SCL'

    # cast image to int16
    image16 = inImage.astype('int16')
    # Set mask
    mask = image16 == intClass
    # Get shapes that match mask
    shapes = features.shapes(image16, mask=mask)
    #Convert Array to polygons
    print ("coverting array result to polygons...")
    transform1 = rasterio.Affine.translation(x, y)* Affine.scale(10.0, -10.0)
    shapes = features.shapes(image16, mask=mask, transform = transform1)

    polygon_shapes = [shapely.geometry.Polygon(shape[0]["coordinates"][0]) for shape in shapes if shape[1] == intClass]
    #merged_shapes = DissolveGeoms(polygon_shapes)
    
    # write features to cloud mask
    print ("Writing features to Cloud Mask")
    maskAdded = writeFeatures ( inLayer, inImageName, imageDate, strClass, polygon_shapes )
    print (maskAdded, intClass, strClass," : Added")



## LOAD GRID FEATURES ###############################################

def load_grid(shapefile_path, target_crs):
    grid_shape = gpd.read_file(shapefile_path)
    return grid_shape.to_crs(target_crs)

## LOAD CONSENT FEATURES ###############################################

def load_consents(shapefile_path, target_crs):
    consent_shape = gpd.read_file(shapefile_path)
    return grid_shape.to_crs(target_crs)

## GET INTERSECT WITH CONSENT FEATURES ###############################################

def grid_intersect(grid_path, consents_path):
    df_grid = gpd.read_file(grid_path)
    df_consents = gpd.read_file(consents_path)
    joined_df =gpd.sjoin(df_consents, df_grid, how='left', predicate='within')
    grid_keys = list(joined_df["GRID_ID"])
    values = list(joined_df["GRID_ID"])
    grid_id = dict(zip(grid_keys, values))
    grid_id = list(grid_id.keys())
    return grid_id 

def get_grid_info(grid_shape, grid_id):
    grid_dict = dict(zip(list(grid_shape.GRID_ID), list(grid_shape.geometry)))
    polygon_str = grid_dict[grid_id].wkt
    polygon = loads(polygon_str)
    coordinates_list = list(polygon.exterior.coords)
    return coordinates_list

def get_bbox(coordinates_list):
    xmin, ymin = coordinates_list[2][0], coordinates_list[2][1]
    xmax, ymax = coordinates_list[0][0], coordinates_list[0][1]
    return [xmin, ymin, xmax, ymax]

#########################################################################################