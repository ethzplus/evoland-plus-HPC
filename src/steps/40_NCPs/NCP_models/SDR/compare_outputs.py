"""Script to compare two outputs of the SDR mode with using the non filled
and the prefilled DEM

We got the following outputs:
sed_deposition.tif  sed_export.tif  sed_retention_index.tif  sed_retention.tif
watershed_results_sdr.dbf  watershed_results_sdr.prj  watershed_results_sdr.shp  watershed_results_sdr.shx

One time in
ncp_output/1/SDR/2020/compare/
and one time in
ncp_output/1/SDR/2020/

We want to find out if the output files have the same results
"""

import os
from os.path import join, basename

import pygeoprocessing as pgp
from osgeo import gdal

output_dir_1 = "ncp_output/1/SDR/2020/compare/"
output_dir_2 = "ncp_output/1/SDR/2020/"

output_files = [
    "sed_deposition.tif",
    "sed_export.tif",
    "sed_retention_index.tif",
    "sed_retention.tif",
    "watershed_results_sdr.shp",
]

def compare_tifs(file1, file2):
    print(f"Comparing {basename(file1)}")
    # compare the files
    os.system(f"gdalinfo {file1} > {file1}.txt")
    os.system(f"gdalinfo {file2} > {file2}.txt")
    os.system(f"diff {file1}.txt {file2}.txt > {file1}_diff.txt")
    with open(f"{file1}_diff.txt", "r") as f:
        diff = f.read()
    if diff:
        print("Metadata of the files are different:")
        print(diff)
    else:
        print("Metadata of the files are the same")

    # load files to arrays and compare them
    array1 = pgp.raster_to_numpy_array(file1)
    array2 = pgp.raster_to_numpy_array(file2)
    if not (array1 == array2).all():
        print("Content of the files are different\n")
    else:
        print("Content of the files are the same\n")

def compare_shp(file1, file2):
    print(f"Comparing {basename(file1)}")
    # compare the shapefiles
    os.system(f"ogrinfo {file1} > {file1}.txt")
    os.system(f"ogrinfo {file2} > {file2}.txt")
    os.system(f"diff {file1}.txt {file2}.txt > {file1}_diff.txt")
    with open(f"{file1}_diff.txt", "r") as f:
        diff = f.read()
    if diff:
        print("Metadata of the files are different:")
        print(diff)
    else:
        print("Metadata of the files are the same")

    # compare the content of the shapefiles
    shp1 = gdal.OpenEx(file1, gdal.OF_VECTOR)
    shp2 = gdal.OpenEx(file2, gdal.OF_VECTOR)
    # print number of layers
    print(f"Amount of layers: {shp1.GetLayerCount()} vs {shp2.GetLayerCount()}")
    print("Comparing the first layer of the shapefiles")
    layer1 = shp1.GetLayer()
    layer2 = shp2.GetLayer()
    if layer1.GetFeatureCount() != layer2.GetFeatureCount():
        print(f"Number of features are different: {layer1.GetFeatureCount()} vs {layer2.GetFeatureCount()}")
        exit()
    else:
        print(f"Number of features are the same: {layer1.GetFeatureCount()}")

    # As order might differ, convert to WKT (Well Known Text) and sort
    wkt1 = sorted([feature.GetGeometryRef().ExportToWkt() for feature in layer1])
    wkt2 = sorted([feature.GetGeometryRef().ExportToWkt() for feature in layer2])

    if wkt1 != wkt2:
        print("Content of the files are different\n")

    print("Content of the files are the same\n")



if __name__ == "__main__":
    tif_files = [f for f in os.listdir(output_dir_1) if f.endswith(".tif")]
    for output_file in tif_files:
        compare_tifs(join(output_dir_1, output_file), join(output_dir_2, output_file))

    shp_files = [f for f in os.listdir(output_dir_1) if f.endswith(".shp")]
    for output_file in shp_files:
        compare_shp(join(output_dir_1, output_file), join(output_dir_2, output_file))

