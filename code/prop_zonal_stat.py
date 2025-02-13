## Calculates zonal statistics with ArcPy extension

import arcpy
from arcpy import env
from arcpy.sa import *
import time

env.workspace = "m:/Users/uqfcho/Documents/remp_models_2023" # Set the workspace to where the REMP model outputs were located
outpath = "output/"
arcpy.CheckOutExtension("Spatial")

climateProjList = ["Avg", "CCCMA_R1", "CCCMA_R2", "CCCMA_R3", "CSIRO_R1", "CSIRO_R2", "CSIRO_R3", "ECHAM_R1", "ECHAM_R2", "ECHAM_R3", "MIROC_R1", "MIROC_R2", "MIROC_R3"]
timestepList = ["t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7"]

#khabThres = Raster("data/koala_habitat/khab_thresh.tif")
prop = "data/Properties_NewPropID.shp"

# Loop over all climateProj and timesteps
# Because ArcGIS should already do multicore processing, not sure if it would improve performance if this is parallelised
# Also there are overhead to distributing data files prop and khabMask to all cores

for climateProj in climateProjList:
    for timestep in timestepList:
        if climateProj == "Avg":
            in_path = f"data/Combine_Inland_Coastal/Koala_{climateProj}_Stable_Coupled_Combined_Pi_{timestep}.tif"
        else:
            in_path = f"data/Combine_Inland_Coastal/{climateProj}/Occupancy/Koala_{climateProj}_Stable_Coupled_Combined_Pi_{timestep}.tif"

        start_time = time.time()
        # Load raster file into memory
        modelRaster = Raster(in_path)

        # Multiply climate data with mask, cells with no data in khabMask will return NODATA, thus will not be factored in zonal statistics
        modelKhabMasked = modelRaster # * khabMask

        # Conduct zonal statistics operation, saved as dbf
        ZonalStatisticsAsTable(prop, "NewPropID", modelKhabMasked, f"output/kitl_prop_unmasked.gdb/kitl_prop_{climateProj}_{timestep}", "NODATA", "MEAN", "CURRENT_SLICE")

        # Monitor progress
        print(f"Completed: {climateProj} {timestep} in {time.time() - start_time} seconds")