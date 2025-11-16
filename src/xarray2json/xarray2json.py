# SPDX-FileCopyrightText: 2025 Louis Pauchet louis.pauchet@insa-rouen.fr
#
# SPDX-License-Identifier: EUPL-1.2
#
# This file is part of Svalbard Weather Information - Arome Arctic Caching.
#
# Svalbard Weather Information - Arome Arctic Caching is free software:
# you can redistribute it and/or modify it under the terms of the
# European Union Public License (EUPL), version 1.2 or later.


import xarray as xr
from typing import Union
import numpy as np
import json
import pandas as pd
import rioxarray


class Xarray2Json:
    def __init__(
        self,
        path: str,
        processing_period: int = 25,
    ):
        self.ds = xr.open_dataset(path)
        self.processing_period = processing_period

    def generate_json_leaflet_velocity(
        self,
        variable_u: str,
        variable_v: str,
        output_path: str = None,
        time: int = 0,
    ):
        if not output_path:
            output_path = f"leaflet_velocity_{variable_u}_{variable_v}_{pd.to_datetime(self.ds.time.values[time]).isoformat().replace(':', '') + 'Z'}.json"

        u = self.ds[variable_u].isel(time=time)
        v = self.ds[variable_v].isel(time=time)

        header = {
            "discipline": 0,
            "disciplineName": "Meteorological products",
            "gribEdition": 2,
            "gribLength": 76420,
            "center": 7,
            "centerName": self.ds.attrs.get(
                "institution", self.ds.attrs.get("publisher_institution", "Unknown")
            ),
            "subcenter": 0,
            "refTime": pd.to_datetime(self.ds.time.values[time]).isoformat() + "Z",
            "significanceOfRT": 1,
            "significanceOfRTName": "Start of forecast",
            "productStatus": 0,
            "productStatusName": "Operational products",
            "productType": 1,
            "productTypeName": "Forecast products",
            "productDefinitionTemplate": 0,
            "productDefinitionTemplateName": "Analysis/forecast at horizontal level/layer at a point in time",
            "parameterCategory": 2,
            "parameterCategoryName": "Momentum",
            "parameterNumber": 2,
            "parameterNumberName": "U-component_of_wind",
            "parameterUnit": "m.s-1",
            "genProcessType": 2,
            "genProcessTypeName": "Forecast",
            "forecastTime": 0,
            "surface1Type": 103,
            "surface1TypeName": "Specified height level above ground",
            "surface1Value": 10.0,
            "surface2Type": 255,
            "surface2TypeName": "Missing",
            "surface2Value": 0.0,
            "gridDefinitionTemplate": 0,
            "gridDefinitionTemplateName": "Latitude_Longitude",
            "numberPoints": 65160,
            "shape": 6,
            "shapeName": "Earth spherical with radius of 6,371,229.0 m",
            "gridUnits": "degrees",
            "resolution": 48,
            "winds": "true",
            "scanMode": 0,
            "nx": 360,
            "ny": 181,
            "basicAngle": 0,
            "subDivisions": 0,
            "lo1": 0.0,
            "la1": 90.0,
            "lo2": 359.0,
            "la2": -90.0,
            "dx": 1.0,
            "dy": 1.0,
        }

        results = [
            {
                "header": {
                    **header,
                    "parameterNumber": 2,
                    "parameterNumberName": variable_u,
                    "nx": u.longitude.size,
                    "ny": u.latitude.size,
                    "lo1": u.longitude.values.min(),
                    "la1": u.latitude.values.max(),
                    "lo2": u.longitude.values.max(),
                    "la2": u.latitude.values.min(),
                    "dx": np.abs(v.longitude.values[1] - v.longitude.values[0]),
                    "dy": np.abs(v.latitude.values[1] - v.latitude.values[0]),
                },
                "data": np.round(u.values.flatten(), decimals=1).tolist(),
            },
            {
                "header": {
                    **header,
                    "parameterNumber": 3,
                    "parameterNumberName": variable_v,
                    "nx": v.longitude.size,
                    "ny": v.latitude.size,
                    "lo1": v.longitude.values.min(),
                    "la1": v.latitude.values.max(),
                    "lo2": v.longitude.values.max(),
                    "la2": v.latitude.values.min(),
                    "dx": np.abs(v.longitude.values[1] - v.longitude.values[0]),
                    "dy": np.abs(v.latitude.values[1] - v.latitude.values[0]),
                },
                "data": np.round(v.values.flatten(), decimals=1).tolist(),
            },
        ]

        with open(output_path, "w") as f:
            json.dump(results, f)

        pass

    def generate_geojson_grided(
        self, variable: str, output_path: str = None, time: int = 0, conversion_fct=None
    ):
        if not output_path:
            output_path = f"cog_{variable}_{pd.to_datetime(self.ds.time.values[time]).isoformat().replace(':', '') + 'Z'}.tif"

        da = self.ds[variable].isel(time=time)

        if conversion_fct:
            da = conversion_fct(da)

        da.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)
        da.rio.write_crs("epsg:4326", inplace=True)
        da_reprojected = da.rio.reproject("epsg:3857")

        # da_reprojected.rio.to_raster(
        #     output_path,
        #     tiled=True,
        #     windowed=True,
        #     driver="COG",  # Explicitly set the driver to COG
        #     compress="LZW",  # Or "LZW" for better compatibility
        #     overview_levels=[2, 4, 8],  # Add overviews for better performance
        #     dtype="float32",  # Ensure data type is compatible
        # )

        da_reprojected.rio.to_raster(
            output_path,
            driver="COG",
            compress="DEFLATE",
            dtype="float32",
            overview_levels=[2, 4, 8, 10],
            tiled=True,
            windowed=True,
        )

    def run_configuration(self, config: Union[dict, str]):
        if isinstance(config, str):
            import json

            with open(config, "r") as f:
                config = json.load(f)
        pass
