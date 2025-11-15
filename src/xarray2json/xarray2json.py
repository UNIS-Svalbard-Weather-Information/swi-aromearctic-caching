import xarray as xr
from typing import Union
import numpy as np
import json
import pandas as pd


class Xarray2Json:
    def __init__(
        self,
        path: str,
        processing_period: int = 25,
    ):
        self.ds = xr.open_dataset(path)
        # self.bounding_box = bounding_box
        # self.safety_factor = safety_factor
        self.processing_period = processing_period
        # self.model_resolution = model_resolution
        # self.lat_var_name = lat_var_name
        # self.lon_var_name = lon_var_name
        # self.reprojected_data = {}

        # # Calculate latitude step
        # self.lat_step = (model_resolution / 1.852) / 60 * self.safety_factor
        # self.target_lats = np.arange(
        #     self.bounding_box["lat_min"],
        #     self.bounding_box["lat_max"]
        #     + self.lat_step,  # Ensure upper bound is included
        #     self.lat_step,
        # )

        # # Calculate longitude step at the minimum latitude (most restrictive)
        # self.lon_step = (
        #     model_resolution
        #     / (111.320 * np.cos(np.deg2rad(self.bounding_box["lat_min"])))
        # ) * self.safety_factor
        # self.target_lons = np.arange(
        #     self.bounding_box["lon_min"],
        #     self.bounding_box["lon_max"]
        #     + self.lon_step,  # Ensure upper bound is included
        #     self.lon_step,
        # )

    def generate_json_leaflet_velocity(
        self,
        variable_u: str,
        variable_v: str,
        output_path: str = None,
        time: Union[str, int] = 0,
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

    def generate_json_grided(self, variable: str, output_path: str, time: str = None):
        pass

    def run_configuration(self, config: Union[dict, str]):
        if isinstance(config, str):
            import json

            with open(config, "r") as f:
                config = json.load(f)
        pass
