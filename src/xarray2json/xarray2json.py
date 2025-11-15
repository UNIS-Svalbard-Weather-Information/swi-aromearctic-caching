import xarray as xr
from typing import Union
import numpy as np
from scipy.interpolate import griddata
import json
import pandas as pd
from pyproj import CRS, Transformer


class Xarray2Json:
    def __init__(
        self,
        path: str,
        bounding_box: dict = {
            "lat_min": 74,
            "lat_max": 81,
            "lon_min": 9,
            "lon_max": 33,
        },
        safety_factor: float = 0.8,
        processing_period: int = 25,
        model_resolution: float = 2.5,
        lat_var_name: str = "latitude",
        lon_var_name: str = "longitude",
    ):
        self.ds = xr.open_dataset(path)
        self.bounding_box = bounding_box
        self.safety_factor = safety_factor
        self.processing_period = processing_period
        self.model_resolution = model_resolution
        self.lat_var_name = lat_var_name
        self.lon_var_name = lon_var_name
        self.reprojected_data = {}

        # Calculate latitude step
        self.lat_step = (model_resolution / 1.852) / 60 * self.safety_factor
        self.target_lats = np.arange(
            self.bounding_box["lat_min"],
            self.bounding_box["lat_max"]
            + self.lat_step,  # Ensure upper bound is included
            self.lat_step,
        )

        # Calculate longitude step at the minimum latitude (most restrictive)
        self.lon_step = (
            model_resolution
            / (111.320 * np.cos(np.deg2rad(self.bounding_box["lat_min"])))
        ) * self.safety_factor
        self.target_lons = np.arange(
            self.bounding_box["lon_min"],
            self.bounding_box["lon_max"]
            + self.lon_step,  # Ensure upper bound is included
            self.lon_step,
        )

    def reproject_to_target_grid(
        self, var_name: str, sel=None, isel=None, method="linear"
    ) -> xr.DataArray:
        if sel is None:
            sel = {}
        if isel is None:
            isel = {}

        cache_key = (var_name, tuple(sorted(sel.items())), tuple(sorted(isel.items())))
        if cache_key in self.reprojected_data:
            return self.reprojected_data[cache_key]

        var = self.ds[var_name].sel(**sel).isel(**isel)

        # Extract source lats and lons (assuming they are 2D)
        source_lats = var[self.lat_var_name].values
        source_lons = var[self.lon_var_name].values

        # Flatten the source grids and variable values
        points = np.column_stack((source_lats.ravel(), source_lons.ravel()))
        values = var.values.ravel()

        # Create meshgrid for target points
        grid_target_lon, grid_target_lat = np.meshgrid(
            self.target_lons, self.target_lats
        )
        target_points = np.column_stack(
            (grid_target_lat.ravel(), grid_target_lon.ravel())
        )

        # Perform interpolation
        interpolated_values = griddata(points, values, target_points, method=method)

        # Reshape to target grid shape
        interpolated_var = interpolated_values.reshape(grid_target_lon.shape)

        interpolated_da = xr.DataArray(
            interpolated_var,
            dims=["lat", "lon"],
            coords={
                "lat": self.target_lats,
                "lon": self.target_lons,
                "time": var.time,  # Add time coordinate
            },
            name=var.name,
        )
        interpolated_da.attrs = var.attrs

        self.reprojected_data[cache_key] = interpolated_da
        return interpolated_da

    def generate_json_leaflet_velocity(
        self,
        variable_u: str,
        variable_v: str,
        output_path: str = None,
        time: str = None,
        reproject: bool = True,
        crs_source: str = CRS.from_proj4(
            "+proj=lcc +lat_1=77.5 +lat_2=77.5 +lat_0=77.5 +lon_0=-25 +R=6371000 +units=m +no_defs"
        ),
        crs_destination: str = CRS.from_epsg(4326),
    ):
        if not output_path:
            output_path = f"leaflet_velocity_{variable_u}_{variable_v}_{pd.to_datetime(time).isoformat().replace(':', '') + 'Z'}.json"

        u = self.reproject_to_target_grid(variable_u, sel={"time": time})
        v = self.reproject_to_target_grid(variable_v, sel={"time": time})

        # if reproject:
        #     transformer = Transformer.from_crs(lambert_crs, wgs84_crs, always_xy=True)

        #     u =

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
            "refTime": pd.to_datetime(time).isoformat() + "Z",
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
                    "nx": u.lon.size,
                    "ny": u.lat.size,
                    "lo1": u.lon.values.min(),
                    "la1": u.lat.values.max(),
                    "lo2": u.lon.values.max(),
                    "la2": u.lat.values.min(),
                    "dx": u.lon.values[1] - u.lon.values[0],
                    "dy": u.lat.values[1] - u.lat.values[0],
                },
                "data": np.round(u.values.flatten(), decimals=1).tolist(),
            },
            {
                "header": {
                    **header,
                    "parameterNumber": 3,
                    "parameterNumberName": variable_v,
                    "nx": v.lon.size,
                    "ny": v.lat.size,
                    "lo1": v.lon.values.min(),
                    "la1": v.lat.values.max(),
                    "lo2": v.lon.values.max(),
                    "la2": v.lat.values.min(),
                    "dx": v.lon.values[1] - v.lon.values[0],
                    "dy": v.lat.values[1] - v.lat.values[0],
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
