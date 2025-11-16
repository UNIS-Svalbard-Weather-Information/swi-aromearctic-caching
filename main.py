# SPDX-FileCopyrightText: 2025 Louis Pauchet louis.pauchet@insa-rouen.fr
#
# SPDX-License-Identifier: EUPL-1.2
#
# This file is part of Svalbard Weather Information - Arome Arctic Caching.
#
# Svalbard Weather Information - Arome Arctic Caching is free software:
# you can redistribute it and/or modify it under the terms of the
# European Union Public License (EUPL), version 1.2 or later.


from src.aa.downloader import download_latest_near_surface_field_data
from src.xarray2json import xarray2json
import os

PATH_VELOCITY = "./data/forecast/aa/velocity"
PATH_COG = "./data/forecast/aa/cog"
HOURS_TO_PROCESS = 35

VELOCITY_VARIABLES = [
    {"variable_u": "eastward_wind_10m", "variable_v": "northward_wind_10m"},
]

COG_VARIABLES = [
    {"variable": "air_temperature_2m", "conversion_fct": lambda k: k - 273.15},
    {"variable": "surface_air_pressure"},
    {"variable": "precipitation_amount"},
    {"variable": "wind_speed_10m"},
    {"variable": "wind_direction_10m"},
    {"variable": "air_temperature_0m", "conversion_fct": lambda k: k - 273.15},
    {"variable": "surface_air_pressure", "conversion_fct": lambda k: k / 100},
]


def main():
    os.makedirs(PATH_VELOCITY, exist_ok=True)
    os.makedirs(PATH_COG, exist_ok=True)

    nc_file = download_latest_near_surface_field_data(time_steps=HOURS_TO_PROCESS)

    xarray2json_obj = xarray2json.Xarray2Json(
        path=nc_file,
    )

    for h in range(0, HOURS_TO_PROCESS):
        for p in VELOCITY_VARIABLES:
            xarray2json_obj.generate_json_leaflet_velocity(
                **p, time=h, output_path=PATH_VELOCITY
            )

        for v in COG_VARIABLES:
            xarray2json_obj.generate_geojson_grided(time=h, **v)


if __name__ == "__main__":
    main()
