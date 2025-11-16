# SPDX-FileCopyrightText: 2025 Louis Pauchet louis.pauchet@insa-rouen.fr
#
# SPDX-License-Identifier: EUPL-1.2
#
# This file is part of Svalbard Weather Information - Arome Arctic Caching.
#
# Svalbard Weather Information - Arome Arctic Caching is free software:
# you can redistribute it and/or modify it under the terms of the
# European Union Public License (EUPL), version 1.2 or later.

from loguru import logger
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
    # Configure loguru to write logs to a file and console
    logger.add("app.log", rotation="10 MB", level="INFO")
    logger.info("Starting the main process")

    try:
        os.makedirs(PATH_VELOCITY, exist_ok=True)
        os.makedirs(PATH_COG, exist_ok=True)
        logger.info(f"Created directories: {PATH_VELOCITY}, {PATH_COG}")

        logger.info("Downloading latest near surface field data...")
        nc_file = download_latest_near_surface_field_data(time_steps=HOURS_TO_PROCESS)
        logger.info(f"Downloaded data to: {nc_file}")

        xarray2json_obj = xarray2json.Xarray2Json(path=nc_file)
        logger.info("Initialized Xarray2Json object")

        for h in range(0, HOURS_TO_PROCESS):
            logger.info(f"Processing hour {h}/{HOURS_TO_PROCESS}")

            for p in VELOCITY_VARIABLES:
                logger.info(f"Generating velocity JSON for variables: {p}")
                xarray2json_obj.generate_json_leaflet_velocity(
                    **p, time=h, output_path=PATH_VELOCITY
                )

            for v in COG_VARIABLES:
                logger.info(f"Generating COG GeoJSON for variable: {v['variable']}")
                xarray2json_obj.generate_geojson_grided(time=h, **v)

        logger.success("Processing completed successfully")

    except Exception as e:
        logger.exception(f"An error occurred: {e}")
        raise


if __name__ == "__main__":
    main()
