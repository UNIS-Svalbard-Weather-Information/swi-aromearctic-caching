# SPDX-FileCopyrightText: 2025 Louis Pauchet louis.pauchet@insa-rouen.fr
#
# SPDX-License-Identifier: EUPL-1.2
#
# This file is part of Svalbard Weather Information - Arome Arctic Caching.
#
# Svalbard Weather Information - Arome Arctic Caching is free software:
# you can redistribute it and/or modify it under the terms of the
# European Union Public License (EUPL), version 1.2 or later.
#
# A significant portion of this code is adapted from:
#   https://github.com/UNISvalbard/unisacsi
#   Copyright (c) 2022 lfrankunis
#   Licensed under the MIT License (see below for full text).
#
# MIT License:
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import xarray as xr
import numpy as np
import os
import rioxarray
import pyproj


def download_MET_model_static_fields(
    model="AA", resolution="2p5km", out_path="./statics/"
):
    """
    Download and save the static fields (longitude, latitude, orography, and land-sea mask)
    from a model of the Norwegian Meteorological Institute (MET Norway).
    This function is adapted from the UNISACSI toolbox developed by Lukas Frank (lfrankunis) and PSelleunis.

    Parameters
    ----------
    model : str, optional
        Name of the model (e.g., "AA" for AROME-Arctic). Default is "AA".
    resolution : str, optional
        Spatial resolution of the model (e.g., "2p5km"). Default is "2p5km".
    out_path : str, optional
        Directory where the static fields NetCDF file will be saved. Default is "./statics/".

    Returns
    -------
    str
        Path to the saved NetCDF file containing the static fields.

    Notes
    -----
    The static fields include: x, y, longitude, latitude, projection_lambert,
    surface_geopotential, and land_area_fraction.
    The function creates the output directory if it does not exist.
    """
    os.makedirs(out_path, exist_ok=True)
    out_file = os.path.join(out_path, f"{model}_static_fields_{resolution}.nc")

    file = "https://thredds.met.no/thredds/dodsC/aromearcticarchive/2022/06/03/arome_arctic_det_2_5km_20220603T00Z.nc"

    with xr.open_dataset(file) as static_fields:
        selected_fields = static_fields.isel(time=0)[
            [
                "x",
                "y",
                "longitude",
                "latitude",
                "projection_lambert",
                "surface_geopotential",
                "land_area_fraction",
            ]
        ].squeeze()

        # Save to NetCDF
        selected_fields.to_netcdf(out_file, unlimited_dims=())

    print(
        f"Static fields for the model {model} resolution {resolution} were successfully saved into {out_file}."
    )
    return os.path.abspath(out_file)


def download_latest_near_surface_field_data(
    filename="https://thredds.met.no/thredds/dodsC/aromearcticlatest/archive/arome_arctic_det_2_5km_latest.nc",
    model="AA",
    resolution="2p5km",
    static_path="./statics/",
    lon_lims=[9, 33],
    lat_lims=[74, 81],
    int_x=1,
    int_y=1,
    int_h=1,
    start_h=0,
    time_steps=30,
    output="./temp/latest.nc",
):
    """
    Download and process 2D fields of near-surface meteorological variables from a specified model dataset.
    The function subsets the data spatially and temporally, reprojects it to a geographic coordinate system,
    and performs unit conversions and derived calculations (e.g., wind direction, potential temperature).

    This function is adapted from the UNISACSI toolbox developed by Lukas Frank (lfrankunis) and PSelleunis.

    Parameters
    ----------
    filename : str, optional
        URL or path to the input NetCDF file containing the model data.
        Default: "https://thredds.met.no/thredds/dodsC/aromearcticlatest/archive/arome_arctic_det_2_5km_latest.nc"
    model : str, optional
        Model identifier (e.g., "AA" for AROME-Arctic).
        Default: "AA"
    resolution : str, optional
        Model resolution (e.g., "2p5km" for 2.5 km).
        Default: "2p5km"
    static_path : str, optional
        Path to the directory containing static model fields.
        Default: "./statics/"
    lon_lims : list of float, optional
        Longitude limits for spatial subsetting: [min, max].
        Default: [9, 33]
    lat_lims : list of float, optional
        Latitude limits for spatial subsetting: [min, max].
        Default: [74, 81]
    int_x : int, optional
        Horizontal grid interval in the x-direction.
        Default: 1
    int_y : int, optional
        Horizontal grid interval in the y-direction.
        Default: 1
    int_h : int, optional
        Temporal interval for flux calculations (hours).
        Default: 1
    start_h : int, optional
        Initial time step for flux calculations.
        Default: 0
    time_steps : int, optional
        Number of time steps to extract from the dataset.
        Default: 30
    output : str, optional
        Path to the output NetCDF file.
        Default: "./temp/latest.nc"

    Returns
    -------
    xarray.Dataset
        A reprojected and processed xarray Dataset containing the near-surface fields.
        The dataset is also saved to the specified output path.

    Notes
    -----
    - The function downloads static model fields if not already present.
    - It subsets the data spatially using `lon_lims` and `lat_lims`.
    - For flux variables, it converts accumulated values to instantaneous rates.
    - Wind components are rotated from grid- to earth-related coordinates.
    - Potential temperature and relative humidity are calculated if the required variables are available.
    - The output is reprojected to WGS84 (EPSG:4326) and saved as a NetCDF file.

    """
    static_file = os.path.join(static_path, f"{model}_static_fields_{resolution}.nc")
    if not os.path.isfile(static_file):
        static_file = download_MET_model_static_fields(
            model=model, resolution=resolution, out_path=static_path
        )

    static_fields = xr.open_dataset(static_file)

    model_varnames = {
        "T": "air_temperature_2m",
        "RH": "relative_humidity_2m",
        "q": "specific_humidity_2m",
        "ur": "x_wind_10m",
        "vr": "y_wind_10m",
        "p_surf": "surface_air_pressure",
        "T_surf": "air_temperature_0m",
        # "MSLP": "air_pressure_at_sea_level",
        # "tau_x": "downward_eastward_momentum_flux_in_air",
        # "tau_y": "downward_northward_momentum_flux_in_air",
        # "SW_net": "integral_of_surface_net_downward_shortwave_flux_wrt_time",
        # "SW_down": "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time",
        # "LW_net": "integral_of_surface_net_downward_longwave_flux_wrt_time",
        # "LW_down": "integral_of_surface_downwelling_longwave_flux_in_air_wrt_time",
        # "LHF": "integral_of_surface_downward_latent_heat_flux_wrt_time",
        # "SHF": "integral_of_surface_downward_sensible_heat_flux_wrt_time",
        # "cloud_cover": "cloud_area_fraction",
        # "ABL_height": "atmosphere_boundary_layer_thickness",
        "precip": "precipitation_amount_acc",
    }

    model_varis = list(model_varnames.values()) + [
        "x",
        "y",
        "longitude",
        "latitude",
        "projection_lambert",
        "surface_geopotential",
        "land_area_fraction",
    ]
    source_crs = static_fields.projection_lambert.attrs["proj4"]

    with xr.open_dataset(filename) as full_file:
        ds = (
            full_file.isel(time=slice(0, time_steps))[model_varis]
            .where(
                (full_file.latitude >= lat_lims[0])
                & (full_file.latitude <= lat_lims[1])
                & (full_file.longitude >= lon_lims[0])
                & (full_file.longitude <= lon_lims[1]),
                drop=True,
            )
            .squeeze()
        )

    for vari in model_varis:
        if vari[:8] == "integral":
            ds[vari][dict(time=range(1, len(ds["time"])))] -= ds[vari][
                dict(time=range(0, len(ds["time"]) - 1))
            ].values
            ds[vari][dict(time=0)] /= start_h
            ds[vari] /= 3600.0 * int_h
            ds[vari].attrs["standard_name"] = ds[vari].attrs["standard_name"][12:-9]
            ds[vari].attrs["units"] = "W/m^2"
            ds[vari].attrs["long_name"] = ds[vari].attrs["long_name"][12:]
            ds = ds.rename({vari: vari[12:-9]})
        elif vari[-3:] == "acc":
            ds[vari][dict(time=range(1, len(ds["time"])))] -= ds[vari][
                dict(time=range(0, len(ds["time"]) - 1))
            ].values
            ds[vari][dict(time=0)] /= start_h
            ds[vari].attrs["long_name"] = ds[vari].attrs["long_name"][12:]
            ds = ds.rename({vari: vari[:-4]})

    print(f"Done downloading {filename}.")

    if "integral_of_surface_net_downward_shortwave_flux_wrt_time" in model_varis:
        # Converting radiative fluxes from net into up
        ds["surface_upwelling_shortwave_flux_in_air"] = (
            ds["surface_downwelling_shortwave_flux_in_air"]
            - ds["surface_net_downward_shortwave_flux"]
        )
        ds["surface_upwelling_longwave_flux_in_air"] = (
            ds["surface_downwelling_longwave_flux_in_air"]
            - ds["surface_net_downward_longwave_flux"]
        )

        del ds["surface_net_downward_shortwave_flux"]
        del ds["surface_net_downward_longwave_flux"]

    if "x_wind_10m" in model_varis:
        # Wind u and v components in the original ds are grid-related.
        # Therefore, we rotate here the wind components from grid- to earth-related coordinates.
        cone = np.sin(
            np.abs(
                np.deg2rad(ds.projection_lambert.attrs["latitude_of_projection_origin"])
            )
        )  # cone factor

        ds["diffn"] = (
            ds.projection_lambert.attrs["longitude_of_central_meridian"] - ds.longitude
        )
        ds["diffn"].values[ds["diffn"].values > 180.0] -= 360
        ds["diffn"].values[ds["diffn"].values < -180.0] += 360

        ds["alpha"] = np.deg2rad(ds["diffn"]) * cone

        ds["eastward_wind_10m"] = ds["x_wind_10m"] * np.cos(ds["alpha"]) - ds[
            "y_wind_10m"
        ] * np.sin(ds["alpha"])
        ds["northward_wind_10m"] = ds["y_wind_10m"] * np.cos(ds["alpha"]) + ds[
            "x_wind_10m"
        ] * np.sin(ds["alpha"])

        # Calculating wind direction
        ds["wind_direction_10m"] = (
            np.rad2deg(np.arctan2(-ds["eastward_wind_10m"], -ds["northward_wind_10m"]))
            + 360.0
        ) % 360.0

        # Calculating wind speed
        ds["wind_speed_10m"] = np.sqrt(
            (ds["eastward_wind_10m"] ** 2.0) + (ds["northward_wind_10m"] ** 2.0)
        )

        del ds["x_wind_10m"]
        del ds["y_wind_10m"]
        del ds["diffn"]
        del ds["alpha"]

    # Calculating potential temperature
    if ("air_temperature_2m" in model_varis) & ("surface_air_pressure" in model_varis):
        ds["air_potential_temperature_2m"] = (ds["air_temperature_2m"]) * (
            (1.0e5 / (ds["surface_air_pressure"])) ** (287.0 / 1005.0)
        )

        if "specific_humidity_2m" in model_varis:
            T = ds["air_temperature_2m"] - 273.15
            e = (ds["specific_humidity_2m"] * ds["surface_air_pressure"]) / (
                0.622 + 0.378 * ds["specific_humidity_2m"]
            )
            ds["relative_humidty_2m"] = e / (611.2 * np.exp((17.62 * T) / (243.12 + T)))

    target_crs = "+proj=longlat +datum=WGS84 +no_defs"
    ds.rio.write_crs(source_crs, inplace=True)
    ds_reprojected = ds.rio.reproject(target_crs)
    ds_reprojected = ds_reprojected.rename({"x": "longitude", "y": "latitude"})

    ds_reprojected.sel(
        latitude=slice(lat_lims[1], lat_lims[0]),
        longitude=slice(lon_lims[0], lon_lims[1]),
    ).to_netcdf(output)

    return ds_reprojected
