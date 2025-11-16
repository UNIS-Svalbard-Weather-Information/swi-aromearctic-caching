# SWI Arome Arctic Caching

## Overview
This service prepares Arome Arctic data for use in the frontend.

## Folder Structure
The `/swi/data` directory should be mapped to the generic observation/weather volume. The expected structure is as follows:

| Path | Description |
|------|-------------|
| `/swi/data/forecast/aa/velocity` | Contains `grib2json`-ready files used by **Leaflet Velocity**. |
| `/swi/data/forecast/aa/cog` | Contains GeoTIFF files for static field variables. |

## Usage
To run the service, use the following Docker command:
```bash
docker run -v ./data:/swi/data ghcr.io/unis-svalbard-weather-information/swi-aromearctic-caching:latest
```