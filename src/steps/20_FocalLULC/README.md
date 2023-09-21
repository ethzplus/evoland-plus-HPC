# Future-EI focal statistics tool

A Python script used to produce focal layers from the simulated LULC maps
produced for each simulation time point. Based on the tool produced for the
[SWECO25 data repository](https://github.com/NKulling/SWECO25/tree/main/focal_statistics_toolbox).

It computes focal statistics for each raster layer within the input folder. It
generates a set of six raster layers, each with focal statistics calculated using
radii ranging from 25m to 3000m.

The toolbox requires the following inputs:

1) A folder containing one or multiple raster layers.
2) An output folder to store the results.

The **mean** value is calculated within circular areas with varying radii.
Radii, focal window type and function are configurable in the global
[`config.json`](config.json).