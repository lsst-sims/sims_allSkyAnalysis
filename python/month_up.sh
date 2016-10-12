#!/bin/bash

# Copy the month of cloud data up to 

rsync -av --progress month_o_clouds.npz lsst-dev.ncsa.illinois.edu:"/lsst/sim-data/CloudMaps/."
