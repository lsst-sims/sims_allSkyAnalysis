#!/bin/sh

find MoviePlots -name "*.png" -delete
./makemovie.py
rm MoviePlots/out.mp4
ffmpeg -framerate 10 -pattern_type glob -i 'MoviePlots/*.png'  MoviePlots/out.mp4
