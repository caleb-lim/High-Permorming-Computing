#!/bin/bash
#update this as necessary 
exe=/opt/homebrew/bin/ffmpeg 
#base name is likely GOL.grid-*-by-*.step-
basename=$1
$exe -framerate 5 -i ${basename}%d.png -c:v libx264 -r 30 -pix_fmt yuv420p GOL.mp4
