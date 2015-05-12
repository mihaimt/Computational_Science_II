#!/bin/bash
# execute in main folder with ./plot.sh

if [ -x grav_refined ]; then
    ./grav_refined 1
    ./grav_refined 2
    ./grav_refined 3
    ./grav_refined 4
    ./grav_refined 5
fi 

if [ -x grav_pure_lvlx ]; then
    ./grav_pure_lvlx 1
    ./grav_pure_lvlx 2
    ./grav_pure_lvlx 3
    ./grav_pure_lvlx 4
fi

if [ -x grav_force_lvlx ]; then
    ./grav_force_lvlx 1
    ./grav_force_lvlx 2
    ./grav_force_lvlx 3
    ./grav_force_lvlx 4
fi

if [ -x grav_force ]; then
    ./grav_force
fi

if [ -x oscillation_force ]; then
    ./oscillation_force 1
    ./oscillation_force 2
    ./oscillation_force 3
    ./oscillation_force 4
fi

if [ -x grav_force_lvl1 ]; then
    ./grav_force_lvl1
fi

if [ -x grav_force_lvl2 ]; then
    ./grav_force_lvl2
fi

cd bin/
./force_plot.py
cd ../
