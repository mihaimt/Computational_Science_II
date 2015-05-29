#!/bin/bash

./diff_gridlvl ../data/f_radial.data ../data/radial_pure_lvl1.data ../data/diff/radial_pure_diff1.data
./diff_gridlvl ../data/f_radial.data ../data/radial_pure_lvl2.data ../data/diff/radial_pure_diff2.data
./diff_gridlvl ../data/f_radial.data ../data/radial_pure_lvl3.data ../data/diff/radial_pure_diff3.data
./diff_gridlvl ../data/f_radial.data ../data/radial_pure_lvl4.data ../data/diff/radial_pure_diff4.data

./diff_gridlvl ../data/f_angular.data ../data/angular_pure_lvl1.data ../data/diff/angular_pure_diff1.data
./diff_gridlvl ../data/f_angular.data ../data/angular_pure_lvl2.data ../data/diff/angular_pure_diff2.data
./diff_gridlvl ../data/f_angular.data ../data/angular_pure_lvl3.data ../data/diff/angular_pure_diff3.data
./diff_gridlvl ../data/f_angular.data ../data/angular_pure_lvl4.data ../data/diff/angular_pure_diff4.data


./diff_gridlvl ../data/f_radial.data ../data/radial_refined_lvl1.data ../data/diff/radial_refined_diff1.data
./diff_gridlvl ../data/f_radial.data ../data/radial_refined_lvl2.data ../data/diff/radial_refined_diff2.data
./diff_gridlvl ../data/f_radial.data ../data/radial_refined_lvl3.data ../data/diff/radial_refined_diff3.data
./diff_gridlvl ../data/f_radial.data ../data/radial_refined_lvl4.data ../data/diff/radial_refined_diff4.data

./diff_gridlvl ../data/f_angular.data ../data/angular_refined_lvl1.data ../data/diff/angular_refined_diff1.data
./diff_gridlvl ../data/f_angular.data ../data/angular_refined_lvl2.data ../data/diff/angular_refined_diff2.data
./diff_gridlvl ../data/f_angular.data ../data/angular_refined_lvl3.data ../data/diff/angular_refined_diff3.data
./diff_gridlvl ../data/f_angular.data ../data/angular_refined_lvl4.data ../data/diff/angular_refined_diff4.data

./residual_plot.py
