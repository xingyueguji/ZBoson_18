#!/bin/bash
numSegments=10
opt=1
numberofsamples=10


shiftlowbin=-0.15
shifthighbin=-0.1
smearlowbin=0.007
smearhighbin=0.014

# Define the directories and macros
dir1="/Users/zhenghuang/ZBoson_18/"

macro1="cd $dir1; root -b -q './get_all_bk_mc.C(0, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro2="cd $dir1; root -b -q './get_all_bk_mc.C(1, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro3="cd $dir1; root -b -q './get_all_bk_mc.C(2, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro4="cd $dir1; root -b -q './get_all_bk_mc.C(3, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro5="cd $dir1; root -b -q './get_all_bk_mc.C(4, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro6="cd $dir1; root -b -q './get_all_bk_mc.C(5, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro7="cd $dir1; root -b -q './get_all_bk_mc.C(6, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro8="cd $dir1; root -b -q './get_all_bk_mc.C(7, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro9="cd $dir1; root -b -q './get_all_bk_mc.C(8, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"
macro10="cd $dir1; root -b -q './get_all_bk_mc.C(9, $numSegments, $opt, $numberofsamples, $shiftlowbin, $shifthighbin, $smearlowbin, $smearhighbin)'"

# Use osascript to open two new Terminal windows, change directories, and run the scripts
osascript <<EOF
tell application "Terminal"
    do script "$macro1"
    do script "$macro2"
    do script "$macro3"
    do script "$macro4"
    do script "$macro5"
    do script "$macro6"
    do script "$macro7"
    do script "$macro8"
    do script "$macro9"
    do script "$macro10"
end tell
EOF

# Exit the script to avoid creating an extra Terminal window
exit 0