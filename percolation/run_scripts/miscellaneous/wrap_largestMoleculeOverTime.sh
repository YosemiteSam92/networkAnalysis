#!/bin/bash

# python3 ../src/largestMoleculeOverTime.py
# rm dummy

function processSystem() {
  data=$1
  baseConn=$2
  endConn=$3
  baseTraj=$4
  endTraj=$5
  outputLargest=$6
  outputOthers=$7
  shift 7
  arr=("$@")
  for i in ${arr[@]}
  do
    pathToConnectivityDump=$baseConn$i$endConn
    pathToTrajectoryDump=$baseTraj$i$endTraj
    python3 ../src/largestMoleculeOverTime.py -c $pathToConnectivityDump \
    -t $pathToTrajectoryDump
  done
}

# dump indices
large1=(2 4 6 8 10 12 14 16 18)

basePathToConnectivityDump="/home/mattia/Documents/emmy_backup/large_503K_1/dump/dump."
endPathToConnecvitiyDump="_btype_batom1_batom2_nptCuring_gelation"

basePathToTrajectoryDump="/home/mattia/Documents/emmy_backup/large_503K_1/dump/dump."
endPathToTrajectoryDump="_atomId_type_scaledCoords_nptCuring_gelation"

processSystem "foo" $basePathToConnectivityDump $endPathToConnecvitiyDump \
$basePathToTrajectoryDump $endPathToTrajectoryDump "foo" "foo" ${large1[@]}

rm dummy # clean up dummy log file
