#!/bin/bash

function processSystem() {
  data=$1
  baseConn=$2
  endConn=$3
  baseOutput=$4
  endOutput=$5
  shift 5
  arr=("$@")
  for i in ${arr[@]}
  do 
    pathToConnectivityDump=$baseConn$i$endConn
    pathToOutputFile=$baseOutput$i$endOutput
    python3 greatestTwoMasses_RMW_numSecondaryCycles.py -d $data -c $pathToConnectivityDump \
    -o $pathToOutputFile
  done
}

# dump indices
large1=(2 4 6 8 10 12 14 16 18)
large2=(2 4 6 8 10 12 14 16 18 20)
large3=(2 4 6 8 10 12 14 16)
large4=(2 4 6 8 10 12 14 16 18 20)
large5=()

medium1=(2 4 6 8 10 12 14 16)

small1=(2 4 6 8 10 12)

pathToDataFile=("/home/mattia/Documents/scriptDeveloper/anlammps/lammps/lammpsData/data.Min_large3"
                "/home/mattia/Documents/scriptDeveloper/anlammps/lammps/lammpsData/data.Min_medium1"
                "/home/mattia/Documents/scriptDeveloper/anlammps/lammps/lammpsData/data.Min_small1"
                )

basePathToConnectivityDump=("/home/mattia/Documents/emmy_backup/large_503K_1/dump/dump."
                            "/home/mattia/Documents/emmy_backup/large_503K_2/dump/dump."
                            "/home/mattia/Documents/emmy_backup/large_503K_3/dump/dump."
                            "/home/mattia/Documents/emmy_backup/large_503K_4/dump/dump."
                            "/home/mattia/Documents/emmy_backup/large_503K_5/dump/dump."
                            "/home/mattia/Documents/emmy_backup/medium_503K_1/dump/dump."
                            "/home/mattia/Documents/emmy_backup/small_503K_1/dump/dump."
                            )
                            
endPathToConnecvitiyDump="_btype_batom1_batom2_nptCuring_gelation"
                        
basePathToOutputFile=("/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large1/gelPoint/gelPointMeasurements"
                      "/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large2/gelPoint/gelPointMeasurements"
                      "/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large3/gelPoint/gelPointMeasurements"
                      "/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large4/gelPoint/gelPointMeasurements"
                      "/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/large5/gelPoint/gelPointMeasurements"
                      "/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/medium1/gelPoint/gelPointMeasurements"
                      "/home/mattia/Documents/epoxy/crossLinkingAndGelationPoint/systems/small1/gelPoint/gelPointMeasurements"
                      )
                      
endPathToOutputFile="_2nsNptCuring"

# small1
processSystem ${pathToDataFile[2]} ${basePathToConnectivityDump[6]} $endPathToConnecvitiyDump \
${basePathToOutputFile[6]} $endPathToOutputFile ${small1[@]} 

# medium1
processSystem ${pathToDataFile[1]} ${basePathToConnectivityDump[5]} $endPathToConnecvitiyDump \
${basePathToOutputFile[5]} $endPathToOutputFile ${medium1[@]}

# large1-5
processSystem ${pathToDataFile[0]} ${basePathToConnectivityDump[0]} $endPathToConnecvitiyDump \
${basePathToOutputFile[0]} $endPathToOutputFile ${large1[@]}

processSystem ${pathToDataFile[0]} ${basePathToConnectivityDump[1]} $endPathToConnecvitiyDump \
${basePathToOutputFile[1]} $endPathToOutputFile ${large2[@]}

processSystem ${pathToDataFile[0]} ${basePathToConnectivityDump[2]} $endPathToConnecvitiyDump \
${basePathToOutputFile[2]} $endPathToOutputFile ${large3[@]}

processSystem ${pathToDataFile[0]} ${basePathToConnectivityDump[3]} $endPathToConnecvitiyDump \
${basePathToOutputFile[3]} $endPathToOutputFile ${large4[@]}

processSystem ${pathToDataFile[0]} ${basePathToConnectivityDump[4]} $endPathToConnecvitiyDump \
${basePathToOutputFile[4]} $endPathToOutputFile ${large5[@]}



