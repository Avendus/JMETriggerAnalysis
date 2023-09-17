#!/bin/bash

# Reco directories path which you want to investigate
DIR_PATHS=(
  "/eos/cms/store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/NoPUGTv4_GTv4_126X_mcRun3_2023_forPU65_v4-v2/2540000/"
  "/eos/cms/store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/NoPUGTv4_GTv4_126X_mcRun3_2023_forPU65_v4-v2/2810000/"
  "/eos/cms/store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/NoPUGTv4_GTv4_126X_mcRun3_2023_forPU65_v4-v2/30000/"
  "/eos/cms/store/mc/Run3Winter23Reco/SinglePionGun_E0p2to200/GEN-SIM-RECO/NoPUGTv4_GTv4_126X_mcRun3_2023_forPU65_v4-v2/30001/"
)

# Loop through each directory in the list.
for DIR_PATH in "${DIR_PATHS[@]}"; do
    echo "Hello world"
    # Find .root files within the current directory.
    ROOT_FILES=($(find $DIR_PATH -type f -name "*.root"))

    # Print the current directory path.
    echo "Set RECO sample directory path= $DIR_PATH"

    # Execute the dasgoclient command for each .root file.
    for FILE in "${ROOT_FILES[@]}"; do
        echo "Command: dasgoclient -json -query="parent file=$FILE""
        FILE_NAME=$(basename $FILE)
        RESULT=$(dasgoclient -json -query="parent file=$FILE")
        echo "Finding RECO's parents for [ $FILE_NAME ] : [ $RESULT ]"
    done
done
