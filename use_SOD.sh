#!/bin/bash

# Autor Maciej Szyszko
# Skrypt do pre-processingu - uruchamiania SOD, przelatuje przez wszystkie foldery po prepare_for_SOD.py (można je będzie połączyć w jedno)

export PATH=$PATH:/home/maciej/Documents/ROOTSOD/sod-master/bin


#python3 /home/maciej/Desktop/SOD-PROGRAM/workflow/prepare_for_SOD.py 

cd /home/maciej/Documents/2D_Materials/4x4-ClusterExpansion/4x4x1_cell/komórki_3_150
# Outer loop: Find all subdirectories for example 2x2x1_cell
find . -maxdepth 1 -mindepth 1 -type d | while read -r subsubdirectory; do

    cd "$subsubdirectory" || exit

    # Execute your command here in the sub-subdirectory
    # For example, executing sod_comb.sh
    /home/maciej/Documents/ROOTSOD/sod-master/bin/sod_comb.sh
    
    echo "  Finished processing sub-subdirectory: $subsubdirectory"

    cd - || exit

done


