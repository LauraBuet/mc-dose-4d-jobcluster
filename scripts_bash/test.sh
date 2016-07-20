#!/bin/bash

# Script to perform registration-based motion field estimation in 4D CT data
# using the itk remote module VariationalRegistration

# -----------------------------------------
# Executable/binary files:

# -----------------------------------------
# Test data:


TEST_DATA_INPUT_DIR="/data_l63/tsothman/NAS_RT/4D-CT/Klinische_Endpunkte"

SUBDIR_idx=0


for dir in "$TEST_DATA_INPUT_DIR"/*
	do
    if [[ "$dir" =~ '+motionStat' ]]; then
		SUBDIR_DATA_PATH[SUBDIR_idx]="$dir"
		SUBDIR_idx=$((SUBDIR_idx+1))
	fi
done

NUMBER_OF_PATIENTS=${#SUBDIR_DATA_PATH[@]}
echo '----------------------------'
echo 'Number of folders:' ${NUMBER_OF_PATIENTS}
echo '----------------------------'

for ((i=0;i<=$NUMBER_OF_PATIENTS;i++)); do

TEST_DATA_MASTER_DIR="${SUBDIR_DATA_PATH[i]}"
TEST_4DCT="${TEST_DATA_MASTER_DIR}/4DCT"

echo '----------------------------'
#echo 'Patient folder:' ${TEST_DATA_MASTER_DIR}
echo '----------------------------'
#echo '4DCT folder:' ${TEST_4DCT}

for dir in "$TEST_4DCT"/*
	do
	if [[ "$dir" =~ 'HERZPHASE' ]]; then
	echo 'move CT DATA'
	#find ${TEST_4DCT} -type f -exec mv --backup=numbered -t ${TEST_4DCT} {} +
	echo ${TEST_DATA_MASTER_DIR}
	fi
done
done

exit
