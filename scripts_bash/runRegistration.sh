#!/bin/bash

# Script to perform registration-based motion field estimation in 4D CT data
# using the itk remote module VariationalRegistration

# -----------------------------------------
# Executable/binary files:

# BUILD PATH TS
BUILD_PATH="../../build/bin"
VARIATIONAL_REGISTRATION="${BUILD_PATH}/icnsVariationalRegistration"
VECTOR_FIELD_STATISTICS="${BUILD_PATH}/icnsVectorFieldStatistics"

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
echo 'Patient folder:' ${TEST_DATA_MASTER_DIR}
echo '----------------------------'
echo '4DCT folder:' ${TEST_4DCT}

for dir in "$TEST_4DCT"/*
	do
	if [[ "$dir" =~ 'HERZPHASE' ]]; then
	echo 'move CT DATA'
	find ${TEST_4DCT} -type f -exec mv --backup=numbered -t ${TEST_4DCT} {} +
	echo 'done moving CT DATA'
	fi
done


TEST_ITV="${TEST_DATA_MASTER_DIR}/Motion_Statistics"

TEST_DATA_OUTPUT_DIR="${TEST_DATA_MASTER_DIR}/Motion_Statistics"

for file in "$TEST_DATA_OUTPUT_DIR"/*
	do
	if [[ "$file" =~ 'MotionStatistics.txt' ]]; then
	echo 'deleting old motion statistics'
	rm "${TEST_DATA_OUTPUT_DIR}/MotionStatistics.txt"
	echo 'done deleting'
	fi
done

TEST_MOTION_FIELD="${TEST_DATA_OUTPUT_DIR}/Displ.mha"
TEST_WARPED_IMAGE="${TEST_DATA_OUTPUT_DIR}/Warped.mha"

TEST_LOGFILE="${TEST_DATA_OUTPUT_DIR}/MotionStatistics.txt"

	for file in "$TEST_ITV"/*
	do
		if [[ "$file" =~ 'mhd' ]]; then
			ITV_PATH="${file}"
			echo '----------------------------'
			echo 'ITV data path:' ${ITV_PATH}
			echo '----------------------------'

# -----------------------------------------
# Program calls:

# First: run registration:
if [ ! -f $TEST_MOTION_FIELD ]; then

REGISTRATION_CALL="${VARIATIONAL_REGISTRATION} -D ${TEST_4DCT} -y 0 -z 5 -O ${TEST_MOTION_FIELD} -W ${TEST_WARPED_IMAGE}"
echo $REGISTRATION_CALL
$REGISTRATION_CALL

fi

# Then analyze displacement field:
MOTION_ANALYSIS_CALL="${VECTOR_FIELD_STATISTICS} -I ${TEST_MOTION_FIELD} -S ${ITV_PATH} -L ${TEST_LOGFILE}"
echo $MOTION_ANALYSIS_CALL
$MOTION_ANALYSIS_CALL
		fi
	done

done

exit

