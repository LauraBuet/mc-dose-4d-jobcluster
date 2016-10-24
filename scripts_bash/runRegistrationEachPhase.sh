#!/bin/bash

# Script to perform registration-based motion field estimation in 4D CT data
# using the itk remote module VariationalRegistration

# -----------------------------------------
# Executable/binary files:

# BUILD PATH TS
BUILD_PATH="../../build/bin"
VARIATIONAL_REGISTRATION="${BUILD_PATH}/icnsVariationalRegistration"
VECTOR_FIELD_STATISTICS="${BUILD_PATH}/icnsVectorFieldStatistics"

echo '===================================================================='
echo '===================================================================='
echo ' Starting bash script'
echo '===================================================================='
echo '===================================================================='

# -----------------------------------------
# Test data:


TEST_DATA_INPUT_DIR="/data_l63/tsothman/NAS_RT/4D-CT/Klinische_Endpunkte"
TEST_DATA_INPUT_DIR_LOKAL="/data_l63/tsothman/ICNS_core/TestData"

SUBDIR_idx=0


for dir in "$TEST_DATA_INPUT_DIR"/*
	do
        #if [[ "$dir" =~ '+motionStat' ]]; then
    if [[ "$dir" =~ 'R2014010' ]]; then
		SUBDIR_DATA_PATH[SUBDIR_idx]="$dir"
		SUBDIR_idx=$((SUBDIR_idx+1))
	fi
done

NUMBER_OF_PATIENTS=${#SUBDIR_DATA_PATH[@]}
echo '----------------------------'
echo 'Number of folders:' ${NUMBER_OF_PATIENTS}
echo '----------------------------'

for ((iNumberOfPatients=0;iNumberOfPatients<${NUMBER_OF_PATIENTS};iNumberOfPatients++));
	do

	TEST_DATA_MASTER_DIR="${SUBDIR_DATA_PATH[${iNumberOfPatients}]}"
	TEST_4DCT="${TEST_DATA_MASTER_DIR}/4DCT"
	TEST_DATA_OUTPUT_DIR="${TEST_DATA_MASTER_DIR}/Motion_Statistics"

	echo '----------------------------'
	echo 'Patient folder:' ${TEST_DATA_MASTER_DIR}
	echo '----------------------------'
	echo '4DCT folder:' ${TEST_4DCT}
	echo '----------------------------'

		for (( i=0; i<=9; i++ ));
		do
			TEST_MOTION_FIELD="${TEST_DATA_OUTPUT_DIR}/Displ_2_to_${i}.mha"
			TEST_WARPED_IMAGE="${TEST_DATA_OUTPUT_DIR}/Warped_2_to_${i}.mha"

			REGISTRATION_CALL="${VARIATIONAL_REGISTRATION} -D ${TEST_4DCT} -y 2 -z ${i} -O ${TEST_MOTION_FIELD} -W ${TEST_WARPED_IMAGE}"

			echo $REGISTRATION_CALL
			$REGISTRATION_CALL
		done
done

exit

