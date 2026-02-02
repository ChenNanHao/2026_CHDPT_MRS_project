#!/bin/bash

input_folder="./T1_nifti/achd"
output_folder="./Sienax_output/achd"

mkdir -p "${output_folder}"

for nii_file in "${input_folder}"/*.nii; do
    filename=$(basename "${nii_file}" .nii)
    sienax "${nii_file}" -o "${output_folder}/${filename}"
done

