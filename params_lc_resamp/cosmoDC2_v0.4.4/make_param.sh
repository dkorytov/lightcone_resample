#!/bin/bash 

z_ranges=("0_1")
#steps_list=("499 487 475 464 453 442 432 421 411 401 392 382 373 365 355 338 331 323 315 307 300 293 286 279 272 266 253 247"  "247 241 235 230 224 219 213 208 203 198 194 189 184 180 176 171 167" "167 163 159 155 151 148 144 141 137 134 131 127 121")
steps_list=("499 487 475 464 453 442 432 421 411 401 392 382 373 365 355 347 338 331 323 315 307 300 293 286 279 272 266 259 253 247")
gltcs_files=("/gpfs/mira-fs0/projects/DarkUniverse_esp/dkorytov/data/Galacticus/low_z/galaxy_library/\${step}_mod.hdf5")

z_types=("low_z")
v1=0
v2=4
v3=4
mkdir logs
mkdir logs/old
rm *cosmo*
echo "mv logs/*.err logs/old/." > "run_all.sh"
echo "mv logs/*.out logs/old/." >> "run_all.sh"
for i in "${!z_ranges[@]}";do
    echo "${i}"
    z_range=${z_ranges[$i]}
    steps=${steps_list[$i]}
    gltcs_file=${gltcs_files[$i]}
    z_type=${z_types[$i]}
    echo ${z_range}
    echo ${steps}
    echo ${gltcs_file}
    echo ${z_type}
    sed "s/#z_range#/${z_range}/g; s/#step_list#/${steps}/g; s@#gltcs_file#@${gltcs_file}@g; s/#z_type#/${z_type}/g; s/#v1#/${v1}/g; s/#v2#/${v2}/g; s/#v3#/${v3}/g"<template.param >cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}.param
    sed "s/#z_range#/${z_range}/g; s/#step_list#/${steps}/g; s@#gltcs_file#@${gltcs_file}@g; s/#z_type#/${z_type}/g; s/#v1#/${v1}/g; s/#v2#/${v2}/g; s/#v3#/${v3}/g"<template.sh >run_cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}.sh
    chmod +x run_cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}.sh

    echo "qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/logs/${z_range}.out -e params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/logs/${z_range}.err  --debuglog=params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/logs/${z_range}.colbalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}.param" >> "run_all.sh"
done
chmod +x "run_all.sh"
