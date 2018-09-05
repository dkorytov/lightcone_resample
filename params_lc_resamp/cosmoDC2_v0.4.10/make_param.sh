#!/bin/bash 

z_ranges=("0_1" "1_2" "2_3")
steps_list=("499 487 475 464 453 442 432 421 411 401 392 382 373 365 355 347 338 331 323 315 307 300 293 286 279 272 266 259 253 247"  "247 241 235 230 224 219 213 208 203 198 194 189 184 180 176 171 167" "167 163 159 155 151 148 144 141 137 134 131 127 124 121")
gltcs_files=("/gpfs/mira-fs0/projects/DarkUniverse_esp/dkorytov/data/Galacticus/low_z/galaxy_library/\${step}_mod.hdf5" "/gpfs/mira-fs0/projects/DarkUniverse_esp/dkorytov/data/Galacticus/high_z/galaxy_library/\${step}_mod.hdf5" "/gpfs/mira-fs0/projects/DarkUniverse_esp/dkorytov/data/Galacticus/high_z/galaxy_library/\${step}_mod.hdf5")
z_types=("low_z" "high_z" "high_z")
#healpix_groups=("533 534 564 565 566 596 597" "598 599 628 629 630 660 661")
healpix_groups=("533")
healpix_group_names=("a" "b" "c" "d" "e" "f")
v1=0
v2=4
v3=10
mkdir logs
mkdir logs/old
rm *cosmo*
echo "#!/bin/bash" >"run_all.sh"
echo "#!/bin/bash" >"run_all_login.sh"
echo "mv logs/*.err logs/old/." >> "run_all.sh"
echo "mv logs/*.out logs/old/." >> "run_all.sh"
for i in "${!z_ranges[@]}";do
    echo "z ${i}"
    z_range=${z_ranges[$i]}
    steps=${steps_list[$i]}
    gltcs_file=${gltcs_files[$i]}
    z_type=${z_types[$i]}
    echo ${z_range}
    echo ${steps}
    echo ${gltcs_file}
    echo ${z_type}
    for j in "${!healpix_groups[@]}";do
	healpix_name=${healpix_group_names[j]}
	healpix_group=${healpix_groups[j]}
	sed "s/#z_range#/${z_range}/g; s/#step_list#/${steps}/g; s@#gltcs_file#@${gltcs_file}@g; s/#z_type#/${z_type}/g; s/#healpix_group#/${healpix_group}/g; s/#v1#/${v1}/g; s/#v2#/${v2}/g; s/#v3#/${v3}/g"<template.param >cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}_hp:${healpix_name}.param
	sed "s/#z_range#/${z_range}/g; s/#step_list#/${steps}/g; s@#gltcs_file#@${gltcs_file}@g; s/#z_type#/${z_type}/g; s/#healpix_name#/${healpix_name}/g; s/#v1#/${v1}/g; s/#v2#/${v2}/g; s/#v3#/${v3}/g"<template.sh >run_cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}_hp:${healpix_name}.sh
	chmod +x run_cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}_hp:${healpix_name}.sh
    echo "qsub -n 1 -t 720 -o params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/logs/${z_range}_${healpix_name}.out -e params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/logs/${z_range}_${healpix_name}.err  --debuglog=params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/logs/${z_range}_${healpix_name}.cobalt  ./lc_resample.py params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}_hp:${healpix_name}.param" >> "run_all.sh"
    echo "./lc_resample.py params_lc_resamp/cosmoDC2_v${v1}.${v2}.${v3}/cosmoDC2_v${v1}.${v2}.${v3}_z_${z_range}_hp:${healpix_name}.param" >> "run_all_login.sh"
    done
done
chmod +x "run_all.sh"
chmod +x "run_all_login.sh"
mv logs/*.err logs/old/.
mv logs/*.out logs/old/.
mv logs/*.cobalt logs/old/.
mv logs/*.colbalt logs/old/.
