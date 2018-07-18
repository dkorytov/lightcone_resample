#!/bin/bash 

z_ranges=("0_1"  "1_2"  "2_3")
#steps_list=("499 487 475 464 453 442 432 421 411 401 392 382 373 365 355 338 331 323 315 307 300 293 286 279 272 266 253 247"  "247 241 235 230 224 219 213 208 203 198 194 189 184 180 176 171 167" "167 163 159 155 151 148 144 141 137 134 131 127 121")
steps_list=("499 487 475 464 453 442 432 421 411 401 392 382 373 365 355 338 331 323 315 307 300 293 286 279 272 266 253 247"  "247 241 235 230 224 219 213 208 203 198 194 189 184 180 176 171 167" "167 163 159 155 151 148 144 141 137 134 131 127")
gltcs_files=("/cosmo/homes/dkorytov/proj/protoDC2/output/ANL_box_v2.1.3_nocut_steps/\${step}.hdf5" "/cosmo/homes/dkorytov/data/Galacticus/galaxy_data/high_z/\${step}.hdf5" "/cosmo/homes/dkorytov/data/Galacticus/galaxy_data/high_z/\${step}.hdf5")
z_types=("low_z" "high_z" "high_z")
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
    #sed "s/#z_range#/${z_range}/g; s/#step_list#/${steps}/g; s/#gltcs_file#/${gltcs_file}/g" <cosmoDC2_v0.0.0_z_template.param >cosmoDC2_v0.0.0_z_${z_range}.param
    sed "s/#z_range#/${z_range}/g; s/#step_list#/${steps}/g; s@#gltcs_file#@${gltcs_file}@g; s/#z_type#/${z_type}/g" <protoDC2_v4.14.1_template.param >protoDC2_v4.14.1_z_${z_range}.param
    cat protoDC2_v4.14.1_z_${z_range}.param
done
