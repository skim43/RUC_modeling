#!/bin/bash

#work_dir=`head working_dir.txt -n 1`

work_dir=${PWD}
Rocpack_dir=Rocpack
t3d_dir=t3d

../do_rp2t3d.sh ${work_dir}/input/Rocpack ${work_dir}/input/t3d
../make_matrix_property_1.sh ${work_dir}/input/t3d
rm -f temp.out
