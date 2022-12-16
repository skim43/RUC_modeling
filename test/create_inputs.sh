#!/bin/bash

#work_dir=`head working_dir.txt -n 1`

work_dir=${PWD}/input
t3d_dir=t3d
PGFem3D_input=PGFem3D_input

if [ ! -d ${work_dir}/${t3d_dir} ]; then
  echo "[${work_dir}/${t3d_dir}] doesn't exist."
  exit 1
fi

if [ ! -d ${work_dir}/${PGFem3D_input} ]; then
  echo "[${work_dir}/${PGFem3D_input}] doesn't exist."
  echo "Create [${work_dir}/${PGFem3D_input}]"
  mkdir ${work_dir}/${PGFem3D_input}
fi

dist=0.01
NP=16

cp json_files/* ${work_dir}/${PGFem3D_input}
cp model_params.in ${work_dir}/${PGFem3D_input}
cp ../local_makeset.pl ${work_dir}/${PGFem3D_input}

cd $work_dir/${PGFem3D_input}
cp ../$t3d_dir/*.t3d .
cp ../$t3d_dir/*.json .

fno=`ls *.t3d | wc -l`
for ((ia=1; ia<=$fno; ia++)); do
  in=`ls *.t3d | head -n $ia | tail -n 1`
  filename_extention="${in##*/}"
  filename="${filename_extention%.*}"
  cp ${filename}_bc.json Mechanical_bc.json
  cp ${filename}_mat.json material_mat.json
  rm -rf *.out partitions.*
  ./local_makeset.pl -np $NP -d $dist -f $filename
  rm Mechanical_bc.json material_mat.json
done

