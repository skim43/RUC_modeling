#!/bin/bash

if [ $# -lt 1 ]; then
  echo "Usage: ./make_matrix_property_1.sh [t3d_path]"
  echo "       [t3d_path] : directory containing t3d files from RoP2T3d"
  exit 1
fi

t3d_dir=$1

if [ ! -d $t3d_dir ]; then
  echo "[${t3d_dir}] doesn't exist."
  exit 1
fi

fno=`ls ${work_dir}/${t3d_dir}/*.t3d | wc -l`

for ((ia=1; ia<=fno; ia++)); do

  in=`ls ${t3d_dir}/*.t3d | head -n $ia | tail -n 1`
  filename_extention="${in##*/}"
  filename="${filename_extention%.*}"
#  extension="${filename_extention##*.}"

  region_id=`tail ${t3d_dir}/${filename}.out.regions -n 1 | grep -o '[0-9]*' | head -n 2 | tail -n 1`
  
  t3d_file=${t3d_dir}/${filename}.t3d
  echo ""
  echo "**************************************************"
  echo "change matrix property 0 -> 1: [${t3d_file}]"
  line=`grep "region $region_id" ${t3d_file}`
  m=0
  while true; do
    line=`echo $line | grep -o "property"`
    if [[ !  -z  $line ]]; then
      break
    else
      m=$[$m + 1]
      line=`grep -A ${m} "region $region_id" ${t3d_file} | tail -n 1`
    fi
  done
  line=`grep -A ${m} "region $region_id" ${t3d_file} | tail -n 1`
  old=${line%"property 0"}
  new="${old} property 1"
  sed -e "s|$line|${new}|" ${t3d_file} > ${t3d_file}_temp
  mv ${t3d_file}_temp ${t3d_file}
done
