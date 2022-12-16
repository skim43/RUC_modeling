#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage: ./do_rp2t3d.sh [Rocpack_path] [t3d_path]"
  echo "       [Rocpack_path] : directory containig Locpack files (e.g. Rocpack for directory Rocpack/*.rp2t3d)"
  echo "       [filebase]     : filebase name"
  echo "       [t3d_path]     : directory for saving t3d files"
  exit 1
fi

rp2t3d=/cswarm/tools/bin/RoP2T3d
rp2t3d_dir=$1
t3d_dir=$2

if [ ! -d $t3d_dir ]; then
  echo "[${t3d_dir}] doesn't exist."
  echo "Create [${t3d_dir}]"
  mkdir ${t3d_dir}
fi

fno=`ls $rp2t3d_dir/*.rp2t3d | wc -l`
for ((ia=1; ia<=$fno; ia++)); do
  in=`ls $rp2t3d_dir/*.rp2t3d | head -n $ia | tail -n 1`
  filename_extention="${in##*/}"
  filename="${filename_extention%.*}"
#  extension="${filename_extention##*.}"
  out=${t3d_dir}/${filename}
  $rp2t3d  -i $in -o $out -bc -npr
  echo $out
done
