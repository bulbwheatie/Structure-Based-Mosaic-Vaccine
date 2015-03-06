#! /bin/bash
# Master script to qsub iterations of rosaicNef.sh

script_dir=$1

qsub ${script_dir}/rosaicNef.sh nef_5_4.0 10
qsub ${script_dir}/rosaicNef.sh nef_5_4.1 10
qsub ${script_dir}/rosaicNef.sh nef_5_4.2 10

qsub ${script_dir}/rosaicNef.sh nef_5_5.0 30
qsub ${script_dir}/rosaicNef.sh nef_5_5.1 30
qsub ${script_dir}/rosaicNef.sh nef_5_5.2 30

qsub ${script_dir}/rosaicNef.sh nef_5_6.0 50
qsub ${script_dir}/rosaicNef.sh nef_5_6.1 50
qsub ${script_dir}/rosaicNef.sh nef_5_6.2 50