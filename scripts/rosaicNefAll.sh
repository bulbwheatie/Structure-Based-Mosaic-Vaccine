#! /bin/bash
# Master script to qsub iterations of rosaicNef.sh

script_dir=$1

qsub ${script_dir}/rosaicNef.sh nef_5_1.0 10
qsub ${script_dir}/rosaicNef.sh nef_5_1.1 10
qsub ${script_dir}/rosaicNef.sh nef_5_1.2 10

qsub ${script_dir}/rosaicNef.sh nef_5_2.0 50
qsub ${script_dir}/rosaicNef.sh nef_5_2.1 50
qsub ${script_dir}/rosaicNef.sh nef_5_2.2 50

qsub ${script_dir}/rosaicNef.sh nef_5_3.0 100
qsub ${script_dir}/rosaicNef.sh nef_5_3.1 100
qsub ${script_dir}/rosaicNef.sh nef_5_3.2 100