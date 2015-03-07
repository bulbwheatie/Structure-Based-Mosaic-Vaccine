#! /bin/bash
# Master script to qsub iterations of ${script_dir}/rosaicV1V2.sh

script_dir=$1

qsub ${script_dir}/rosaicV1V2.sh v1v2_6_3.0 10 "squared"
qsub ${script_dir}/rosaicV1V2.sh v1v2_6_3.1 10 "squared"
qsub ${script_dir}/rosaicV1V2.sh v1v2_6_3.2 10 "squared"

qsub ${script_dir}/rosaicV1V2.sh v1v2_6_4.0 30 "squared"
qsub ${script_dir}/rosaicV1V2.sh v1v2_6_4.1 30 "squared"
qsub ${script_dir}/rosaicV1V2.sh v1v2_6_4.2 30 "squared"

qsub ${script_dir}/rosaicV1V2.sh v1v2_6_5.0 70 "squared"
qsub ${script_dir}/rosaicV1V2.sh v1v2_6_5.1 70 "squared"
qsub ${script_dir}/rosaicV1V2.sh v1v2_6_5.2 70 "squared"

qsub ${script_dir}/rosaicV1V2.sh v1v2_5_10.0 10 "exponential"
qsub ${script_dir}/rosaicV1V2.sh v1v2_5_10.1 10 "exponential"
qsub ${script_dir}/rosaicV1V2.sh v1v2_5_10.2 10 "exponential"

qsub ${script_dir}/rosaicV1V2.sh v1v2_5_11.0 30 "exponential"
qsub ${script_dir}/rosaicV1V2.sh v1v2_5_11.1 30 "exponential"
qsub ${script_dir}/rosaicV1V2.sh v1v2_5_11.2 30 "exponential"

qsub ${script_dir}/rosaicV1V2.sh v1v2_5_12.0 70 "exponential"
qsub ${script_dir}/rosaicV1V2.sh v1v2_5_12.1 70 "exponential"
qsub ${script_dir}/rosaicV1V2.sh v1v2_5_12.2 70 "exponential"
