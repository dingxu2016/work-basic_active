#!/bin/bash

dt=0.02
v0=0.1
iseed=66
natom=4096
phi=0.45
step=110000

run=/data1/dingxu/active_shear/normal/build/test
cd $run

/data1/dingxu/active_shear/normal/build/main<<EOF

$dt
$v0
$iseed
$natom
$phi
$step

EOF

exit 0;
