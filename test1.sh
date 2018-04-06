#!/bin/bash

dt=0.1
v0=0.5
iseed=91
natom=2
box=1.2e1
step=2000

run=/data1/dingxu/active_shear/normal/build/test
cd $run

/data1/dingxu/active_shear/normal/build/main<<EOF

$dt
$v0
$iseed
$natom
$box
$step

EOF

exit 0;
