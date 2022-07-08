#!/bin/bash

source /sphenix/u/shuhang98/setup.sh~

echo "------------------setting up environment--------------------"
export Cur_dir=$(pwd)
echo "running area:" ${Cur_dir}
echo "-------------------------running----------------------------"
cd ${Cur_dir}
ls
root '/sphenix/user/shuhangli/hcalfullsim/macros/Fun4All_G4_sPHENIX.C(1e6)' > notes.log

rm G4sPHENIX_qa.root

echo "JOB COMPLETE!"
