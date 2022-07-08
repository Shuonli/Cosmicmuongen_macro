#!/bin/bash

export TargetDir="$PWD"/condorRun


if [ -d ${TargetDir} ]; then
  rm -rf ${TargetDir}/OutDir*
else
  mkdir ${TargetDir}
fi


njob=500

for((i=0;i<$njob;i++));
do

  mkdir ${TargetDir}/OutDir$i
  export WorkDir="${TargetDir}/OutDir$i"
  echo "WorkDir:" ${WorkDir}

  pushd ${WorkDir}

  cp -v "$PWD"/../../CondorRun.sh CondorRunTC$i.sh
  cp "$PWD"/../../G4Setup_sPHENIX.C .
  cp "$PWD"/../../DisplayOn.C .
  cp "$PWD"/../../init_gui_vis.mac .
  cp "$PWD"/../../vis.mac .
  chmod +x CondorRunTC$i.sh

  cat >>ff.sub<< EOF
+JobFlavour                   = "workday"
transfer_input_files          = ${WorkDir}/CondorRunTC$i.sh,${WorkDir}/G4Setup_sPHENIX.C,${WorkDir}/DisplayOn.C,,${WorkDir}/init_gui_vis.mac,,${WorkDir}/vis.mac
Executable                    = CondorRunTC$i.sh
Universe                      = vanilla
Notification                  = Never
GetEnv                        = True
Priority                      = +20
Output                        = test.out
Error                         = test.err
Log                           = test.log
Notify_user                   = blair.daniel.seidlitz@cern.ch

Queue
EOF

  condor_submit ff.sub
  popd
done
