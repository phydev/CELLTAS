!#/bin/bash
SEED=-763928
nohup p085 085 $SEED 0.85 > log.085 &
nohup p090 090 $SEED 0.90 > log.090 &
nohup p095 095 $SEED 0.95 > log.095 &
echo 'running cells3d for the following porosities: 0.85, 0.90, 0.95'
echo $SEED
