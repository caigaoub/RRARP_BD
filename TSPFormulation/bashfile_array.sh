#！/bin/bash
# This is an example
DATPATH=/media/caigao/E:/Project1Linux/Git/RRARP_BD/InstanceGenerator/ret/RRARP_15/RRARP_instance_n_15_E_
DSTNZ=12;
for i in {1..5}
do
  ./bin/MainRun $DATPATH$i.dat $DSTNZ
done

DATPATH=/media/caigao/E:/Project1Linux/Git/RRARP_BD/InstanceGenerator/ret/RRARP_15/RRARP_instance_n_15_M_
for i in {1..5}
do
  ./bin/MainRun $DATPATH$i.dat $DSTNZ
done

DATPATH=/media/caigao/E:/Project1Linux/Git/RRARP_BD/InstanceGenerator/ret/RRARP_15/RRARP_instance_n_15_H_
for i in {1..5}
do
  ./bin/MainRun $DATPATH$i.dat $DSTNZ
done
