#!/bin/bash
echo "STARTING"
#User defined
#N1=(33 65 129 257 513)
N1=(33)
#N2=(33 65 129 257 513)
N2=(33)
#DT=(1e-4 5e-4 1e-5 5e-5 1e-6 5e-6 1e-7 5e-7 1e-8 5e-8 1e-9 5e-9 1e-10)
DT=(1e-8)
#Theta=(0 0.5 1)
Theta=(0.5)
#T=(5e-4 1e-3 1e-2)
T=(1e-3)
#muuu=(1e-3 1e-2 1e-1 1)
muuu=(1e-3)

#Loop over DT and EX
for per in ${!T[@]}
do
for m in ${!muuu[@]}
do
for th in ${!Theta[@]}
do
for t in ${!DT[@]}
do
for n1 in ${!N1[@]}
do
for n2 in ${!N2[@]}
do
    sed -i "/^dir_=/cdir_=cwd+'/Theta${Theta[$th]}-mu${muuu[$m]}-DT${DT[$t]}-${N1[$n1]}-${N2[$n2]}-T${T[$per]}'" escoamentoemduto_df.py
    sed -i "/^theta =/ctheta = ${Theta[$th]}" escoamentoemduto_df.py
    sed -i "/^N1    =/cN1    = ${N1[$n1]}" escoamentoemduto_df.py
    sed -i "/^N2    =/cN2    = ${N2[$n2]}" escoamentoemduto_df.py
    sed -i "/^dt    =/cdt    = ${DT[$t]}" escoamentoemduto_df.py
    sed -i "/^T     =/cT     = ${T[$per]}" escoamentoemduto_df.py
    sed -i "/^mu    =/cmu    = ${muuu[$m]}" escoamentoemduto_df.py
    
    python3 escoamentoemduto_df.py
done
done
done
done
done
done
echo "THE END"
