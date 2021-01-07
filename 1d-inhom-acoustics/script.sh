#!/bin/sh
sigma[0]=0.1
sigma[1]=0.2
sigma[2]=0.3
sigma[3]=0.4
sigma[4]=0.5
sigma[5]=0.6
sigma[6]=0.7
sigma[7]=0.8
sigma[8]=0.9
sigma[9]=1.0

epsilon[0]=0.01
epsilon[1]=0.02
epsilon[2]=0.03
epsilon[3]=0.04
epsilon[4]=0.05
k[0]=1
k[1]=2
k[2]=3
k[3]=4
k[4]=5
k[5]=6
k[6]=7
k[7]=8
k[8]=9
k[9]=10

sed -i -e "7s/.*/sigma0=${sigma[0]}/" avSlopePlot.py
sed -i -e "8s/.*/sigmaN=${sigma[9]}/" avSlopePlot.py
sed -i -e "9s/.*/epsilon0=${epsilon[0]}/" avSlopePlot.py
sed -i -e "10s/.*/epsilonN=${epsilon[4]}/" avSlopePlot.py
sed -i -e "11s/.*/k0=${k[0]}/" avSlopePlot.py
sed -i -e "12s/.*/kN=${k[9]}/" avSlopePlot.py
sed -i -e "13s/.*/Nx=10/" avSlopePlot.py
sed -i -e "14s/.*/Ny=5/" avSlopePlot.py
sed -i -e "15s/.*/Nz=10/" avSlopePlot.py

for i in {0..9}
do
  sed -i -e "34s/.*/    sigma=${sigma[i]}/" alfven.py
  for j in {0..4}
  do
    sed -i -e "36s/.*/    epsilon=${epsilon[j]}/" alfven.py
    for kappa in {0..9}
    do
      sed -i -e "37s/.*/    k=${k[kappa]}/" alfven.py
      for iter in {0..100}
      do
        python alfven.py && python energy.py
      done
      sed -i -e "4s/.*/sigma=${sigma[i]}/" avSlope.py
      sed -i -e "5s/.*/epsilon=${epsilon[j]}/" avSlope.py
      sed -i -e "6s/.*/k=${k[kappa]}/" avSlope.py
      python avSlope.py
      rm Slope.dat
    done
  done
done

python avSlopePlot.py
 
