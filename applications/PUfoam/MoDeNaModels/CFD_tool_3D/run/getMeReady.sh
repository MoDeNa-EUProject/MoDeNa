#!/bin/bash

# delete everything but the good stuff!
ls > list
egrep -v 'constant|system|results|0$|getMeReady.sh|log*|README|run.sh|prep.sh|case.foam' list > list2
mv list2 list
rm -rf $(<list)

rm -rf 1??
rm -rf 2??
rm -rf 3??
rm -rf ?0
rm -rf 4??
rm -rf 5??
rm -rf 6??

# find everything but blockMeshDict and remove them

find ./constant/polyMesh/ ! -name "blockMeshDict*" -type f -exec rm -r '{}' \;

# prepare the 0 directory for run
cd ./0
rm -v alpha.gas 
cp -v alphagas.org alpha.gas

rm -v mZero 
cp -v mZero.org mZero
rm -v M0 
cp -v M0.org M0

rm -v mOne 
cp -v mOne.org mOne
rm -v M1 
cp -v M1.org M1

rm -v mTwo 
cp -v mTwo.org mTwo
rm -v M2 
cp -v M2.org M2

rm -v mThree 
cp -v mThree.org mThree
rm -v M3 
cp -v M3.org M3

rm -v mFour 
cp -v mFour.org mFour
rm -v M4 
cp -v M4.org M4

rm -v mFive 
cp -v mFive.org mFive
rm -v M5 
cp -v M5.org M5

rm -v wBA_l
cp -v wBA_l.org wBA_l

rm -v rho_gas 
cp -v rho_gas.org rho_gas
rm -v rho_foam
cp -v rho_foam.org rho_foam

rm -v muFoam
cp -v muFoam.org muFoam
rm -v muMixture
cp -v muMixture.org muMixture

rm -v muFoamCorr
cp -v muFoamCorr.org muFoamCorr
rm -v XNCO
cp -v XNCO.org XNCO

rm -v rhoPU
cp -v rhoPU.org rhoPU

rm -v weight0
cp -v weight0.org weight0
rm -v weight1
cp -v weight1.org weight1
rm -v weight2
cp -v weight2.org weight2

rm -v node0
cp -v node0.org node0
rm -v node1
cp -v node1.org node1
rm -v node2
cp -v node2.org node2

rm -v TS
cp -v TS.org TS

# cp -v mZero.org mZero
cd ./..
