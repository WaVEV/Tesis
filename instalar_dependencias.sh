## System Libraries: GFORTRAN, BLAS, LAPACK, ARPACK
sudo apt-get update -qq
sudo apt-get install -y gfortran libopenblas-dev liblapack-dev libarpack2-dev
cd arpackpp
## BLAS:
./install-openblas.sh

## ARPACK
./install-arpack-ng.sh

## SUPERLU (version 5.0):
./install-superlu.sh

./install-suitesparse.sh