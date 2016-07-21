## Wall drainage
Simulates wall drainage in PU foaming process

## Installation
### 1. Install fson library
The library is installed to `{HOME}/lib` and `{HOME}/include`
```
cd where-you-want-source-files
git clone https://github.com/josephalevin/fson.git
cd fson
make
make install
```
### 2. Install bspline library
The library is installed to `{HOME}/lib` and `{HOME}/include`
```
cd where-you-want-source-files
git clone https://github.com/japaf/bspline-fortran.git
cd bspline-fortran
cmake .
make
make install
```
### 3. Compile the model
```
cmake .
make
```

## Run
```
./walldrainage
```

## Postprocess
visualize results using plot.py
