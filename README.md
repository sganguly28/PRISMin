# PRISMin

Example repo for performing minimal PRISM fits with DUNE flux

## Prerequisits

Requires ROOT and a C++17-capable compiler

## Build like:

```
git clone https://github.com/luketpickering/PRISMin.git
cd PRISMin; mkdir build
cd build
cmake ..
make
```

## Run like;

```
cd /path/to/PRISMin/build
./PRISMin LBNF_TDRFlux_Nov19.root
```

This will write out a `Coeffs.pdf` file with a few example PRISM fits.