#!/usr/bin/env bash
# Reproduce the focused Apple ARM/GNU Fortran instrumented build from a CP2K source root.
set -eu
root=$(pwd)
mkdir -p build-serial/affine-debug
flags=(-cpp -O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined
       -fno-sanitize-recover=all -fopenmp -ffree-form -ffree-line-length-none)
includes=(-Ibuild-serial/affine-debug -Ibuild-serial/src/mod_files -Jbuild-serial/affine-debug)
libraries=(-Lbuild-serial/src -lcp2k.2026.2 "-Wl,-rpath,$root/build-serial/src"
           -L/opt/homebrew/opt/openblas/lib -lopenblas -lstdc++)
/opt/homebrew/bin/gfortran "${flags[@]}" "${includes[@]}" -c src/topology_wyckoff.F \
    -o build-serial/affine-debug/topology_wyckoff.o
for name in site wyckoff; do
    /opt/homebrew/bin/gfortran "${flags[@]}" "${includes[@]}" \
        "src/topology_${name}_unittest.F" build-serial/affine-debug/topology_wyckoff.o \
        "${libraries[@]}" -o "build-serial/affine-debug/${name}-no-database"
done
