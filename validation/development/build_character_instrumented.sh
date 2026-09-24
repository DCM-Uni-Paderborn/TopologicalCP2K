#!/usr/bin/env bash
# Focused Apple ARM/GNU Fortran bounds, floating-point and undefined-behavior checks.
set -eu
root=$(pwd)
target=build-serial/character-debug
mkdir -p "$target"
flags=(-cpp -O1 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -fsanitize=undefined
       -fno-sanitize-recover=all -fopenmp -ffree-form -ffree-line-length-none)
includes=(-I"$target" -Ibuild-serial/src/mod_files -J"$target")
objects=()
for module in little_group corepresentations character_match; do
    object="$target/topology_${module}.o"
    /opt/homebrew/bin/gfortran "${flags[@]}" "${includes[@]}" \
        -c "src/topology_${module}.F" -o "$object"
    objects+=("$object")
done
/opt/homebrew/bin/gfortran "${flags[@]}" "${includes[@]}" \
    src/topology_character_unittest.F "${objects[@]}" \
    -Lbuild-serial/src -lcp2k.2026.2 "-Wl,-rpath,$root/build-serial/src" \
    -L/opt/homebrew/opt/openblas/lib -lopenblas -lstdc++ \
    -o "$target/character-reference"
