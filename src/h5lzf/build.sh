#!/bin/sh
# Build the h5py LZF filter as an HDF5 plugin (libh5lzf.dylib) against the
# HDF5 library used by the Julia project src/jscan. macOS only for now.
# Usage: sh src/h5lzf/build.sh        (from anywhere)
set -eu
here=$(cd "$(dirname "$0")" && pwd)
jscan="$here/../jscan"

lib=$(julia --project="$jscan" -e 'using HDF5; print(HDF5.API.libhdf5)')
hdf5=$(dirname "$(dirname "$lib")")
echo "HDF5 artifact: $hdf5"

# The artifact is an MPI-enabled build: its header includes mpi.h, which lives
# in the MPI artifact that libhdf5 links against (libmpi.12 = MPICH).
mpi_inc=""
for d in "$HOME"/.julia/artifacts/*/; do
    if [ -f "$d/include/mpi.h" ] && [ -e "$d/lib/libmpi.12.dylib" ]; then
        mpi_inc="$d/include"; break
    fi
done
[ -n "$mpi_inc" ] || { echo "mpi.h not found in the Julia depot"; exit 1; }
echo "MPI headers:   $mpi_inc"

cc -O2 -shared -fPIC -Wno-expansion-to-defined \
   -I"$hdf5/include" -I"$mpi_inc" -I"$here" -I"$here/lzf" \
   -L"$hdf5/lib" -lhdf5 -Wl,-rpath,"$hdf5/lib" \
   -o "$here/libh5lzf.dylib" \
   "$here/lzf_filter.c" "$here/lzf/lzf_c.c" "$here/lzf/lzf_d.c"
echo "built $here/libh5lzf.dylib"
