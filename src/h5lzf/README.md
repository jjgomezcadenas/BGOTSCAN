LZF filter of h5py (filter id 32000) as an HDF5 plugin, so that Julia's HDF5
library can read the LZF-compressed run files. Sources from
https://github.com/h5py/h5py/tree/master/lzf (BSD, see LICENSE.txt), unchanged.
`sh build.sh` produces `libh5lzf.dylib` (not in git); `JScan.register_lzf()`
points HDF5 to this directory. Rebuild after the Julia HDF5 artifact changes.
