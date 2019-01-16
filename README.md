:copyright: Tom Swinburne 2019
swinburne at cinam.univ-mrs.fr
Simple lattice KMC implementation for surface migration studies
- Simple python front end for easy prototyping
- C++ backend for reasonably fast execution

# Installation
- Compile C++ backend
```
cd lib/c++
make
```
# Execution
- Edit `kmc_run.py` for initial geometry, bonding energies, driving force and plotting options
```
python kmc_run.py
```

