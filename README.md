# Parameter Partition Algorithm
This provides a program to compute a, c slices of the RNA polytope using viennaRNA. 

## Install
This requires the ViennaRNA python package v2.7.0. It comes either from a full install of ViennaRNA or can be
installed on its own with pip [here](https://pypi.org/project/ViennaRNA/).

It also requires [numpy](https://pypi.org/project/numpy/) and [scipy](https://pypi.org/project/scipy/).

## Running
The main algorithm is param_partition.py. These are the options:

```
usage: TL_HPI [-h] [-b B] [--geometry] [--timing] [--unskewed] [--bounds n n n n] [--LP] [-d n]
              FASTA_PATH SAVE_PATH

Computes a,c slice for an RNA polytope

positional arguments:
  FASTA_PATH        Path to fasta file.
  SAVE_PATH         Directory to save output files.

options:
  -h, --help        show this help message and exit
  -b B              b value for computation. Default 0.
  --geometry        Save geometry of regions as SAVE_PATH/seqname_geometry.txt.
  --timing          Save timing data as SAVE_PATH/seqname_timing.txt
  --unskewed        Compute in unskewed space (i.e. excess branching)
  --bounds n n n n  Integer bounds in dckals/mol for slice in format a_min a_max c_min c_max. Default a>=0 plane (i.e. 0, 10000, -10000, 10000).
  --LP              Compute without --noLP option in viennaRNA.
  -d n              Set dandle mode (1 or 2).
```

By default param_partition will save a file seqname_sig_structs.txt consiting of all the signatures and corresponding structures found within the designated frame. 

Adding --geometry and --timing saves additional files. For instance,
```
python param_partition.py path/to/RNASEQ.fasta save/files/here --geometry --timing
```
will compute the a,c slice of the a >= plane  and save files `save/files/here/RNASEQ_geometry.txt`, `save/files/here/RNASEQ_timing.txt`, and the default `save/files/here/RNASEQ_sig_structs.txt`

# Additional Files
The directory partition_algo_testing contains the runTL2ggb.sage file which uses [sagemath](https://www.sagemath.org/) to create geogebra files for visualizing the computed partiton.