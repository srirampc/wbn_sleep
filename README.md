## Generate Network for Whole brain Modelling

### Build Requisites

1. Linux Distribution (Tested in Ubuntu 24.04 and Redhat Linux v9)
2. make  (Tested with v4.3)
3. C++ 17 compiler (Tested with gcc versions: 11.4.1, 12.3 and 13.2)
4. CMake (version >3.2; Optional, if required to build tests with cmake)


### Building the Source
Clone the repository with the following command 

    git clone https://github.com/srirampc/wbn_sleep.git


To build the executable using make, run `make` as below

    make

### Usage 

The executable is generated in the build/ directory. Run as below

    ./build/generate_network -c ./config/mriNetworkNew.json -d ./config/data_config.json -s output/connect_summary.txt -o output/random_connection_info.txt

Usage of generate_network is as shown below:

    Usage: ./build/generate_network [OPTIONS]
    
    Options:
      -h,--help                   Print this help message and exit
      -c,--config TEXT:FILE REQUIRED
                                  Path to the Input Config File
      -s,--connect_summary TEXT [connSummaryFileexample.txt]
                                  Path to Connection Summary File
      -d,--data_config TEXT:FILE  Path to Data Config File
      -o,--output_file TEXT [network_output.txt]  REQUIRED
                                  Path to Output File
      -p,--parallel               Flag if Parallel File Read and Build Network
      --parallel-read             Flag if Parallel File Read only
      --parallel-build             Flag if Parallel Build Only

Parallel implementation consists of two parts:

1. Reading the input data file (specifically the weights files).
2. Building the network (i.e., construction of the edge list).

The weights input file is read in parallel if `--parallel-read` flag is set.
The edge list construction is done in parallel if `--parallel-build` flag is set.
If the `-p` or `--parallel` flag is set, both the parallellizaton is run.
As of this now, parallel network building always produces a random network. 
If the network construction needed to be sequential, please use `--parallel-read`
flag.

#### Input/Config Files

The following two input config files are required. Examples are provided in the config/ directory.
1. mriNetworkNew.json : Defines the network neuron types and connections.
2. data_config.json : Defines the paths to data files, whose format is as below; Make sure the data files pointed in the files are available 
.

```json
     {
         "left_subnet" : "./data/Map_642_To_10242_30-Aug-2017_LH.txt",
         "right_subnet": "./data/Map_642_To_10242_30-Aug-2017_RH.txt",
         "left_right_map": "./data/Map_LH_RH_10242_30-Aug-2017.txt",
         "left_right_map_small": "./data/Map_LH_RH_642_30-Aug-2017.txt",
         "dist3d_probability": "./data/dist_inMeters_prob_LH_onesOnDiag.txt",
         "weight_factor": "./data/weights_LH.txt",
         "weight_factor_inpy": "./data/weights_INPY_both_hemi.txt",
         "thalamus_cortex_distance": "./data/thalCort_dist_inMeters_LH.txt",
         "intra_thalamus_distance":  "./data/dist_TC_LH_zeroDiag.txt"
     }
```

The program generates two ouput files

1. Summary file: Location provided with the -s argument, contains the summary.
2. Output network files: Location provided with the -o argument.
