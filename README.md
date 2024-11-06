## Generate Network for Whole brain Modelling

Clone the repository with the following command 

    git clone https://github.com/srirampc/wbn_sleep.git


To build the executable using make, run `make` as below

    make release

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

 The program generate two ouput files

1. Summary file: Provied with the -s argument, Outputs the summary.
2. Output netowrk files: Provied with the -o argument.
