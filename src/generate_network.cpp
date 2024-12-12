#include <cassert>
#include <fmt/format.h>
#include <iostream>
#include <string>
//
#include <nlohmann/json.hpp>
using njson = nlohmann::json;
//
#include <cell_types.hpp>
#include <connections.hpp>
#include <data_config.hpp>
#include <edge.hpp>
#include <edge_list.hpp>
#include <net_config.hpp>
#include <source_data.hpp>
#include <utils.hpp>
#include <args.hpp>

int generate_network(const AppArguments& in_args) {
    timer run_timer;
    NetConfig net_cfg(in_args.config_file_path, in_args.connect_summary_path,
                      in_args.out_file_path);
    if (!net_cfg.is_valid) {
        std::cout << "Failed to parse Network Config File: "
                  << in_args.config_file_path << std::endl;
        return -1;
    }
    PRINT_RUNTIME_MEMUSED(run_timer, "Init. Config : ", std::cout);

    run_timer.reset();
    if (in_args.data_cfg_path.length() > 0) {
        assert(net_cfg.both_hemispheres_detected);

        SourceData net_sd(net_cfg);  // Init. Source Data
        net_sd.load_source_data(in_args.data_cfg, in_args.parallel_read);

        // Reduction of number neurons, if ratio in config
        net_cfg.CT.MRI_reduce();

        EdgeList net_ce(net_cfg, net_sd, in_args.parallel_build);      // Init. Edge List
        net_ce.generate_multilayer_edges();

        timer w_timer;
        net_ce.write_connections(net_cfg.out_fstream, 0);
        PRINT_RUNTIME_MEMUSED(w_timer, "Write Connections : ", std::cout);
        net_ce.clear();
    } else {
        // Default Run
        SourceData net_sd(net_cfg);  // Init. Source Data
        EdgeList net_ce(net_cfg, net_sd, in_args.parallel);  // Init. Edge List
        net_ce.generate_radom_edges(false);  //  random connectivity
        // CE.dump_from_neuron_edges(4,20,20);
        // CE.dump_from_neuron_edges(0,0,0,2);
        net_ce.write_connections(net_cfg.out_fstream);
        // CT.dump(); CN.dump();
    }
    PRINT_RUNTIME_MEMUSED(run_timer, "Gen. Network (complete) : ", std::cout);
    return 0;
}

int main(int argc, char* argv[]) {
    AppArguments in_args;
    CLI11_PARSE(in_args.app, argc, argv);  // Parse command line arguments
    if(in_args.parallel) {
         in_args.parallel_read = in_args.parallel_build = in_args.parallel;
    }
    in_args.print();
    if (in_args.data_cfg_path.length() > 0) {
        try {
            in_args.data_cfg = DataConfig::parse(in_args.data_cfg_path);
        } catch (njson::exception& e) {
            std::cerr << "Failed to parse Data Config File "
                      << in_args.data_cfg_path << ":" << e.what() << std::endl;
            return -1;
        }
        in_args.data_cfg.print();
    }
    return generate_network(in_args);
}
