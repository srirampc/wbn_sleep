#ifndef ARGS_HPP
#define ARGS_HPP

#include <algorithm>
#include <string>
#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <data_config.hpp>

struct AppArguments {
    CLI::App app;
    std::string config_file_path;
    std::string connect_summary_path;
    std::string data_cfg_path;
    std::string out_file_path;
    bool parallel;
    DataConfig data_cfg;

    AppArguments()
        : app(), config_file_path("mriNetworkNew.cfg"),
          connect_summary_path("connSummaryFileexample.txt"),
          out_file_path("network_output.txt"),
          parallel(false) {
        app.add_option("-c,--config", config_file_path,
                       "Path to the Input Config File")
            ->mandatory()
            ->check(CLI::ExistingFile);
        app.add_option("-s,--connect_summary", connect_summary_path,
                       "Path to Connection Summary File")
            ->capture_default_str();
        app.add_option("-d,--data_config", data_cfg_path,
                       "Path to Data Config File")
            ->check(CLI::ExistingFile);
        app.add_option("-o,--output_file", out_file_path, "Path to Output File")
            ->capture_default_str()
            ->mandatory();
        app.add_flag("-p,--parallel", parallel, "Flag if Parallel Run")
            ->capture_default_str();
    }

    void print() const {
        // Display arguments
        std::string arg1 =
            fmt::format(" Network Config. File       : {} ", config_file_path);
        std::string arg2 =
            fmt::format(" Output File                : {} ", out_file_path);
        std::string arg3 =
            fmt::format(" Connection Summary File    : {} ", connect_summary_path);
        std::string arg4 =
            fmt::format(" Data Config. File          : {} ", data_cfg_path);
        std::string arg5 =
            fmt::format(" Parallel Constr.           : {} ", parallel ? "true" : "false");
        std::size_t argsz =
            std::max({arg1.size(), arg2.size(), arg3.size(), arg4.size(), arg5.size()});
        fmt::print(" ┌{0:─^{1}}┐\n"
                   " │{2: <{1}}│\n"
                   " │{3: <{1}}│\n"
                   " │{4: <{1}}│\n"
                   " │{5: <{1}}│\n"
                   " │{6: <{1}}│\n"
                   " └{0:─^{1}}┘\n",
                   "", argsz, arg1, arg2, arg3, arg4, arg5);
    }
};

#endif // ARGS_HPP