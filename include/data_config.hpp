#ifndef DATA_CONFIG_HPP
#define DATA_CONFIG_HPP

#include <string>
#include <algorithm>
#include <fstream>
//
#include <fmt/format.h>
#include <nlohmann/json.hpp>
using njson = nlohmann::json;

struct DataConfig {
    std::string left_subnet;
    std::string right_subnet;
    std::string left_right_map;
    std::string left_right_map_small;
    std::string dist3d_probability;
    std::string weight_factor;
    std::string weight_factor_inpy;
    std::string thalamus_cortex_distance;
    std::string intra_thalamus_distance;

    void print() {
        std::string arg1 =
            fmt::format(" left_subnet               : {} ", left_subnet);
        std::string arg2 =
            fmt::format(" right_subnet              : {} ", right_subnet);
        std::string arg3 =
            fmt::format(" left_right_map            : {} ", left_right_map);
        std::string arg4 = fmt::format(" left_right_map_small      : {} ",
                                       left_right_map_small);
        std::string arg5 =
            fmt::format(" dist3d_probability        : {} ", dist3d_probability);
        std::string arg6 =
            fmt::format(" weight_factor             : {} ", weight_factor);
        std::string arg7 =
            fmt::format(" weight_factor_inpy        : {} ", weight_factor_inpy);
        std::string arg8 = fmt::format(" thalamus_cortex_distance  : {} ",
                                       thalamus_cortex_distance);
        std::string arg9 = fmt::format(" intra_thalamus_distance   : {} ",
                                       intra_thalamus_distance);
        std::size_t argsz = std::max({arg1.size(), arg2.size(), arg3.size(),
                                      arg4.size(), arg5.size(), arg6.size(),
                                      arg7.size(), arg8.size(), arg9.size()});
        fmt::print(" ┌{0:─^{1}}┐\n"
                   " │{2: <{1}}│\n"
                   " │{3: <{1}}│\n"
                   " │{4: <{1}}│\n"
                   " │{5: <{1}}│\n"
                   " │{6: <{1}}│\n"
                   " │{7: <{1}}│\n"
                   " │{8: <{1}}│\n"
                   " │{9: <{1}}│\n"
                   " │{10: <{1}}│\n"
                   " └{0:─^{1}}┘\n",
                   "", argsz, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                   arg9);
    }

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(DataConfig, left_subnet, right_subnet,
                                   left_right_map, left_right_map_small,
                                   dist3d_probability, weight_factor, 
                                   weight_factor_inpy, thalamus_cortex_distance,
                                   intra_thalamus_distance);

    // initialize
    static DataConfig parse(const std::string& pathToDataConfig){
        std::ifstream d_stream(pathToDataConfig);
        njson c_jdata = njson::parse(d_stream, nullptr, false, true);
        return c_jdata.template get<DataConfig>();
    }
};

#endif // DATA_CONFIG_HPP 
