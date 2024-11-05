#ifndef NET_CONFIG_HPP
#define NET_CONFIG_HPP
//
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
//
#include <nlohmann/json.hpp>
using njson = nlohmann::json;
//
#include <cell_types.hpp>
#include <connections.hpp>
#include <utils.hpp>
//
struct NetConfig {
    std::ofstream conn_summary_fstream;
    std::ofstream out_fstream;

    CellTypes CT;
    Connections CN;
    bool both_hemispheres_detected = false;
    bool is_valid = false;

private:
    std::string get_input_line(FILE* f, std::string& line) {  // NOLINT
        char buf[4096];
        std::string tok1;

        while (1) {
            if (!fgets(buf, sizeof(buf), f)) {
                return line = "";
            }  // EOF
            line = buf;
            if (line.empty())
                continue;
            std::stringstream parser(buf);
            if (!(parser >> tok1))
                continue;  // whitespace line
            if (tok1[0] == '#')
                continue;  // just comment
            break;
        }
        return tok1;
    }

    // in case we load ucsd data for both hemispheres. signalized in
    // config via existence of "r" prefix in cell type name
    int parse_config_cfg(char const* file) {
        std::string line, tok;
        double x, y, ratio_x, ratio_y;
        FILE* f = fopen(file, "r");
        assert(f);

        // Types
        tok = get_input_line(f, line);
        assert(tok == "[Types]");
        // ratios
        tok = get_input_line(f, line);
        std::stringstream parser(line);
        assert(parser >> ratio_x);
        assert(parser >> ratio_y);
        // types enumeration
        tok = get_input_line(f, line);
        while (tok != "[Connections]") {
            assert(CT.size() < CellTypes::MAXTYPES);
            //
            parser.str(line);
            CellTypes::cell_t pcx;
            assert(parser >> pcx.name >> pcx.x >> pcx.y);
            //
            if (!both_hemispheres_detected && pcx.name.find("r") == 0) {
                std::cerr
                    << "Configuration for both-hemisphere connections detected "
                       "(rXXX "
                       "cell type found).\n";
                both_hemispheres_detected = true;
                // only x will be used, because MRI data are 1 dimensionally
                // indexed, and it must be set _after_ loading MRI data
                if (ratio_x != 1) {
                    std::cerr
                        << "Ratio will be used in different way because both "
                           "hemispheres detected "
                        << std::endl;
                    CT.set_mri_ratio(ratio_x);
                    ratio_x = 1;
                }
            }
            CT.add_type(pcx);
            tok = get_input_line(f, line);
        }  // types

        // ratios after loop because MRI scenario
        CT.apply_reduction(ratio_x, ratio_y);

        // Short range connections
        std::string from, to, syn, dist;
        double rad, cs, p, poc, str, mf, ms;
        assert(tok == "[Connections]");
        tok = get_input_line(f, line);
        while (tok != "[Longrange]") {
            assert(CN.n < CellTypes::MAXTYPES * CellTypes::MAXTYPES);
            parser.str(line);
            assert(parser >> from >> to >> syn >> rad >> cs >> p >> poc >>
                   dist >> str >> ms >> mf);

            Connections::conn_t& c = CN.C[CN.n];  // current connection
            c.from = CT.get_type(from);
            c.to = CT.get_type(to);
            c.synapse = syn;
            c.radius_min = 0;
            c.radius_max = rad;
            c.cs = cs;
            c.probab = p;
            c.probab_oc = poc;
            c.distribution = dist;
            c.strength = str;
            c.mini_strength = ms;
            c.mini_freq = mf;
            c.range = 1;
            CN.n++;

            tok = get_input_line(f, line);
        }  // connections

        // Long range connections
        double rad_min, rad_max;
        assert(tok == "[Longrange]");
        tok = get_input_line(f, line);
        while (!feof(f)) {
            assert(CN.n < CellTypes::MAXTYPES * CellTypes::MAXTYPES);
            parser.str(line);
            assert(parser >> from >> to >> syn >> rad_min >> rad_max >> p >>
                   dist >> str >> ms >> mf);

            Connections::conn_t& c = CN.C[CN.n];  // current connection
            c.from = CT.get_type(from);
            c.to = CT.get_type(to);
            c.synapse = syn;
            c.radius_min = rad_min;
            c.radius_max = rad_max;
            c.probab = p;
            c.range = 0;
            CN.n++;

            // TODO(x) we should change rad_min for small (1D) networks to be
            // bigger or use specific function

            tok = get_input_line(f, line);
        }  // connections

        fclose(f);
        return 0;
    }

    int parse_config_json(char const* config_file) {
        std::ifstream cfg_ptr(config_file);
        njson cfg_data = njson::parse(cfg_ptr, nullptr, false, true);
        //
        // Types
        double ratio_x, ratio_y;
        ratio_x = cfg_data["ratio_x"];
        ratio_y = cfg_data["ratio_y"];
        // std::cout << "CTN : " << cfg_data["types"].size() << std::endl;
        for (auto& ct_entry : cfg_data["types"]) {
            assert(CT.size() < CellTypes::MAXTYPES);
            std::string type_name = ct_entry["name"];
            if (type_name[0] == '#')
                continue;
            CT.add_type(CellTypes::cell_t(ct_entry["dims"]["x"],
                                          ct_entry["dims"]["y"], type_name));
        }
        // assert(CT.n <= cfg_data["types"].size());
        // ratios after loop because MRI scenario
        CT.apply_reduction(ratio_x, ratio_y);

        // CT.dump(&std::cout);
        //
        // Connections
        for (auto& cn_entry : cfg_data["connections"]) {
            assert(CN.n < CellTypes::MAXTYPES * CellTypes::MAXTYPES);
            std::string from_t = cn_entry["from"], to_t = cn_entry["to"];
            if (from_t[0] == '#' || to_t[0] == '#')
                continue;
            Connections::conn_t& c = CN.C[CN.n];  // current connection
            c.from = CT.get_type(from_t);
            c.to = CT.get_type(to_t);
            c.synapse = cn_entry["type"];
            auto& prop_map = cn_entry["property_map"];
            c.radius_min = 0;
            c.radius_max = prop_map["radius"];
            c.cs = prop_map["col_sci"];
            c.probab = prop_map["prob"];
            c.probab_oc = prop_map["prob_oc"];
            c.distribution = prop_map["distribution"];
            c.strength = prop_map["strength"];
            c.mini_strength = prop_map["mini_s"];
            c.mini_freq = prop_map["mini_f"];
            c.range = 1;
            CN.n++;
        }
        // CN.dump(&std::cout);
        //
        // Long Range Connections
        for (auto& cn_entry : cfg_data["long_range"]) {
            assert(CN.n < CellTypes::MAXTYPES * CellTypes::MAXTYPES);
            std::string from_t = cn_entry["from"], to_t = cn_entry["to"];
            if (from_t[0] == '#' || to_t[0] == '#')
                continue;
            Connections::conn_t& c = CN.C[CN.n];  // current connection
            c.from = CT.get_type(from_t);
            c.to = CT.get_type(to_t);
            c.synapse = cn_entry["type"];
            auto& prop_map = cn_entry["property_map"];
            c.radius_min = prop_map["radius_min"];
            c.radius_max = prop_map["radius_max"];
            c.probab = prop_map["prob"];
            c.range = 0;
            CN.n++;
        }
        //
        if (CT.has_both_hemispheres()) {
            std::cerr
                << "Configuration for both-hemisphere connections detected "
                   "(rXXX "
                   "cell type found).\n";
            both_hemispheres_detected = true;
            // only x will be used, because MRI data are 1 dimensionally
            // indexed, and it must be set _after_ loading MRI data
            if (ratio_x != 1) {
                std::cerr << "Ratio will be used in different way because both "
                             "hemispheres "
                             "detected\n";
                CT.set_mri_ratio(ratio_x);
                ratio_x = 1;
            }
        }
        return 0;
    }

    int parse_config(std::string cfg_file) {
        int rcx = 0;
        if (StringUtil::ends_with(cfg_file, ".cfg")) {
            rcx = parse_config_cfg(cfg_file.c_str());
            is_valid = true;
        } else if (StringUtil::ends_with(cfg_file, ".json")) {
            try {
                rcx = parse_config_json(cfg_file.c_str());
                is_valid = true;
            } catch (njson::exception& e) {
                std::cout << "Failed to parse Config File " << cfg_file << ":"
                          << e.what() << std::endl;
                rcx = -1;
            }
        }
        return rcx;
    }
public:
    explicit NetConfig(std::string cfg_file, std::string connect_summary_file,
                       std::string output_file)
        : conn_summary_fstream(connect_summary_file), out_fstream(output_file) {
        // Config Parse
        parse_config(cfg_file);
    }

  private:
    NetConfig() {}
};

#endif  // NET_CONFIG_HPP
