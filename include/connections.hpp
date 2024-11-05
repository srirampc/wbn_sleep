#ifndef CONNECTIONS_HPP
#define CONNECTIONS_HPP

#include <string>
#include <iostream>
#include <cassert>

#include <cell_types.hpp>

// connection types (not particular connections between neurons)
class Connections {
  public:
    // connection types between layers
    struct conn_t {
        int from, to;
        std::string synapse;
        double radius_min, radius_max, cs;  // radius; column size
        double probab, probab_oc;  // probability within  and outside of column
        std::string distribution;
        double strength,
            mini_strength,  // minimal strength of miniature synaptic potentials
            mini_freq;  // frequency of spontaneous minis (miniature synaptic
                        // potentials)
        int range;  // long/short range; this value has different meaning for
                    // generated (1==short range) and MRI network
    };  // conn
    int n;                        // number of connections
    conn_t C[CellTypes::MAXTYPES * CellTypes::MAXTYPES];  // connections
    int stats[CellTypes::MAXTYPES * CellTypes::MAXTYPES];  // connection statistics for debug purposes

    int find(int from, int to, std::string synapse, int range_type) {
        for (int i = 0; i < n; i++)
            if (C[i].from == from && C[i].to == to && C[i].synapse == synapse &&
                range_type == C[i].range)
                return i;
        std::cerr << "Missing connection " << from << "->" << to << " "
                  << synapse << ", range " << range_type << "\n";
        // exit(0); //connection type not entered in config file
        return -1;
    }  // find
    void dump(std::ostream& ConnSummaryFile) const {
        for (int i = 0; i < n; i++)
            dump(i, ConnSummaryFile);
    }  // dump
    void dump(int i, std::ostream& ConnSummaryFile) const {
        assert(i < n);
        ConnSummaryFile << C[i].from << " " << C[i].to << " " << C[i].synapse
                        << " " << C[i].radius_min << " " << C[i].radius_max
                        << " " << C[i].probab << "\n";
    }
    void dump_stats(const CellTypes& CT, std::ostream& ConnSummaryFile) const {
        for (int c = 0; c < n; c++) {
            ConnSummaryFile << c + 1 << " of " << n
                            << " , edges: " << stats[c] << " \ttype: (";
            CT.dump(C[c].from, ConnSummaryFile);
            ConnSummaryFile << ") -> (";
            CT.dump(C[c].to, ConnSummaryFile);
            ConnSummaryFile << ") ";
            dump(c, ConnSummaryFile);
        }
    }  // dump stats
    Connections() : n(0) {
        for (int i = 0; i < CellTypes::MAXTYPES * CellTypes::MAXTYPES; i++)
            stats[i] = 0;
    }
};  // Connections

#endif
