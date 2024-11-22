#ifndef EDGE_HPP
#define EDGE_HPP

#include <iostream>

#include <connections.hpp>

// representing single synaptic connection when generating network
struct Edge {
    int from_t, from_x, from_y, to_t, to_x, to_y;  // out==from, in==to
    const Connections::conn_t* syntype;
    double delay, weightFactor;

    // edge(int ft, int fx, int fy, int tt, int tx, int ty, Connections::conn*
    // s,
    //      double d)
    //     : from_t(ft), from_x(fx), from_y(fy), to_t(tt), to_x(tx), to_y(ty),
    //       syntype(s), delay(d), weightFactor(1) {}
    Edge(int ft, int fx, int fy, int tt, int tx, int ty,
         const Connections::conn_t* s, double d, double weightF = 1)
        : from_t(ft), from_x(fx), from_y(fy), to_t(tt), to_x(tx), to_y(ty),
          syntype(s), delay(d), weightFactor(weightF) {}
    // static int ct, cx, cy;  // cache for last to_ printed -> we want to print
    //                        // "In: " line sometimes

    Edge() : from_t(-1), from_x(-1), from_y(-1), to_t(-1), to_x(-1), to_y(-1),
          syntype(nullptr), delay(0.0), weightFactor(1) {}

    Edge& operator=(const Edge& other) {
        from_t = other.from_t;
        from_x = other.from_x;
        from_y = other.from_y;
        to_t = other.to_t;
        to_x = other.to_x;
        to_y = other.to_y;
        syntype = other.syntype;
        delay = other.delay;
        weightFactor = other.weightFactor;
        return *this;
    }

    inline bool same_to(const Edge& other) const {
        return (to_t == other.to_t && to_x == other.to_x &&
                to_y == other.to_y);
    }

    inline void header(std::ostream& out_stream) const {
        out_stream << "In: " << to_t << " " << to_x << " " << to_y << "\n";
    }

    inline void output(std::ostream& out_stream) const {
        if (weightFactor > 0) {
            out_stream << from_t << " " << from_x << " " << from_y << " "
                       << syntype->synapse << " "
                       << syntype->strength * weightFactor << " "
                       << syntype->mini_strength << " " << syntype->mini_freq
                       << " " << syntype->range << " " << delay << "\n";
        }
    }

    // inline int output0(std::ofstream& out_stream) const {
    //     int new_neuron = 0;  // did we enter new neurons while printing
    //                          // connections; for total neurons count.
    //     if (to_t != ct || to_x != cx || to_y != cy) {
    //         ct = to_t;
    //         cx = to_x;
    //         cy = to_y;
    //         new_neuron = 1;
    //     }
    //     if (weightFactor > 0) {
    //         out_stream << from_t << " " << from_x << " " << from_y << " "
    //                    << syntype->synapse << " "
    //                    << syntype->strength * weightFactor << " "
    //                    << syntype->mini_strength << " " << syntype->mini_freq
    //                    << " " << syntype->range << " " << delay << "\n";
    //     }
    //     return new_neuron;
    // }
    
    void dump() {
        std::cerr << from_t << " " << from_x << " " << from_y << " -> " << to_t
                  << " " << to_x << " " << to_y << " " << syntype->synapse
                  << "\n";
    }
   
    inline static bool cmp(const Edge& a, const Edge& b) {  // for sorting
        if (a.to_t != b.to_t)
            return (a.to_t < b.to_t);
        if (a.to_y != b.to_y)
            return (a.to_y < b.to_y);
        if (a.to_x != b.to_x)
            return (a.to_x < b.to_x);
        if (a.from_t != b.from_t)
            return (a.from_t < b.from_t);
        if (a.from_y != b.from_y)
            return (a.from_y < b.from_y);
        if (a.from_x != b.from_x)
            return (a.from_x < b.from_x);
        return false;
    }  // cmp
};  // edge

#endif  // !EDGE_HPP
