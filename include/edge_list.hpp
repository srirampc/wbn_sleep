#ifndef EDGE_LIST_HPP
#define EDGE_LIST_HPP

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <numeric>
#include <omp.h>
#include <random>
#include <string>
#include <vector>

#include <cell_types.hpp>
#include <connections.hpp>
#include <edge.hpp>
#include <net_config.hpp>
#include <source_data.hpp>

struct ConnectPair {
    int FX;
    int FY;
    int TX;
    int TY;
    std::string fromCellType;  // fr type
    std::string toCellType;    // to type
    int fromLayer;
    int toLayer;  // Layer
    int from_IN;
    int to_IN;  // IN PY

    ConnectPair(const NetConfig& ncfg, const SourceData& sdt, int conn_idx) {
        const CellTypes& CT = ncfg.CT;
        const Connections& CN = ncfg.CN;
        const auto& mapLayerNameId_PY = sdt.mapLayerNameId_PY_Ref();
        const auto& mapLayerNameId_IN = sdt.mapLayerNameId_IN_Ref();
        //
        FX = CT[CN.C[conn_idx].from].x;
        FY = CT[CN.C[conn_idx].from].y;  // input layer dimension
        TX = CT[CN.C[conn_idx].to].x;
        TY = CT[CN.C[conn_idx].to].y;  // output layer dimension
        assert(FY == 1);
        assert(TY == 1);  // MRI data are stored as 1D indexes
        //
        fromCellType = CT[CN.C[conn_idx].from].name;
        toCellType = CT[CN.C[conn_idx].to].name;
        fromLayer = -1;
        toLayer = -1;
        from_IN = 0;
        to_IN = 0;
        auto mplFromIter = mapLayerNameId_PY.find(CT[CN.C[conn_idx].from].name);
        if (mplFromIter != mapLayerNameId_PY.end()) {
            fromLayer = mplFromIter->second;
        } else if (mapLayerNameId_IN.find(CT[CN.C[conn_idx].from].name) !=
                   mapLayerNameId_IN.end()) {
            from_IN = 1;
        }
        auto mplToIter = mapLayerNameId_PY.find(CT[CN.C[conn_idx].to].name);
        if (mplToIter != mapLayerNameId_PY.end()) {
            toLayer = mplToIter->second;
        } else if (mapLayerNameId_IN.find(CT[CN.C[conn_idx].to].name) !=
                   mapLayerNameId_IN.end()) {
            to_IN = 1;
        }
    }
};

struct MultiLayerConfig {
    bool right_hemisphere;
    int indOffset;
    int large;
    int small;

    MultiLayerConfig(const NetConfig& net_cfg, bool right) {
        const CellTypes& CT = net_cfg.CT;
        right_hemisphere = right;
        small = CT[CT.get_type("TC")].x;
        large = CT[CT.get_type("CX")].x;
        indOffset = 0;
        if (net_cfg.both_hemispheres_detected && right) {
            small = CT[CT.get_type("rTC")].x;
            large = CT[CT.get_type("rCX")].x;
            indOffset = CT[CT.get_type("CX")].x;
        }
    }
};

class EdgeList {
  public:
    const Edge NullEdge = Edge(-1, -1, -1, -1, -1, -1, nullptr, 0.0);

  private:
    NetConfig& net_cfg;
    const SourceData& src_data;
    std::vector<Edge> edges;
    bool parallel_construction;

    inline int add_edge(Edge an_edge) {
        if (an_edge.syntype == nullptr)
            return -1;
        edges.emplace_back(an_edge);
        return edges.size() - 1;
    }

    inline int add_edge(int ft, int fx, int fy, int tt, int tx, int ty,
                        const Connections::conn_t* s, double d,
                        double weightF = 1) {
        if (s == nullptr)
            return -1;
        edges.emplace_back(ft, fx, fy, tt, tx, ty, s, d, weightF);
        return edges.size() - 1;
    }

    inline int add_edge_to(std::vector<Edge>* vEdgePtr, int ft, int fx, int fy,
                           int tt, int tx, int ty, const Connections::conn_t* s,
                           double d, double weightF = 1) {
        if (s == nullptr)
            return -1;
        vEdgePtr->emplace_back(ft, fx, fy, tt, tx, ty, s, d, weightF);
        return vEdgePtr->size() - 1;
    }

    void generate_3D_connections(bool right_hemisphere) {
        //
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        auto& Subnet = src_data.SubnetRef();
        bool both_hemispheres_detected = net_cfg.both_hemispheres_detected;

        // generates connection using external distance matrix
        std::cerr << "Generating connection file from pre-set 3D MRI data, "
                     "hemisphere: "
                  << right_hemisphere << "\n";

        int small = CT[CT.get_type("TC")].x;
        int large = CT[CT.get_type("CX")].x;

        if (both_hemispheres_detected && right_hemisphere)
            small = CT[CT.get_type("rTC")].x, large = CT[CT.get_type("rCX")].x;

        std::cerr << "Full set(CX): " << large
                  << ", Subnet(TC,RE,IN...): " << small << "\n";

        std::vector<int> nedges(CN.n, 0);  // edges counter for statistics
        for (int c = 0; c < CN.n; c++) {   // layer connections
            if (right_hemisphere == !CT.right_hemisphere(CN.C[c].from))
                continue;  // must be intrahemispheric
            if (right_hemisphere == !CT.right_hemisphere(CN.C[c].to))
                continue;
            int FX = CT[CN.C[c].from].x;
            int FY = CT[CN.C[c].from].y;  // input layer dimension
            int TX = CT[CN.C[c].to].x;
            int TY = CT[CN.C[c].to].y;  // output layer dimension
            assert(FY == 1);
            assert(TY == 1);  // MRI data are stored as 1D indexes

            for (int fxi = 0,
                     fx = (FX == large) ? 0 : Subnet[right_hemisphere][0];
                 (FX == large) ? fx < FX : Subnet[right_hemisphere][fxi] != -1;
                 (FX == large) ? fx++ : fx = Subnet[right_hemisphere][++fxi]) {
                // scaling is actually hidden in Subnet array, identity in large
                // case, map in small case
                int fy = 0;
                int sc_x = fx;
                int sc_y = fy;

                for (int txi = 0,
                         tx = (TX == large) ? 0 : Subnet[right_hemisphere][0];
                     (TX == large) ? tx < TX
                                   : Subnet[right_hemisphere][txi] != -1;
                     (TX == large)
                         ? tx++
                         : tx = Subnet[right_hemisphere]
                                      [++txi]) {  // particular connections
                                                  // between neurons
                    int ty = 0;
                    if (sc_x == tx && sc_y == ty &&
                        CT[CN.C[c].from].name == CT[CN.C[c].to].name)
                        continue;  // never connect neuron to itself
                    // distance is measured between 'From' neuron projected into
                    // 'TO' layer
                    //  printf("sc_x: %d, tx: %d, right_hem:%d ",
                    //  sc_x,tx,right_hemisphere);
                    double d = src_data.dist3D(sc_x, tx, right_hemisphere);
                    double delayTime = src_data.getScaledDelay(d);
                    if (d >= CN.C[c].radius_min && d <= CN.C[c].radius_max) {
                        // column is no more checked
                        if (CN.C[c].range == 1) {
                            if (MathUtil::rand01() <= CN.C[c].probab_oc)
                                check(CT, CN.C[c].from,
                                      (FX == large) ? fx : fxi) &&
                                    check(CT, CN.C[c].to,
                                          (TX == large) ? tx : txi) &&
                                    (add_edge(CN.C[c].from,
                                              (FX == large) ? fx : fxi, fy,
                                              CN.C[c].to,
                                              (TX == large) ? tx : txi, ty,
                                              &CN.C[c], delayTime),
                                     nedges[c]++);
                        } else {  // longrange
                            if (MathUtil::rand01() <= CN.C[c].probab)
                                check(CT, CN.C[c].from,
                                      (FX == large) ? fx : fxi) &&
                                    check(CT, CN.C[c].to,
                                          (TX == large) ? tx : txi) &&
                                    (add_edge(CN.C[c].from,
                                              (FX == large) ? fx : fxi, fy,
                                              CN.C[c].to,
                                              (TX == large) ? tx : txi, ty,
                                              &CN.C[c], delayTime),
                                     nedges[c]++);
                        }  // longrange
                    }  // radius
                }  // tx
            }  // fx
        }  // c

        dump_summary(net_cfg.conn_summary_fstream, nedges);
    }

    Edge multilayer_edge(int conn_idx,
                         const ConnectPair& cpair,            // Connection Pair
                         const MultiLayerConfig& mlcx,        // ML Config
                         const int& from_x, const int& to_x,  // From & to
                         std::mt19937* pgen) const {
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        const Vector4d<float>& weight_factor = src_data.weight_factor_Ref();
        const Vector2d<float>& weight_factor_INPY =
            src_data.weight_factor_INPY_Ref();
        const Vector1d<float>& TC_Dist = src_data.TC_Dist_Ref();
        const float maxDelay_TCm = src_data.maxDelay_TCm;
        const float maxDelay_TCc = src_data.maxDelay_TCc;

        std::mt19937& gen = *pgen;

        const int from_y = 0;  // all y are zero
        const int to_y = 0;

        if (from_x == to_x && from_y == to_y &&
            CT[CN.C[conn_idx].from].name == CT[CN.C[conn_idx].to].name)
            return NullEdge;  // never connect neuron to itself
        // distance is measured between 'From' neuron projected into
        // 'TO' layer
        //  printf("sc_x: %d, tx: %d, right_hem:%d ",
        //  sc_x,tx,right_hemisphere);
        double d =
            src_data.dist_only(from_x + mlcx.indOffset, to_x + mlcx.indOffset);
        // if(d<0) continue;
        double delayTime = 0.0;  // = getScaledDelay(d);

        // use spherical Thalamus
        size_t foundFromTC = cpair.fromCellType.find("TC");
        size_t foundFromRE = cpair.fromCellType.find("RE");
        size_t foundToTC = cpair.toCellType.find("TC");
        size_t foundToRE = cpair.toCellType.find("RE");

        if (((foundFromTC != std::string::npos) ||
             (foundFromRE != std::string::npos)) &&
            ((foundToTC != std::string::npos) ||
             (foundToRE != std::string::npos))) {
            d = src_data.intraThalamicDistOnSphere(from_x, to_x);
            // no indOffsetfor the right hemi within thalamus
        }

        if (d < 0)
            return NullEdge;

        // use dealyes only on CX->CX and CX->IN connections
        size_t foundFromCX = cpair.fromCellType.find("CX");
        size_t foundToCX = cpair.toCellType.find("CX");
        size_t foundToIN = cpair.toCellType.find("IN");
        if (!((foundFromCX != std::string::npos) &&
              ((foundToCX != std::string::npos) ||
               (foundToIN != std::string::npos)))) {
            delayTime = 0.0;
        } else {
            delayTime = src_data.getScaledDelay(d);
        }

        double currEdgeProb =
            src_data.prob_only(from_x + mlcx.indOffset, to_x + mlcx.indOffset);
        double clampProb = (CN.C[conn_idx].probab_oc * currEdgeProb < 1)
                               ? CN.C[conn_idx].probab_oc * currEdgeProb
                               : 1;

        if (cpair.fromLayer > -1 && cpair.toLayer > -1) {  // PY->PY connections
            // Random :: Draw from a bernoulli_distribution
            std::bernoulli_distribution distribution(clampProb);
            if (distribution(gen)) {
                double currEdgeWeightFactor = weight_factor(
                    from_x + mlcx.indOffset, to_x + mlcx.indOffset,
                    cpair.fromLayer, cpair.toLayer);
                if (d <= 0.0001) {
                    currEdgeWeightFactor *= 5;  // local weights scaled up to
                                                // test their influence
                } else if (d > CN.C[conn_idx].radius_max) {
                    currEdgeWeightFactor = 0;
                }
                if (currEdgeWeightFactor > 0) {
                    if (check(CT, CN.C[conn_idx].from, from_x) &&
                        check(CT, CN.C[conn_idx].to, to_x))
                        return Edge(CN.C[conn_idx].from, from_x, from_y,
                                    CN.C[conn_idx].to, to_x, to_y,
                                    &CN.C[conn_idx], delayTime,
                                    currEdgeWeightFactor);
                }
            }
        } else if (cpair.fromLayer > -1 &&
                   cpair.to_IN > 0) {  // PY->IN connections
            // Random :: Draw from a bernoulli_distribution
            std::bernoulli_distribution distribution(clampProb);
            if (distribution(gen)) {
                double currEdgeWeightFactor =
                    weight_factor(from_x + mlcx.indOffset,
                                  to_x + mlcx.indOffset, cpair.fromLayer,
                                  cpair.fromLayer);  // Since PY->IN are only
                                                     // within the same layer,
                                                     // use 'fromLayer' twice
                if (currEdgeWeightFactor > 0) {
                    if (check(CT, CN.C[conn_idx].from, from_x) &&
                        check(CT, CN.C[conn_idx].to, to_x))
                        return Edge(CN.C[conn_idx].from, from_x, from_y,
                                    CN.C[conn_idx].to, to_x, to_y,
                                    &CN.C[conn_idx], delayTime,
                                    currEdgeWeightFactor);
                }
            }
        } else if (cpair.from_IN > 0 &&
                   cpair.toLayer > -1) {  // IN->PY connections
            if (d >= CN.C[conn_idx].radius_min &&
                d <= CN.C[conn_idx].radius_max) {
                if (CN.C[conn_idx].range == 1) {
                    if (MathUtil::rand01() <= CN.C[conn_idx].probab_oc) {
                        double currEdgeWeightFactor = weight_factor_INPY(
                            from_x + mlcx.indOffset, to_x + mlcx.indOffset);
                        if (currEdgeWeightFactor > 0) {
                            if (check(CT, CN.C[conn_idx].from, from_x) &&
                                check(CT, CN.C[conn_idx].to, to_x))
                                return Edge(CN.C[conn_idx].from, from_x, from_y,
                                            CN.C[conn_idx].to, to_x, to_y,
                                            &CN.C[conn_idx], delayTime,
                                            currEdgeWeightFactor);
                        }
                    }
                }
            }  // radius
        } else {  // all other types of connections, all are based
                  // on the disc model
            if (d >= CN.C[conn_idx].radius_min &&
                d <= CN.C[conn_idx].radius_max) {
                // use TC_Dist only on thalamocortical connections
                // (CX,TC,IN,RE)
                size_t foundFromCX = cpair.fromCellType.find("CX");
                size_t foundFromIN = cpair.fromCellType.find("IN");
                size_t foundToCX = cpair.toCellType.find("CX");
                size_t foundToIN = cpair.toCellType.find("IN");

                size_t foundFromTC = cpair.fromCellType.find("TC");
                size_t foundFromRE = cpair.fromCellType.find("RE");
                size_t foundToTC = cpair.toCellType.find("TC");
                size_t foundToRE = cpair.toCellType.find("RE");

                // different maxDealy on Core and Matrix connections
                // be careful with condition since there is CX5a
                // which is not matrix (TCa)
                size_t foundFromMatrix = cpair.fromCellType.find("a");
                size_t foundToMatrix = cpair.toCellType.find("a");

                if (((foundFromCX != std::string::npos) ||
                     (foundFromIN != std::string::npos)) &&
                    ((foundToTC != std::string::npos) ||
                     (foundToRE != std::string::npos))) {
                    if ((foundToMatrix != std::string::npos)) {
                        // TODO(x) for both hemis use: sc_x +
                        // indOffset, tx
                        // + indOffset
                        delayTime = src_data.getScaled_TC_Delay(TC_Dist[from_x],
                                                                maxDelay_TCm);
                    } else {
                        delayTime = src_data.getScaled_TC_Delay(TC_Dist[from_x],
                                                                maxDelay_TCc);
                    }
                }

                if (((foundToCX != std::string::npos) ||
                     (foundToIN != std::string::npos)) &&
                    ((foundFromTC != std::string::npos) ||
                     (foundFromRE != std::string::npos))) {
                    if ((foundFromMatrix != std::string::npos)) {
                        // TODO(x) for both hemis use: sc_x +
                        // indOffset, tx
                        // + indOffset
                        delayTime = src_data.getScaled_TC_Delay(TC_Dist[to_x],
                                                                maxDelay_TCm);
                    } else {
                        delayTime = src_data.getScaled_TC_Delay(TC_Dist[to_x],
                                                                maxDelay_TCc);
                    }
                }

                if (delayTime < 0) {
                    delayTime = 0;
                }

                double currEdgeWeightFactor = 1;
                if ((CT[CN.C[conn_idx].from].name == "TC") &&
                    (CT[CN.C[conn_idx].to].name == "CX4")) {
                    currEdgeWeightFactor = 3;
                }
                // column is no more checked
                if (CN.C[conn_idx].range == 1) {
                    if (MathUtil::rand01() <= CN.C[conn_idx].probab_oc) {
                        if (check(CT, CN.C[conn_idx].from, from_x) &&
                            check(CT, CN.C[conn_idx].to, to_x))
                            return Edge(CN.C[conn_idx].from, from_x, from_y,
                                        CN.C[conn_idx].to, to_x, to_y,
                                        &CN.C[conn_idx], delayTime,
                                        currEdgeWeightFactor);
                    }
                } else {  // longrange
                    if (MathUtil::rand01() <= CN.C[conn_idx].probab) {
                        if (check(CT, CN.C[conn_idx].from, from_x) &&
                            check(CT, CN.C[conn_idx].to, to_x))
                            return Edge(CN.C[conn_idx].from, from_x, from_y,
                                        CN.C[conn_idx].to, to_x, to_y,
                                        &CN.C[conn_idx], delayTime,
                                        currEdgeWeightFactor);
                    }
                }  // longrange
            }  // radius
        }
        return NullEdge;
    }

    void
    par_generate_multilayer_connections(const ConnectPair& cpx,
                                        const MultiLayerConfig& mlcx, int c,
                                        std::vector<std::mt19937*> vrGenPtr,
                                        std::vector<int>::iterator nEdgesItr) {
        const Vector1d<float>& Subnet =
            src_data.Subnet_Ref(mlcx.right_hemisphere);
        std::vector<int> vt_nedges, vt_starts;
        int nedges = 0;
        int edge_start = edges.size();
        std::size_t ncandidates = std::size_t(cpx.FX) * std::size_t(cpx.TX);
#pragma omp parallel default(none)                                             \
    shared(cpx, mlcx, c, vrGenPtr, Subnet, vt_nedges, vt_starts, nedges,       \
               edge_start, ncandidates)
        {
            int n_threads = omp_get_num_threads(),
                thread_id = omp_get_thread_num(),
                t_start(block_low(thread_id, n_threads, ncandidates)),
                t_end(block_high(thread_id, n_threads, ncandidates));
            std::vector<Edge> vlt_edges;  // edges local to this thread
#pragma omp single
            {
                vt_nedges.resize(n_threads, 0);
                vt_starts.resize(n_threads, 0);
            }
            // thread local random generator
            for (int tidx = t_start; tidx <= t_end; tidx++) {
                int from_x = tidx / cpx.TX;
                int to_x = tidx % cpx.TX;
                if ((cpx.FX == mlcx.large) ? from_x >= cpx.FX
                                           : Subnet[from_x] == -1)
                    continue;
                //
                if ((cpx.TX == mlcx.large) ? to_x >= cpx.TX
                                           : Subnet[to_x] == -1)
                    continue;
                Edge edg = multilayer_edge(c, cpx, mlcx, from_x, to_x,
                                           vrGenPtr[thread_id]);
                if (edg.syntype != nullptr) {
                    vlt_edges.emplace_back(edg);
                }
            }
            vt_nedges[thread_id] = vlt_edges.size();
            //
#pragma omp barrier
#pragma omp single
            {
                nedges = std::accumulate(vt_nedges.begin(), vt_nedges.end(), 0);
                edges.resize(edge_start + nedges);
                std::exclusive_scan(vt_nedges.begin(), vt_nedges.end(),
                                    vt_starts.begin(), 0);
            }
            std::copy(vlt_edges.begin(), vlt_edges.end(),
                      edges.begin() + edge_start + vt_starts[thread_id]);
        }
        *nEdgesItr = nedges;
    }

    void
    par_generate_multilayer_connections(const MultiLayerConfig& mlcx,
                                        std::vector<int>::iterator nEdgeItr) {
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        std::vector<std::mt19937*> vrGenPtr;
#pragma omp parallel default(none) shared(vrGenPtr)
        {
            int n_threads = omp_get_num_threads();
            int thread_id = omp_get_thread_num();
#pragma omp single
            {
                vrGenPtr.resize(n_threads);
            }
            std::random_device rd;
            vrGenPtr[thread_id] = new std::mt19937(rd());
        }

        for (int c = 0; c < CN.n; c++) {  // layer connections
            if (mlcx.right_hemisphere == !CT.right_hemisphere(CN.C[c].from))
                continue;  // must be intrahemispheric
            if (mlcx.right_hemisphere == !CT.right_hemisphere(CN.C[c].to))
                continue;
            ConnectPair cpx(net_cfg, src_data, c);
            // std::cout << "R:" << c << " :: (" << cpx.FX << ", " << cpx.FY
            //           << ") ; (" << cpx.TX << ", " << cpx.TY << ")"
            //           << std::endl;
            assert((cpx.fromLayer < 6) && (cpx.toLayer < 6));
            par_generate_multilayer_connections(cpx, mlcx, c, vrGenPtr,
                                                nEdgeItr + c);
        }  // c
        //
        for (auto ix = 0u; ix < vrGenPtr.size(); ix++)
            delete vrGenPtr[ix];
    }

    void generate_multilayer_connections(const ConnectPair& cpx,
                                         const MultiLayerConfig& mlcx, int c,
                                         std::mt19937* pgen,
                                         std::vector<int>::iterator nedgeItr) {
        const Vector1d<float>& Subnet =
            src_data.Subnet_Ref(mlcx.right_hemisphere);
        for (int fx = 0;
             (cpx.FX == mlcx.large) ? fx < cpx.FX : Subnet[fx] != -1; fx++) {
            // scaling is actually hidden in Subnet array, identity in
            // large case, map in small case
            // int fy = 0;
            // int sc_x = fx;
            // int sc_y = fy;
            // if ((cpx.FX == mlcx.large) ? fx < cpx.FX : Subnet[fx] == -1)
            //    break;
            for (int tx = 0;
                 (cpx.TX == mlcx.large) ? tx < cpx.TX : Subnet[tx] != -1;
                 tx++) {
                // particular connections between neurons
                if (add_edge(multilayer_edge(c, cpx, mlcx, fx, tx, pgen)) != -1)
                    (*nedgeItr)++;
            }  // tx
        }  // fx
    }

    void generate_multilayer_connections(const MultiLayerConfig& mlcx,
                                         std::vector<int>::iterator nEdgeItr) {
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        std::vector<int> nedges(CN.n, 0);  // edges counter for statistics
        // // uncomment to use a non-deterministic seed:
        //    std::random_device rd;
        //    std::mt19937 gen(rd());
        std::mt19937 gen(1937);
        for (int c = 0; c < CN.n; c++) {  // layer connections
            if (mlcx.right_hemisphere == !CT.right_hemisphere(CN.C[c].from))
                continue;  // must be intrahemispheric
            if (mlcx.right_hemisphere == !CT.right_hemisphere(CN.C[c].to))
                continue;
            ConnectPair cpx(net_cfg, src_data, c);
            assert((cpx.fromLayer < 6) && (cpx.toLayer < 6));
            // std::cout << "R:" << c << "/" << CN.n << " :: (" << cpx.FX << ",
            // "
            //           << cpx.FY << ") ; (" << cpx.TX << ", " << cpx.TY << ")"
            //           << std::endl;
            generate_multilayer_connections(cpx, mlcx, c, &gen, nEdgeItr + c);
        }  // c
    }  // generate_MultiLayer_connections

    void generate_multilayer_connections(bool right_hemisphere) {
        std::cerr << "Generating connection file from pre-set 3D MRI data, new "
                     "MultiLayer connections, hemisphere: "
                  << right_hemisphere << std::endl;  //
        MultiLayerConfig mlcfg(net_cfg, right_hemisphere);
        std::cerr << "Full set(CX): " << mlcfg.large
                  << "; Subnet(TC,RE,IN...): " << mlcfg.small << std::endl;
        std::vector<int> nedges(net_cfg.CN.n, 0);  // edges counter
        if (parallel_construction) {
            par_generate_multilayer_connections(mlcfg, nedges.begin());
        } else {
            generate_multilayer_connections(mlcfg, nedges.begin());
        }
        dump_summary(net_cfg.conn_summary_fstream, nedges);
    }  // generate_multilayer_connections

    // load network from the file (fixed structure) we
    // get from ucsd folks approximating MRI topology
    void load_MRI_network(const char* file) {
        //
        // Tables for mapping specific synapses to MAP types for the hybrid
        // model Tables are derived from the loading code, but the whole thing
        // can be explained much simpler at the end: CX & IN are Map based
        // neurons and _any_ synapse on them must be converted into MAP type;
        // lazy to simplify this now. The mechanism needs to be checked in case
        // short X long range connection differ in synapses used. All to all
        // connections, e.g. CX* -> CX*
        const std::string rule_all[][4] = {
            {"CX", "CX", "AMPA_D2", "AMPAMap_D1"},
            {"CX", "CX", "NMDA_D1", "NMDAMap_D1"},
            {"TC", "CX", "AMPA_D1", "AMPAMap_D1"},
            {"TC", "IN", "AMPA_D1", "AMPAMap_D1"}};

        // Paired connections, eg. CX/a/6 -> IN/a/6 (not CX6->INa)
        const std::string rule_paired[][4] = {
            {"CX", "IN", "AMPA_D2", "AMPAMap_D1"},
            {"CX", "IN", "NMDA_D1", "NMDAMap_D1"},
            {"IN", "CX", "GABA_A_D2", "GABAMap_A_D1"},
            {"CX", "CX", "NMDA_D1", "NMDAMap_D1"}};

        const int rules_all = sizeof(rule_all) / sizeof(rule_all[0]);
        const int rules_paired = sizeof(rule_paired) / sizeof(rule_paired[0]);

        std::cerr << "Loading network\n";
        CellTypes& CT = net_cfg.CT;
        Connections& CN = net_cfg.CN;
        // cleanup x,y header information, we will fill it with what we read
        // from here.
        for (int ii = 0; ii < CT.size(); ii++)
            CT[ii].x = CT[ii].y = 0;

        FILE* fptr = fopen(file, "r");
        assert(fptr);
        // now we go through the connections file and fill in the info we need
        while (1) {
            int type, x_loc, y_loc;
            if (fscanf(fptr,
                       "Incoming connections for cell: type: %d x: %d y: %d\n",
                       &type, &x_loc, &y_loc) != 3)
                break;

            CT[type].x = std::max(x_loc, CT[type].x);
            CT[type].y = std::max(y_loc, CT[type].y);

            // second pass going through all the connections to the cell
            while (1) {
                char* syn_type = new char[256];
                // double strength = 0.0;
                // double mini_s = 0.0;
                // double mini_fre = 0.0;
                int mri_range = 0;
                int in_type, in_x_loc, in_y_loc;

                if (fscanf(fptr,
                           "type: %d x: %d y: %d Syntype: %s range: %d \n",
                           &in_type, &in_x_loc, &in_y_loc, syn_type,
                           &mri_range) != 5)
                    break;

                // transform synaptic type from given file to hybrid model with
                // MAP synapse type currently used.
                std::string syn_type_lit(syn_type);
                std::string from(CT[in_type].name);
                std::string to(CT[type].name);

                // Short/long range connection conversion:
                // 3 used for inter-hemispheric "short" range connections, 0 for
                // intrahemispheric long range connections, 1 forshort range 3
                // is currently ignored, we don't have connections defined
                if (mri_range == 3)
                    continue;
                int short_range = mri_range;  // except for 3 the same meaning
                                              // as in normal generation

                // transform synapse names according to the tables defined above
                for (int i = 0; i < rules_all; i++) {
                    if (rule_all[i][2] != syn_type_lit)
                        continue;
                    if (StringUtil::match(from, rule_all[i][0], to,
                                          rule_all[i][1]))
                        syn_type_lit = rule_all[i][3];
                }
                for (int i = 0; i < rules_paired; i++) {
                    if (rule_paired[i][2] != syn_type_lit)
                        continue;
                    if (StringUtil::match_paired(from, rule_paired[i][0], to,
                                                 rule_paired[i][1]))
                        syn_type_lit = rule_paired[i][3];
                }

                int conn_type =
                    CN.find(in_type, type, syn_type_lit, short_range);
                add_edge(in_type, in_x_loc, in_y_loc, type, x_loc, y_loc,
                         &CN.C[conn_type], 0);
                CN.stats[conn_type]++;  // for stats printing
            }  // while 2nd pass
        }  // while file

        CN.dump_stats(CT, net_cfg.conn_summary_fstream);
        fclose(fptr);
    }  // load_MRI_network

    void dump_summary(std::ostream& ConnSummaryFile,
                      const std::vector<int>& nedges) {
        //
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        //
        for (int c = 0; c < CN.n; c++) {  // layer connections
            // dbg
            ConnSummaryFile << c + 1 << " of " << CN.n
                            << " , edges: " << nedges[c] << " \ttype: (";
            CT.dump(CN.C[c].from, ConnSummaryFile);
            ConnSummaryFile << ") -> (";
            CT.dump(CN.C[c].to, ConnSummaryFile);
            ConnSummaryFile << ") ";
            CN.dump(c, ConnSummaryFile);
        }  // c
    }

    //
    void generate_3D_intrahemispheric_connections() {
        //
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        // auto& Subnet = sd.SubnetRef();
        const Vector1d<int>& map_LH_RH = src_data.map_LH_RH_Ref();
        const Vector1d<int>& map_LH_RH_small = src_data.map_LH_RH_small_Ref();
        //
        // generates connection using external distance matrix
        std::cerr << "Generating intrahemispheric connections from pre-set 3D "
                     "MRI data "
                     "\n";

        int smallLH = CT[CT.get_type("TC")].x;
        int largeLH = CT[CT.get_type("CX")].x;
        int smallRH = CT[CT.get_type("rTC")].x;
        int largeRH = CT[CT.get_type("rCX")].x;

        std::cerr << "Full set(CX): " << largeLH << " (L), " << largeRH
                  << " (R), Subnet(TC,RE,IN...): " << smallLH << " (L), "
                  << smallRH << " (R)\n";

        std::vector<int> nedges(CN.n, 0);  // edges counter for statistics
        for (int c = 0; c < CN.n; c++) {   // layer connections
            if (CT.right_hemisphere(CN.C[c].from) ==
                CT.right_hemisphere(CN.C[c].to))
                continue;  // only connections between different hemispheres
            int FX = CT[CN.C[c].from].x;
            int FY = CT[CN.C[c].from].y;  // input layer dimension
            int TX = CT[CN.C[c].to].x;
            int TY = CT[CN.C[c].to].y;  // output layer dimension
            assert(FY == 1);
            assert(TY == 1);  // MRI data are stored as 1D indexes
            assert(FX == largeLH);
            assert(TX == largeRH);  // at this moment we allow only CX-rCX rules
            // HACK* before we get proper fine mesh file from Eran, at this
            // moment we switch from fine mesh to coarse mesh //TODO remove hack
            FX = smallLH;
            TX = smallRH;

            for (int fx = 0; fx < FX; fx++) {
                int tx = map_LH_RH[fx];  // fx-tx pair of homotopic neurons

                // HACK* use coarse mesh and remap it to fine mesh //TODO remove
                // hack
                tx = map_LH_RH_small[fx];
                int fx2 = src_data.Subnet_Ref(0)[fx];
                int tx2 = src_data.Subnet_Ref(1)[tx];

                // printf("1 c: %d,fx2: %d,tx2: %d tx: %d \n", c,fx2,tx2, tx);
                nedges[c] += connect_local_neighbours(c, fx2, tx2, 1 /*right*/);
                // printf("0 c: %d,fx2: %d,tx2: %d\n", c,fx2,tx2);
                nedges[c] += connect_local_neighbours(c, tx2, fx2, 0 /*left*/);
            }  // fx
        }  // c

        dump_summary(net_cfg.conn_summary_fstream, nedges);
    }  // generate_3D_intrahemi_connections

    // for given pair fx,tx of homotopic neurons create neighbourhood around tx
    // and connect it to fx (works for both hemispheres)
    int connect_local_neighbours(int c /*rule*/, int fx, int tx,
                                 bool neighbourhood_hemisphere) {
        //
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        const Vector1d<int>& map_LH_RH = src_data.map_LH_RH_Ref();
        int ec = 0;  // edges counter
        for (int nbtx = 0; map_LH_RH[nbtx] != -1;
             nbtx++) {  // should be the same as largeLH
            // printf("fx: %d, tx: %d, nbtx: %d, right_hem:%d ", fx, tx,nbtx,
            // neighbourhood_hemisphere);
            double d = src_data.dist3D(tx, nbtx, neighbourhood_hemisphere);
            if (d >= CN.C[c].radius_min && d <= CN.C[c].radius_max) {
                assert(CN.C[c].range == 1);  // no long range done here
                if (MathUtil::rand01() <= CN.C[c].probab_oc) {
                    int from =
                        neighbourhood_hemisphere ? CN.C[c].from : CN.C[c].to;
                    int to =
                        neighbourhood_hemisphere ? CN.C[c].to : CN.C[c].from;
                    double dist =
                        src_data.dist3D(from, to, neighbourhood_hemisphere);
                    double delayTime = src_data.getScaledDelay(d);
                    check(CT, from, fx) && check(CT, to, nbtx) &&
                        (add_edge(from, fx, 0, to, nbtx, 0, &CN.C[c],
                                  delayTime),
                         ec++);  // ty==fy==0, see assert in caller
                }
            }  // radius
        }  // for nbtx
        return ec;
    }  // connect_local_neighbours

    // assign connection attributes for each edge in network (this can be
    // dynamic) this routine is separated from topology generation since we load
    // some topology-only networks
    void apply_synaptic_types() {

        std::cerr << "Sorting connections\n";
        std::sort(edges.begin(), edges.end(), Edge::cmp);

        std::cerr << "Adjusting synaptic types\n";
        std::vector<Edge>::iterator it = edges.begin();
        std::vector<Edge>::iterator end = edges.end();
        for (; it != end; it++) {
            Edge& c = (*it);
            if (c.syntype->distribution == "fixed")
                continue;  // syntype->already contains values we will use in
                           // output

            // TODO(x) decide what to do with gaussian or normal distribution
            //    if (c.syntype->distribution=="uniform") {}// ..
            //    if (c.syntype->distribution=="gauss") {}// ..
        }
    }  // apply

    // this is rough version, performance can be substantially improved in case
    // we still need it by
    // 1) testing/enumerating neurons only in radius - DONE
    // 2) Using 3D fixed array for sorting results instead of 1D vector
    //// generates geometry topology
    void generate_random_connections(bool useEuclideanDelays) {
        std::cerr << "No input file given, generating network on our own\n";
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        std::vector<int> nedges(CN.n, 0);  // edges counter for statistics
        for (int c = 0; c < CN.n; c++) {   // layer connections
            int FX = CT[CN.C[c].from].x;
            int FY = CT[CN.C[c].from].y;  // input layer dimension
            int TX = CT[CN.C[c].to].x;
            int TY = CT[CN.C[c].to].y;  // output layer dimension
            for (int fx = 0; fx < FX; fx++)
                for (int fy = 0; fy < FY; fy++) {
                    // get scaled position of target neuron in TO layer
                    int sc_x = static_cast<double>(TX) / FX * fx;
                    int sc_y = static_cast<double>(TY) / FY * fy;
                    // restrict the local neighbourhood where we look for
                    // possible connections by given radius. this speedup trick
                    // can fail in case someone switch to non-eucledian metric.
                    //[max(sc_x-radius,0)..min(sc_x+radius,TX)]
                    //[max(sc_y-radius,0)..min(sc_y+radius,TY)]
                    int radius = CN.C[c].radius_max + 1;

                    for (int tx = std::max(sc_x - radius, 0);
                         tx < std::min(sc_x + radius, TX);
                         tx++)  // particular connections between neurons
                        for (int ty = std::max(sc_y - radius, 0);
                             ty < std::min(sc_y + radius, TY); ty++) {

                            if (sc_x == tx && sc_y == ty &&
                                CT[CN.C[c].from].name == CT[CN.C[c].to].name)
                                continue;  // never connect neuron to itself
                            // distance is measured between 'From' neuron
                            // projected into 'TO' layer
                            double d = MathUtil::dist(sc_x, sc_y, tx, ty);

                            if (d >= CN.C[c].radius_min &&
                                d <= CN.C[c].radius_max) {
                                double delay = 0;
                                if (useEuclideanDelays) {
                                    delay = MathUtil::dist(
                                                static_cast<double>(fx),
                                                static_cast<double>(fy),
                                                static_cast<double>(tx),
                                                static_cast<double>(ty)) *
                                            0.5;
                                    // 2 millisecond delay for every
                                    // index the neurons are apart
                                    // delay = d; //this distance doesn't
                                    // correspond to he indexes used to create
                                    // the synapse (fx, fy, tx, ty) but it may
                                    // be better to use in other circumstances
                                    //            // fx != sc_x and fy != sc_y
                                    //            in all cases
                                }
                                // check if neurons are in same or different
                                // column
                                if (CN.C[c].range == 1) {
                                    if (src_data.is_column(sc_x, sc_y, tx, ty,
                                                           CN.C[c].cs)) {
                                        if (MathUtil::rand01() <=
                                            CN.C[c].probab)
                                            add_edge(CN.C[c].from, fx, fy,
                                                     CN.C[c].to, tx, ty,
                                                     &(CN.C[c]), delay),
                                                nedges[c]++;
                                    } else {  // out of column
                                        if (MathUtil::rand01() <=
                                            CN.C[c].probab_oc)
                                            add_edge(CN.C[c].from, fx, fy,
                                                     CN.C[c].to, tx, ty,
                                                     &CN.C[c], delay),
                                                nedges[c]++;
                                    }
                                } else {  // longrange
                                    if (MathUtil::rand01() <= CN.C[c].probab)
                                        add_edge(CN.C[c].from, fx, fy,
                                                 CN.C[c].to, tx, ty, &CN.C[c],
                                                 delay),
                                            nedges[c]++;
                                }
                            }  // radius
                        }  // ty
                }  // fy
        }  // c

        dump_summary(net_cfg.conn_summary_fstream, nedges);
    }  // generate_connections

  public:
    EdgeList(NetConfig& nc, const SourceData& sd, bool parc)  // NOLINT
        : net_cfg(nc), src_data(sd), parallel_construction(parc) {}

    void clear() { edges.clear(); }

    bool check(const CellTypes& CT, int ctype, int neuronx) const {
        return CT[ctype].x > neuronx;
    }

    void generate_radom_edges(bool useEuclideanDelays) {
        generate_random_connections(
            useEuclideanDelays);  // our own random connectivity
        apply_synaptic_types();   // needed in case we need sample
                                  // parameters for strength/mini_X
                                  // parameters
    }

    int generate_multilayer_edges() {
        //
        timer gen_timer;
        // we do this independently for distinct hemipheres, because of memory
        //(some crazy radius scenarios get up to 200 GB of RAM)
        gen_timer.reset();
        generate_multilayer_connections(0);
        gen_timer.measure_accumulate_print(
            "Generate Connections (Left) : ", std::cout, false);
        print_memory_usage<uint64_t>("; ", std::cout);
        //  CE.apply_synaptic_types();
        //  CE.write_connections(0);
        //  CE.edges.clear();

        //
        gen_timer.reset();
        generate_multilayer_connections(1);
        gen_timer.measure_accumulate_print(
            "Generate Connections (Right) : ", std::cout, false);
        print_memory_usage<uint64_t>("; ", std::cout);
        //  CE.apply_synaptic_types();
        //  CE.write_connections(1);
        //  CE.edges.clear();

        // intrahemispheric connections
        //    CE.generate_3D_intrahemispheric_connections();
        gen_timer.reset();
        apply_synaptic_types();
        gen_timer.measure_accumulate_print("Apply Synaptic Types : ", std::cout,
                                           false);
        print_memory_usage<uint64_t>(";", std::cout);
        return 0;
    }

    // Generate Outputs
    void write_connections(std::ostream& out_stream,
                           bool right_hemisphere = false) const {
        std::cerr << "Writing output file\n";
        // Cell Types
        if (!right_hemisphere)
            net_cfg.CT.output(out_stream);
        int te = 0, tn = 0;  // total edges, total neurons
        // Edges
        std::vector<Edge>::const_iterator it = edges.cbegin(),
                                          end = edges.cend(),
                                          prev = edges.cbegin();
        for (; it != end; prev = it, it++, te++) {
            assert(it->syntype != nullptr);
            // Edge type header
            if (it == edges.cbegin() || it->same_to(*prev) == false) {
                it->header(out_stream);
                tn++;
            }
            it->output(out_stream);
            //(*it).dump();
        }
        std::cerr << "Total counts (hemisphere " << right_hemisphere
                  << "): " << tn << " neurons, " << te << " edges.\n";
    }  // write_connections

    // dump all output connection for a given neuron (we usually write down the
    // opposite)
    void dump_from_neuron_edges(int type, int x, int y, int to_type) {
        std::vector<Edge>::iterator it = edges.begin();
        std::vector<Edge>::iterator end = edges.end();
        for (; it != end; it++) {
            Edge& c = (*it);
            if (c.from_t == type && c.from_x == x && c.from_y == y &&
                (to_type == -1 || c.to_t == to_type))
                std::cerr << type << " " << x << "," << y << " -> " << c.to_t
                          << " " << c.to_x << "," << c.to_y << "\n";
        }
    }

    // dump
    void dump_from_neuron_edges(int type, int x, int y) {
        dump_from_neuron_edges(type, x, y, -1);
    }
};  // Edges

#endif  // EDGE_LIST_HPP
