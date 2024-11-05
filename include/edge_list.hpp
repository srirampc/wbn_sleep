#ifndef EDGE_LIST_HPP
#define EDGE_LIST_HPP

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <cell_types.hpp>
#include <connections.hpp>
#include <edge.hpp>
#include <net_config.hpp>
#include <source_data.hpp>

class EdgeList {
  public:
    const Edge NullEdge = Edge(-1, -1, -1, -1, -1, -1, nullptr, 0.0);

  private:
    NetConfig& net_cfg;
    const SourceData& src_data;
    std::vector<Edge> edges;

    void generate_3D_connections(bool right_hemisphere) {
        //
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        auto& ConnSummaryFile = net_cfg.conn_summary_fstream;
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

        for (int c = 0; c < CN.n; c++) {  // layer connections
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

            int ec = 0;  // edges counter for statistics

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
                                    (edges.push_back(Edge(
                                         CN.C[c].from, (FX == large) ? fx : fxi,
                                         fy, CN.C[c].to,
                                         (TX == large) ? tx : txi, ty, &CN.C[c],
                                         delayTime)),
                                     ec++);
                        } else {  // longrange
                            if (MathUtil::rand01() <= CN.C[c].probab)
                                check(CT, CN.C[c].from,
                                      (FX == large) ? fx : fxi) &&
                                    check(CT, CN.C[c].to,
                                          (TX == large) ? tx : txi) &&
                                    (edges.push_back(Edge(
                                         CN.C[c].from, (FX == large) ? fx : fxi,
                                         fy, CN.C[c].to,
                                         (TX == large) ? tx : txi, ty, &CN.C[c],
                                         delayTime)),
                                     ec++);
                        }  // longrange
                    }  // radius
                }  // tx
            }  // fx

            // dbg
            ConnSummaryFile << c + 1 << " of " << CN.n << " , edges: " << ec
                            << " \ttype: (";
            CT.dump(CN.C[c].from, ConnSummaryFile);
            ConnSummaryFile << ") -> (";
            CT.dump(CN.C[c].to, ConnSummaryFile);
            ConnSummaryFile << ") ";
            CN.dump(c, ConnSummaryFile);
        }  // c
    }

    void generate_MultiLayer_connections(bool right_hemisphere) {
        //
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        auto& ConnSummaryFile = net_cfg.conn_summary_fstream;
        const Vector1d<float>& Subnet_hemisphere =
            src_data.Subnet_hemisphere_Ref(right_hemisphere);
        const std::map<std::string, int>& mapLayerNameId_IN =
            src_data.mapLayerNameId_IN_Ref();
        const std::map<std::string, int>& mapLayerNameId_PY =
            src_data.mapLayerNameId_PY_Ref();
        const Vector4d<float>& weight_factor = src_data.weight_factor_Ref();
        const Vector2d<float>& weight_factor_INPY =
            src_data.weight_factor_INPY_Ref();
        const Vector1d<float>& TC_Dist = src_data.TC_Dist_Ref();
        float maxDelay_TCm = src_data.maxDelay_TCm;
        float maxDelay_TCc = src_data.maxDelay_TCc;
        bool both_hemispheres_detected = net_cfg.both_hemispheres_detected;
        //
        std::cerr << "Generating connection file from pre-set 3D MRI data, new "
                     "MultiLayer "
                     "connections, hemisphere: "
                  << right_hemisphere << "\n";

        // // uncomment to use a non-deterministic seed:
        //    std::random_device rd;
        //    std::mt19937 gen(rd());
        std::mt19937 gen(1937);

        int small = CT[CT.get_type("TC")].x;
        int large = CT[CT.get_type("CX")].x;
        int indOffset = 0;
        if (both_hemispheres_detected && right_hemisphere) {
            small = CT[CT.get_type("rTC")].x, large = CT[CT.get_type("rCX")].x;
            indOffset = CT[CT.get_type("CX")].x;
        }

        std::cerr << "Full set(CX): " << large
                  << ", Subnet(TC,RE,IN...): " << small << "\n";

        for (int c = 0; c < CN.n; c++) {  // layer connections
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

            std::string fromCellType = CT[CN.C[c].from].name;
            std::string toCellType = CT[CN.C[c].to].name;

            int fromLayer = -1;
            int toLayer = -1;
            int from_IN = 0;
            int to_IN = 0;
            if (mapLayerNameId_PY.find(CT[CN.C[c].from].name) !=
                mapLayerNameId_PY.end()) {
                fromLayer =
                    mapLayerNameId_PY.find(CT[CN.C[c].from].name)->second;
            } else if (mapLayerNameId_IN.find(CT[CN.C[c].from].name) !=
                       mapLayerNameId_IN.end()) {
                from_IN = 1;
            }

            if (mapLayerNameId_PY.find(CT[CN.C[c].to].name) !=
                mapLayerNameId_PY.end()) {
                toLayer = mapLayerNameId_PY.find(CT[CN.C[c].to].name)->second;
            } else if (mapLayerNameId_IN.find(CT[CN.C[c].to].name) !=
                       mapLayerNameId_IN.end()) {
                to_IN = 1;
            }

            assert(fromLayer < 6);
            assert(toLayer < 6);

            int ec = 0;  // edges counter for statistics

            for (int fxi = 0, fx = (FX == large) ? 0 : Subnet_hemisphere[0];
                 (FX == large) ? fx < FX : Subnet_hemisphere[fxi] != -1;
                 (FX == large) ? fx++ : fx = Subnet_hemisphere[++fxi]) {
                // scaling is actually hidden in Subnet array, identity in large
                // case, map in small case
                int fy = 0;
                int sc_x = fx;
                int sc_y = fy;

                for (int txi = 0,
                         tx = (TX == large) ? 0 : Subnet_hemisphere[0];
                     (TX == large) ? tx < TX : Subnet_hemisphere[txi] != -1;
                     (TX == large)
                         ? tx++
                         : tx = Subnet_hemisphere[++txi]) {  // particular
                                                             // connections
                                                             // between neurons
                    int ty = 0;
                    if (sc_x == tx && sc_y == ty &&
                        CT[CN.C[c].from].name == CT[CN.C[c].to].name)
                        continue;  // never connect neuron to itself
                    // distance is measured between 'From' neuron projected into
                    // 'TO' layer
                    //  printf("sc_x: %d, tx: %d, right_hem:%d ",
                    //  sc_x,tx,right_hemisphere);
                    double d =
                        src_data.dist_only(sc_x + indOffset, tx + indOffset);
                    // if(d<0) continue;
                    double delayTime = 0.0;  // = getScaledDelay(d);

                    // use spherical Thalamus
                    size_t foundFromTC = fromCellType.find("TC");
                    size_t foundFromRE = fromCellType.find("RE");
                    size_t foundToTC = toCellType.find("TC");
                    size_t foundToRE = toCellType.find("RE");

                    if (((foundFromTC != std::string::npos) ||
                         (foundFromRE != std::string::npos)) &&
                        ((foundToTC != std::string::npos) ||
                         (foundToRE != std::string::npos))) {
                        d = src_data.intraThalamicDistOnSphere(sc_x, tx);
                        // no indOffsetfor the right hemi within thalamus
                    }

                    if (d < 0)
                        continue;

                    // use dealyes only on CX->CX and CX->IN connections
                    size_t foundFromCX = fromCellType.find("CX");
                    size_t foundToCX = toCellType.find("CX");
                    size_t foundToIN = toCellType.find("IN");
                    if (!((foundFromCX != std::string::npos) &&
                          ((foundToCX != std::string::npos) ||
                           (foundToIN != std::string::npos)))) {
                        delayTime = 0.0;
                    } else {
                        delayTime = src_data.getScaledDelay(d);
                    }

                    double currEdgeProb =
                        src_data.prob_only(sc_x + indOffset, tx + indOffset);
                    double clampProb = (CN.C[c].probab_oc * currEdgeProb < 1)
                                           ? CN.C[c].probab_oc * currEdgeProb
                                           : 1;

                    if (fromLayer > -1 && toLayer > -1) {  // PY->PY connections
                        std::bernoulli_distribution distribution(clampProb);
                        if (distribution(gen)) {
                            double currEdgeWeightFactor =
                                weight_factor(sc_x + indOffset, tx + indOffset,
                                              fromLayer, toLayer);
                            if (d <= 0.0001) {
                                currEdgeWeightFactor *=
                                    5;  // local weights scaled up to test their
                                        // influence
                            } else if (d > CN.C[c].radius_max) {
                                currEdgeWeightFactor = 0;
                            }
                            if (currEdgeWeightFactor > 0) {
                                check(CT, CN.C[c].from,
                                      (FX == large) ? fx : fxi) &&
                                    check(CT, CN.C[c].to,
                                          (TX == large) ? tx : txi) &&
                                    (edges.push_back(Edge(
                                         CN.C[c].from, (FX == large) ? fx : fxi,
                                         fy, CN.C[c].to,
                                         (TX == large) ? tx : txi, ty, &CN.C[c],
                                         delayTime, currEdgeWeightFactor)),
                                     ec++);
                            }
                        }
                    } else if (fromLayer > -1 &&
                               to_IN > 0) {  // PY->IN connections
                        std::bernoulli_distribution distribution(clampProb);
                        if (distribution(gen)) {
                            double currEdgeWeightFactor = weight_factor(
                                sc_x + indOffset, tx + indOffset, fromLayer,
                                fromLayer);  // Since PY->IN are only within
                                             // the same layer, use 'fromLayer'
                                             // twice
                            if (currEdgeWeightFactor > 0) {
                                check(CT, CN.C[c].from,
                                      (FX == large) ? fx : fxi) &&
                                    check(CT, CN.C[c].to,
                                          (TX == large) ? tx : txi) &&
                                    (edges.push_back(Edge(
                                         CN.C[c].from, (FX == large) ? fx : fxi,
                                         fy, CN.C[c].to,
                                         (TX == large) ? tx : txi, ty, &CN.C[c],
                                         delayTime, currEdgeWeightFactor)),
                                     ec++);
                            }
                        }
                    } else if (from_IN > 0 &&
                               toLayer > -1) {  // IN->PY connections
                        if (d >= CN.C[c].radius_min &&
                            d <= CN.C[c].radius_max) {
                            if (CN.C[c].range == 1) {
                                if (MathUtil::rand01() <= CN.C[c].probab_oc) {
                                    double currEdgeWeightFactor =
                                        weight_factor_INPY(sc_x + indOffset,
                                                           tx + indOffset);
                                    if (currEdgeWeightFactor > 0) {
                                        check(CT, CN.C[c].from,
                                              (FX == large) ? fx : fxi) &&
                                            check(CT, CN.C[c].to,
                                                  (TX == large) ? tx : txi) &&
                                            (edges.push_back(
                                                 Edge(CN.C[c].from,
                                                      (FX == large) ? fx : fxi,
                                                      fy, CN.C[c].to,
                                                      (TX == large) ? tx : txi,
                                                      ty, &CN.C[c], delayTime,
                                                      currEdgeWeightFactor)),
                                             ec++);
                                    }
                                }
                            }
                        }  // radius
                    } else {  // all other types of connections, all are based
                              // on the disc model
                        if (d >= CN.C[c].radius_min &&
                            d <= CN.C[c].radius_max) {
                            // use TC_Dist only on thalamocortical connections
                            // (CX,TC,IN,RE)
                            size_t foundFromCX = fromCellType.find("CX");
                            size_t foundFromIN = fromCellType.find("IN");
                            size_t foundToCX = toCellType.find("CX");
                            size_t foundToIN = toCellType.find("IN");

                            size_t foundFromTC = fromCellType.find("TC");
                            size_t foundFromRE = fromCellType.find("RE");
                            size_t foundToTC = toCellType.find("TC");
                            size_t foundToRE = toCellType.find("RE");

                            // different maxDealy on Core and Matrix connections
                            // be careful with condition since there is CX5a
                            // which is not matrix (TCa)
                            size_t foundFromMatrix = fromCellType.find("a");
                            size_t foundToMatrix = toCellType.find("a");

                            if (((foundFromCX != std::string::npos) ||
                                 (foundFromIN != std::string::npos)) &&
                                ((foundToTC != std::string::npos) ||
                                 (foundToRE != std::string::npos))) {
                                if ((foundToMatrix != std::string::npos)) {
                                    // TODO(x) for both hemis use: sc_x +
                                    // indOffset, tx
                                    // + indOffset
                                    delayTime = src_data.getScaled_TC_Delay(
                                        TC_Dist[sc_x], maxDelay_TCm);
                                } else {
                                    delayTime = src_data.getScaled_TC_Delay(
                                        TC_Dist[sc_x], maxDelay_TCc);
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
                                    delayTime = src_data.getScaled_TC_Delay(
                                        TC_Dist[tx], maxDelay_TCm);
                                } else {
                                    delayTime = src_data.getScaled_TC_Delay(
                                        TC_Dist[tx], maxDelay_TCc);
                                }
                            }

                            if (delayTime < 0) {
                                delayTime = 0;
                            }

                            double currEdgeWeightFactor = 1;
                            if ((CT[CN.C[c].from].name == "TC") &&
                                (CT[CN.C[c].to].name == "CX4")) {
                                currEdgeWeightFactor = 3;
                            }
                            // column is no more checked
                            if (CN.C[c].range == 1) {
                                if (MathUtil::rand01() <= CN.C[c].probab_oc)
                                    check(CT, CN.C[c].from,
                                          (FX == large) ? fx : fxi) &&
                                        check(CT, CN.C[c].to,
                                              (TX == large) ? tx : txi) &&
                                        (edges.push_back(
                                             Edge(CN.C[c].from,
                                                  (FX == large) ? fx : fxi, fy,
                                                  CN.C[c].to,
                                                  (TX == large) ? tx : txi, ty,
                                                  &CN.C[c], delayTime,
                                                  currEdgeWeightFactor)),
                                         ec++);
                            } else {  // longrange
                                if (MathUtil::rand01() <= CN.C[c].probab)
                                    check(CT, CN.C[c].from,
                                          (FX == large) ? fx : fxi) &&
                                        check(CT, CN.C[c].to,
                                              (TX == large) ? tx : txi) &&
                                        (edges.push_back(
                                             Edge(CN.C[c].from,
                                                  (FX == large) ? fx : fxi, fy,
                                                  CN.C[c].to,
                                                  (TX == large) ? tx : txi, ty,
                                                  &CN.C[c], delayTime,
                                                  currEdgeWeightFactor)),
                                         ec++);
                            }  // longrange
                        }  // radius
                    }

                }  // tx
            }  // fx

            // dbg
            ConnSummaryFile << c + 1 << " of " << CN.n << " , edges: " << ec
                            << " \ttype: (";
            CT.dump(CN.C[c].from, ConnSummaryFile);
            ConnSummaryFile << ") -> (";
            CT.dump(CN.C[c].to, ConnSummaryFile);
            ConnSummaryFile << ") ";
            CN.dump(c, ConnSummaryFile);
        }  // c
    }  // generate_MultiLayer_connections

    // load network from the file (fixed structure) we
    // get from ucsd folks approximating MRI topology
    void load_MRI_network(const char* file) {
        //
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
        auto& ConnSummaryFile = net_cfg.conn_summary_fstream;
        // cleanup x,y header information, we will fill it with what we read
        // from here.
        for (int ii = 0; ii < CT.size(); ii++)
            CT[ii].x = CT[ii].y = 0;

        FILE* f = fopen(file, "r");
        assert(f);
        // now we go through the connections file and fill in the info we need
        while (1) {
            int type, x_loc, y_loc;
            if (fscanf(f,
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

                if (fscanf(f, "type: %d x: %d y: %d Syntype: %s range: %d \n",
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
                edges.push_back(Edge(in_type, in_x_loc, in_y_loc, type, x_loc,
                                     y_loc, &CN.C[conn_type], 0));
                CN.stats[conn_type]++;  // for stats printing
            }  // while 2nd pass
        }  // while file

        CN.dump_stats(CT, ConnSummaryFile);
        fclose(f);
    }  // load_MRI_network

    //
    void generate_3D_intrahemispheric_connections() {
        //
        const CellTypes& CT = net_cfg.CT;
        const Connections& CN = net_cfg.CN;
        auto& ConnSummaryFile = net_cfg.conn_summary_fstream;
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

        for (int c = 0; c < CN.n; c++) {  // layer connections
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

            int ec = 0;  // edges counter for statistics

            // HACK* before we get proper fine mesh file from Eran, at this
            // moment we switch from fine mesh to coarse mesh //TODO remove hack
            FX = smallLH;
            TX = smallRH;

            for (int fx = 0; fx < FX; fx++) {
                int tx = map_LH_RH[fx];  // fx-tx pair of homotopic neurons

                // HACK* use coarse mesh and remap it to fine mesh //TODO remove
                // hack
                tx = map_LH_RH_small[fx];
                int fx2 = src_data.Subnet_hemisphere_Ref(0)[fx];
                int tx2 = src_data.Subnet_hemisphere_Ref(1)[tx];

                // printf("1 c: %d,fx2: %d,tx2: %d tx: %d \n", c,fx2,tx2, tx);
                ec += connect_local_neighbours(c, fx2, tx2, 1 /*right*/);
                // printf("0 c: %d,fx2: %d,tx2: %d\n", c,fx2,tx2);
                ec += connect_local_neighbours(c, tx2, fx2, 0 /*left*/);
            }  // fx

            // dbg
            ConnSummaryFile << c + 1 << " of " << CN.n << " , edges: " << ec
                            << " \ttype: (";
            CT.dump(CN.C[c].from, ConnSummaryFile);
            ConnSummaryFile << ") -> (";
            CT.dump(CN.C[c].to, ConnSummaryFile);
            ConnSummaryFile << ") ";
            CN.dump(c, ConnSummaryFile);
        }  // c
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
                        (edges.push_back(Edge(from, fx, 0, to, nbtx, 0,
                                              &CN.C[c], delayTime)),
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
        auto& ConnSummaryFile = net_cfg.conn_summary_fstream;
        for (int c = 0; c < CN.n; c++) {  // layer connections
            int FX = CT[CN.C[c].from].x;
            int FY = CT[CN.C[c].from].y;  // input layer dimension
            int TX = CT[CN.C[c].to].x;
            int TY = CT[CN.C[c].to].y;  // output layer dimension
            int ec = 0;                 // edges counter for statistics
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
                                            edges.push_back(
                                                Edge(CN.C[c].from, fx, fy,
                                                     CN.C[c].to, tx, ty,
                                                     &(CN.C[c]), delay)),
                                                ec++;
                                    } else {  // out of column
                                        if (MathUtil::rand01() <=
                                            CN.C[c].probab_oc)
                                            edges.push_back(
                                                Edge(CN.C[c].from, fx, fy,
                                                     CN.C[c].to, tx, ty,
                                                     &CN.C[c], delay)),
                                                ec++;
                                    }
                                } else {  // longrange
                                    if (MathUtil::rand01() <= CN.C[c].probab)
                                        edges.push_back(Edge(
                                            CN.C[c].from, fx, fy, CN.C[c].to,
                                            tx, ty, &CN.C[c], delay)),
                                            ec++;
                                }
                            }  // radius
                        }  // ty
                }  // fy

            // dbg
            ConnSummaryFile << c + 1 << " of " << CN.n << " , edges: " << ec
                            << " \ttype: (";
            CT.dump(CN.C[c].from, ConnSummaryFile);
            ConnSummaryFile << ") -> (";
            CT.dump(CN.C[c].to, ConnSummaryFile);
            ConnSummaryFile << ") ";
            CN.dump(c, ConnSummaryFile);
        }  // c
    }  // generate_connections

  public:
    EdgeList(NetConfig& nc, const SourceData& sd)  // NOLINT
        : net_cfg(nc), src_data(sd) {}

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
        //
        // we do this independently for distinct hemipheres, because of memory
        //(some crazy radius scenarios get up to 200 GB of RAM)
        gen_timer.reset();
        generate_MultiLayer_connections(0);
        gen_timer.measure_accumulate_print(
            "Generate Connections (Left) : ", std::cout, false);
        print_memory_usage<uint64_t>("; ", std::cout);
        //  CE.apply_synaptic_types();
        //  CE.write_connections(0);
        //  CE.edges.clear();

        //
        gen_timer.reset();
        generate_MultiLayer_connections(1);
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
    }  // dump
    void dump_from_neuron_edges(int type, int x, int y) {
        dump_from_neuron_edges(type, x, y, -1);
    }
};  // Edges

#endif  // EDGE_LIST_HPP
