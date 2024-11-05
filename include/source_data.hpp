#ifndef SOURCE_DATA_HPP
#define SOURCE_DATA_HPP

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
//
#include <nlohmann/json.hpp>
using njson = nlohmann::json;
#include <mio/mio.hpp>
//
#include <cell_types.hpp>
#include <connections.hpp>
#include <data_config.hpp>
#include <net_config.hpp>
#include <utils.hpp>
#include <vector_nd.hpp>
//
using BYTE = unsigned char;
//

class SourceData {
  public:
    static constexpr int MRI_POINTS = 12500;
    static constexpr int MAX_SUBNET_POINTS = 12500;
    static constexpr int NUM_TC_POINTS = 12500;
    static constexpr int MRI_POINTS_FULL = 20500;  // 20500 //  12500
    static constexpr int NUM_PARAM = 2;
    static constexpr int NUM_LAYERS = 6;
    //
    static constexpr int mriPointsFull = 20500;
    static constexpr int thalPopulations = 2 * 642;
    static constexpr int numParam = 2;
    static constexpr int numLayers = 6;

    // max delay between farthest neurons. all delays will be scaled according
    // to this
    static constexpr float maxDelay = 40;
    static constexpr float maxDelay_TCc = 40;  // core TC
    static constexpr float maxDelay_TCm = 40;  // matrix TC
    static constexpr float avgDist = 0;
    static constexpr float numCons = 0;

  private:
    const NetConfig& nc;
    float maxDistance = 0;  // max distance between two neurons
    float max_TC_Distance = 0;
    bool useDelays = false;
    //
    // distance lookup table given for MRI data; [0][][] for left hemispher
    std::array<Vector2d<float>, 2> Distance3D;
    //
    // The actual number of smaller net is to be read out from network.cfg (TC
    // cells)
    // distance lookup table given for MRI data, [0][] for left hemisphere
    std::array<Vector1d<float>, 2> Subnet;
    //
    // mapping between LH-RH neurons
    Vector1d<int> map_LH_RH = Vector1d<int>(MAX_SUBNET_POINTS, 0);
    //
    // mapping between a subset of LH-RH neurons
    Vector1d<int> map_LH_RH_small = Vector1d<int>(MAX_SUBNET_POINTS, 0);
    //
    //
    Vector1d<float> TC_Dist = Vector1d<float>(NUM_TC_POINTS, 0.0);
    //
    // lookup table given for MRI data; [0][][] for left hemispher
    Vector3d<float> Dist_Prob;
    //
    //
    Vector2d<float> Dist_intraThalamus;
    //
    /// distance lookup table given for MRI data; [0][][] for left hemispher
    Vector4d<float> weight_factor;
    Vector2d<float> weight_factor_INPY;
    //
    std::map<std::string, int> mapLayerNameId_IN;
    std::map<std::string, int> mapLayerNameId_PY;

  public:
    explicit SourceData(const NetConfig& in_nc) : nc(in_nc) {}
    ~SourceData() {}

    // Access functions
    inline const std::array<Vector1d<float>, 2>& SubnetRef() const {
        return Subnet;
    }

    inline const Vector1d<float>&
    Subnet_hemisphere_Ref(bool right_hemisphere) const {
        return Subnet[right_hemisphere];
    }

    inline const Vector1d<int>& map_LH_RH_Ref() const { return map_LH_RH; }

    inline const Vector1d<int>& map_LH_RH_small_Ref() const {
        return map_LH_RH_small;
    }

    inline const std::map<std::string, int>& mapLayerNameId_IN_Ref() const {
        return mapLayerNameId_IN;
    }

    inline const std::map<std::string, int>& mapLayerNameId_PY_Ref() const {
        return mapLayerNameId_PY;
    }

    inline const Vector4d<float>& weight_factor_Ref() const {
        return weight_factor;
    }

    inline const Vector2d<float>& weight_factor_INPY_Ref() const {
        return weight_factor_INPY;
    }

    inline const Vector1d<float>& TC_Dist_Ref() const { return TC_Dist; }

    inline double dist3D(int id_from, int id_to, bool right_hemisphere) const {
        // printf("Distance3D - %lf from %d to %d
        // \n",Distance3D[right_hemisphere][id_from][id_to],id_from,id_to);

        assert(Distance3D[right_hemisphere](id_from, id_to) != -1);
        return Distance3D[right_hemisphere](id_from, id_to);
    }

    double getScaledDelay(double distance) const {
        double dist = useDelays
                          ? distance
                          : 0;  // if you should use delays or not use delays
        return (dist / maxDistance) * maxDelay;
    }

    double getScaled_TC_Delay(double distance, double maxDelay) const {
        double dist = useDelays
                          ? distance
                          : 0;  // if you should use delays or not use delays
        return (dist / max_TC_Distance) * maxDelay;
    }

    inline double dist_only(int id_from, int id_to) const {
        // assert(Dist_Prob[id_from][id_to][0]!=-1);
        return Dist_Prob(id_from, id_to, 0);
    }

    inline double prob_only(int id_from, int id_to) const {
        // assert(Dist_Prob[id_from][id_to][1]!=-1);
        return Dist_Prob(id_from, id_to, 1);
    }

    inline double intraThalamicDistOnSphere(int id_from, int id_to) const {
        return Dist_intraThalamus(id_from, id_to);
    }

    // Utility Functions
    static inline int is_column(int x1, int y1, int x2, int y2, int CS) {
        if ((x1 / CS == x2 / CS) && (y1 / CS == y2 / CS))
            return 1;
        else
            return 0;
    }

  private:
    int handle_error(const std::error_code& error) {
        const auto& errmsg = error.message();
        std::cerr << "error mapping file:  " << errmsg << " exiting..."
                  << std::endl;
        return error.value();
    }

    std::map<std::string, int> fillDict_PY(std::map<std::string, int> m) {
        m.insert(std::pair<std::string, int>("CX", 0));
        m.insert(std::pair<std::string, int>("CX3", 1));
        m.insert(std::pair<std::string, int>("CX4", 2));
        m.insert(std::pair<std::string, int>("CX5a", 3));
        m.insert(std::pair<std::string, int>("CX5b", 4));
        m.insert(std::pair<std::string, int>("CX6", 5));

        m.insert(std::pair<std::string, int>("rCX", 0));
        m.insert(std::pair<std::string, int>("rCX3", 1));
        m.insert(std::pair<std::string, int>("rCX4", 2));
        m.insert(std::pair<std::string, int>("rCX5a", 3));
        m.insert(std::pair<std::string, int>("rCX5b", 4));
        m.insert(std::pair<std::string, int>("rCX6", 5));

        return m;
    }

    std::map<std::string, int> fillDict_IN(std::map<std::string, int> m) {

        m.insert(std::pair<std::string, int>("IN", 0));
        m.insert(std::pair<std::string, int>("IN3", 1));
        m.insert(std::pair<std::string, int>("IN4", 2));
        m.insert(std::pair<std::string, int>("IN5a", 3));
        m.insert(std::pair<std::string, int>("IN5b", 4));
        m.insert(std::pair<std::string, int>("IN6", 5));

        m.insert(std::pair<std::string, int>("rIN", 0));
        m.insert(std::pair<std::string, int>("rIN3", 1));
        m.insert(std::pair<std::string, int>("rIN4", 2));
        m.insert(std::pair<std::string, int>("rIN5a", 3));
        m.insert(std::pair<std::string, int>("rIN5b", 4));
        m.insert(std::pair<std::string, int>("rIN6", 5));

        return m;
    }

    //
    // Data Loading Functions
    void load_3D_distances(const char* file, bool right_hemisphere) {
        if (Distance3D[right_hemisphere].size() == 0) {
            Distance3D[right_hemisphere] =
                Vector2d<float>(MRI_POINTS, MRI_POINTS, -1.0);
        }
        //
        std::cerr << "Loading 3D network distances (hemisphere): " << file
                  << "(" << (right_hemisphere ? "R" : "L") << ")" << std::endl;

        for (int i = 0; i < MRI_POINTS; i++)
            for (int ii = 0; ii < MRI_POINTS; ii++)
                Distance3D[right_hemisphere](i, ii) = -1;

        FILE* f = fopen(file, "r");
        assert(f);

        // go through the distance matrix
        // int i = 0;
        while (1) {
            int from, to;
            float dist;
            if (fscanf(f, "%d %d %f\n", &from, &to, &dist) != 3)
                break;
            from--;
            to--;  // MATLAB vs C indexing
            assert(from < MRI_POINTS);
            assert(to < MRI_POINTS);
            Distance3D[right_hemisphere](from, to) = dist;
            // avgDist += dist;
            // numCons = numCons + 1;

            if (dist > maxDistance)  // get max distance for scaling later
                maxDistance = dist;
        }  // while

        fclose(f);
    }

    // load subset of CX neurons to be flagged as indexes for TC,RE & IN
    void load_3D_subnet(const char* file, bool right_hemisphere) {
        std::cerr << "Loading 3D sub network (hemisphere): " << file << "("
                  << (right_hemisphere ? "R" : "L") << ")" << std::endl;
        //
        if (Subnet[right_hemisphere].size() == 0) {
            Subnet[right_hemisphere] = Vector1d<float>(MAX_SUBNET_POINTS, 0.0);
        }

        for (int i = 0; i < MAX_SUBNET_POINTS; i++)
            Subnet[right_hemisphere](i) = -1;

        FILE* f = fopen(file, "r");
        assert(f);

        // go through the distance matrix
        // int i = 0;
        while (1) {
            int id, map;
            if (fscanf(f, "%d %d\n", &id, &map) != 2)
                break;
            id--;
            map--;  // MATLAB vs C indexing
            //    cerr<<id<<" "<<map<<" " << CT.cell[CT.get_type("TC")].x<<"\n";
            if (right_hemisphere)
                assert(id < nc.CT[nc.CT.get_type("rTC")].x);
            if (!right_hemisphere)
                assert(id < nc.CT[nc.CT.get_type("TC")].x);
            assert(id < MRI_POINTS);
            assert(map < MAX_SUBNET_POINTS);

            Subnet[right_hemisphere](id) = map;
        }  // while

        fclose(f);
    }

    void load_3D_LH_RH_correspondence(const char* file, int Map[]) {
        std::cerr << "Loading 3D LH-RH interhemispheric mapping : " << file
                  << std::endl;

        for (int i = 0; i < MRI_POINTS; i++)
            Map[i] = -1;

        FILE* f = fopen(file, "r");
        assert(f);

        // int i = 0;
        while (1) {
            int id, map;
            if (fscanf(f, "%d %d\n", &id, &map) != 2)
                break;
            id--;
            map--;                    // MATLAB vs C indexing
            assert(id < MRI_POINTS);  // assert(map<MAX_SUBNET_POINTS);

            Map[id] = map;

            // printf("Mapping id: %d to %d \n",id,map);
        }  // while

        fclose(f);
    }

    void load_3D_LH_RH_correspondence(const char* data_file, bool full_data) {
        if (full_data) {
            load_3D_LH_RH_correspondence(data_file, map_LH_RH.data());
        } else {
            load_3D_LH_RH_correspondence(data_file, map_LH_RH_small.data());
        }
    }

    void load_thalCort_dist(const char* file) {
        std::cerr << "Load thalCort_dist : " << file << std::endl;

        for (int ii = 0; ii < NUM_TC_POINTS; ii++)
            TC_Dist(ii) = -1;

        FILE* f = fopen(file, "r");
        assert(f);

        int i = 0;
        while (1) {
            float dist;
            if (fscanf(f, "%f\n", &dist) != 1)
                break;

            if (dist >= 0) {
                TC_Dist(i) = dist;
            }

            if (dist > max_TC_Distance) {  // get max distance for scaling later
                max_TC_Distance = dist;
            }

            i++;
        }  // while

        fclose(f);
    }

    void load_3D_dist_prob(const char* file) {
        Dist_Prob = Vector3d<float>(MRI_POINTS, MRI_POINTS, NUM_PARAM, -1.0);

        std::cerr << "Load 3D_dist_prob : " << file << std::endl;

        // for (int i=0; i<MRI_POINTS_FULL; i++){
        //    for (int ii=0; ii<MRI_POINTS_FULL; ii++){
        //      for (int kk=0; kk<NUM_PARAM; kk++){
        //      // cerr << i << " " << ii << " " << kk << endl;
        //	      Dist_Prob[i][ii][kk]=-1;
        //      }
        //    }
        //  }

        FILE* f = fopen(file, "r");
        assert(f);

        // int i = 0;
        while (1) {
            int from, to;
            float dist, conn_prob;
            if (fscanf(f, "%d %d %f %f\n", &from, &to, &dist, &conn_prob) != 4)
                break;
            from--;
            to--;  // MATLAB vs C indexing
            assert(from < MRI_POINTS);
            assert(to < MRI_POINTS);
            Dist_Prob(from, to, 0) = dist;
            Dist_Prob(from, to, 1) = conn_prob;

            if (dist > maxDistance)  // get max distance for scaling later
                maxDistance = dist;
        }  // while

        fclose(f);
    }

    void load_intraThalamus_dist(const char* file) {
        Dist_intraThalamus =
            Vector2d<float>(thalPopulations, thalPopulations, -1.0);
        std::cerr << "Load intraThalamus_dist : " << file << std::endl;
        FILE* f = fopen(file, "r");
        assert(f);
        while (1) {
            int from, to;
            float dist;
            if (fscanf(f, "%d %d %f\n", &from, &to, &dist) != 3)
                break;
            from--;
            to--;
            assert(from < MRI_POINTS);
            assert(to < MRI_POINTS);
            Dist_intraThalamus(from, to) = dist;
        }
        fclose(f);
    }

    void load_weight_factors(const char* file) {
        std::cerr << "load_weight_factors : " << file << std::endl;
        timer run_timer;
        weight_factor = Vector4d<float>(
            SourceData::MRI_POINTS_FULL, SourceData::MRI_POINTS_FULL,
            SourceData::NUM_LAYERS, SourceData::NUM_LAYERS, -1.0);
        PRINT_RUNTIME_MEMUSED(run_timer, "Init weight_factor : ", std::cout);
        //

        std::error_code err_code;
        mio::mmap_source in_mmap;
        in_mmap.map(file, err_code);
        if (err_code) {
            handle_error(err_code);
            return;
        }

        run_timer.reset();
        const char* c_ptr = in_mmap.data();
        // int i = 0;
        do {
            int fromNeuron, toNeuron, fromLayer, toLayer;  //,factor;
            float factor;
            //
            char* end;
            fromNeuron = std::strtol(c_ptr, &end, 10);
            if (c_ptr == end) {
                break;
            }
            c_ptr = end;
            toNeuron = std::strtol(c_ptr, &end, 10);
            if (c_ptr == end) {
                break;
            }
            c_ptr = end;
            fromLayer = std::strtol(c_ptr, &end, 10);
            if (c_ptr == end) {
                break;
            }
            c_ptr = end;
            toLayer = std::strtol(c_ptr, &end, 10);
            if (c_ptr == end) {
                break;
            }
            c_ptr = end;
            factor = std::strtof(c_ptr, &end);
            if (c_ptr == end) {
                break;
            }
            c_ptr = end;

            // fmt::print("{:d} {:d} {:d} {:d} {:f} \n", fromNeuron, toNeuron,
            // fromLayer, toLayer, factor);
            fromNeuron--;
            toNeuron--;
            fromLayer--;
            toLayer--;  // MATLAB vs C indexing
            assert(fromNeuron < SourceData::MRI_POINTS_FULL);
            assert(toNeuron < SourceData::MRI_POINTS_FULL);
            weight_factor(fromNeuron, toNeuron, fromLayer, toLayer) = factor;

        } while (1);  // while
        PRINT_RUNTIME_MEMUSED(run_timer, "Loaded weight_factor : ", std::cout);
    }

    void load_weightFactors_INPY(const char* file) {
        weight_factor_INPY =
            Vector2d<float>(MRI_POINTS_FULL, MRI_POINTS_FULL, -1);
        std::cerr << "load_weightFactors_INPY" << std::endl;
        //  for (int i=0; i<MRI_POINTS_FULL; i++){
        //    for (int ii=0; ii<MRI_POINTS_FULL; ii++){
        //              weight_factor_INPY[i][ii] =-1;
        //    }
        //  }

        FILE* f = fopen(file, "r");
        assert(f);

        // int i = 0;
        while (1) {
            int fromNeuron, toNeuron;  //,factor;
            float factor;
            if (fscanf(f, "%d %d %f \n", &fromNeuron, &toNeuron, &factor) != 3)
                break;
            fromNeuron--;
            toNeuron--;  // MATLAB vs C indexing
            assert(fromNeuron < MRI_POINTS_FULL);
            assert(toNeuron < MRI_POINTS_FULL);
            weight_factor_INPY(fromNeuron, toNeuron) = factor;

        }  // while

        fclose(f);
    }

    //////////////// Read binary files /////////////////////

    void intReshape(std::vector<BYTE> fileData, int dim1, int dim2, int dim3,
                    int dim4) {
        for (int ll = 0; ll < dim4; ll++) {
            for (int kk = 0; kk < dim3; kk++) {
                for (int jj = 0; jj < dim2; jj++) {
                    for (int ii = 0; ii < dim1; ii++) {
                        weight_factor(ii, jj, kk, ll) =
                            fileData[ll * (dim3 * dim2 * dim1) +
                                     kk * (dim2 * dim1) + jj * dim1 + ii];
                    }
                }
            }
        }
    }

    void load_weightFactors_binary(const char* filename) {

        std::cerr << "load_weightFactors_binary : " << filename << std::endl;
        std::streampos fileSize;
        std::ifstream file(filename, std::ios::binary);

        // get its size:
        file.seekg(0, std::ios::end);
        fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        // read the data:
        std::vector<BYTE> fileData(fileSize);
        file.read(reinterpret_cast<char*>(&fileData[0]), fileSize);

        intReshape(fileData, 20484, 20484, 6, 6);
    }
    ////////////// End - read binary files /////////////////

  public:
    int load_source_data(const DataConfig& data_cfg) {
        timer load_timer;
        //
        useDelays = true;
        load_3D_subnet(data_cfg.left_subnet.c_str(), 0);   // left
        load_3D_subnet(data_cfg.right_subnet.c_str(), 1);  // right
        PRINT_RUNTIME_MEMUSED(load_timer, "Load 3D subnetwork : ", std::cout);

        //
        load_timer.reset();
        load_3D_LH_RH_correspondence(data_cfg.left_right_map.c_str(), true);
        load_3D_LH_RH_correspondence(data_cfg.left_right_map_small.c_str(),
                                     false);
        PRINT_RUNTIME_MEMUSED(load_timer,
                              "Load 3D LH-RH correspondence : ", std::cout);

        mapLayerNameId_IN = fillDict_IN(mapLayerNameId_IN);
        mapLayerNameId_PY = fillDict_PY(mapLayerNameId_PY);

        //
        load_timer.reset();
        load_3D_dist_prob(data_cfg.dist3d_probability.c_str());
        PRINT_RUNTIME_MEMUSED(load_timer,
                              "Load dist3d_probability : ", std::cout);
        //
        load_timer.reset();
        load_weightFactors_INPY(data_cfg.weight_factor_inpy.c_str());
        load_timer.measure_accumulate_print("Load IN PY : ", std::cout, false);
        print_memory_usage<uint64_t>("; ", std::cout);
        //
        //    load_weightFactors_binary(argv[3]);
        load_timer.reset();
        load_weight_factors(data_cfg.weight_factor.c_str());
        PRINT_RUNTIME_MEMUSED(load_timer, "Load Weight Factor : ", std::cout);

        load_timer.reset();
        load_thalCort_dist(data_cfg.thalamus_cortex_distance.c_str());
        PRINT_RUNTIME_MEMUSED(load_timer,
                              "Load Thalamic Cortex Dist : ", std::cout);

        load_timer.reset();
        load_intraThalamus_dist(data_cfg.intra_thalamus_distance.c_str());
        PRINT_RUNTIME_MEMUSED(load_timer,
                              "Load Intra Thalamic Dist : ", std::cout);
        return 0;
    }
};
#endif  // SOURCE_DATA_HPP
