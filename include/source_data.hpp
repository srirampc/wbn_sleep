#ifndef SOURCE_DATA_HPP
#define SOURCE_DATA_HPP

#include <cassert>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <omp.h>
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

template <typename RT> RT handle_mio_error(const std::error_code& error) {
    const auto& errmsg = error.message();
    std::cerr << "error mapping file:  " << errmsg << " exiting..."
              << std::endl;
    return RT(error.value());
}

template <typename VT>
std::vector<std::size_t> par_mio_load_vector(VT* pout_vec, const char* file) {
    std::uintmax_t f_size = std::filesystem::file_size(file);
    std::vector<std::size_t> vt_nrecords;
    //
#pragma omp parallel default(none) shared(pout_vec, file, f_size, vt_nrecords)
    {
        int n_threads = omp_get_num_threads();
        int tid = omp_get_thread_num();
#pragma omp single
        {
            vt_nrecords.resize(n_threads, 0);
        }
        std::uintmax_t m_offst = block_low<uintmax_t>(tid, n_threads, f_size);
        std::uintmax_t m_size = block_size<uintmax_t>(tid, n_threads, f_size);
        //
        std::error_code err_code;
        mio::mmap_source in_mmap =
            mio::make_mmap_source(file, m_offst, m_size, err_code);
        if (err_code) {
            handle_mio_error<int>(err_code);
        } else {
            const char* c_ptr = in_mmap.data();
            if (tid != 0)
                while (*c_ptr != '\n')
                    c_ptr++;
            //
            vt_nrecords[tid] =
                load_vector_until(pout_vec, c_ptr, c_ptr + m_size);
        }
    }
    return vt_nrecords;
}

template <typename VT>
std::size_t mio_load_vector(VT* pout_vec, const char* file) {  // NOLINT
    std::error_code err_code;
    mio::mmap_source in_mmap = mio::make_mmap_source(file, err_code);
    if (err_code) {
        handle_mio_error<int>(err_code);
        return 0;
    } else {
        const char* c_ptr = in_mmap.data();
        return load_vector(pout_vec, c_ptr);
    }
}

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
    float maxDistance;  // max distance between two neurons
    float max_TC_Distance;
    bool useDelays;
    //
    // mapping between LH-RH neurons
    Vector1d<int> map_LH_RH;
    //
    // mapping between a subset of LH-RH neurons
    Vector1d<int> map_LH_RH_small;
    //
    //
    Vector1d<float> TC_Dist;
    //
    // distance lookup table given for MRI data; [0][][] for left hemispher
    std::array<Vector2d<float>, 2> Distance3D;
    //
    // The actual number of smaller net is to be read out from network.cfg (TC
    // cells)
    // distance lookup table given for MRI data, [0][] for left hemisphere
    std::array<Vector1d<float>, 2> Subnet;
    //
    // lookup table given for MRI data; [0][][] for left hemispher
    std::array<Vector2d<float>, NUM_PARAM> Dist_Prob;
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
    explicit SourceData(const NetConfig& in_nc)
        : nc(in_nc), maxDistance(0), max_TC_Distance(0), useDelays(false),
          map_LH_RH(MAX_SUBNET_POINTS, -1),
          map_LH_RH_small(MAX_SUBNET_POINTS, -1), TC_Dist(NUM_TC_POINTS, -1.0) {
    }
    ~SourceData() {}

    // Access functions
    inline const std::array<Vector1d<float>, 2>& SubnetRef() const {
        return Subnet;
    }

    inline const Vector1d<float>& Subnet_Ref(bool right) const {
        return Subnet[right];
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
        return Dist_Prob[0](id_from, id_to);
    }

    inline double prob_only(int id_from, int id_to) const {
        // assert(Dist_Prob[id_from][id_to][1]!=-1);
        return Dist_Prob[1](id_from, id_to);
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
            Subnet[right_hemisphere] = Vector1d<float>(MAX_SUBNET_POINTS, -1);
        }

        // for (int i = 0; i < MAX_SUBNET_POINTS; i++)
        //    Subnet[right_hemisphere](i) = -1;

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

    void load_3D_LH_RH_correspondence(Vector1d<int>* pMap, const char* file) {
        std::cerr << "Loading 3D LH-RH interhemispheric mapping : " << file
                  << std::endl;

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

            pMap->set(id, map);
            // printf("Mapping id: %d to %d \n",id,map);
        }  // while

        fclose(f);
    }

    inline void load_3D_LH_RH_correspondence(const char* data_file,
                                             bool full_data) {
        if (full_data) {
            load_3D_LH_RH_correspondence(&map_LH_RH, data_file);
        } else {
            load_3D_LH_RH_correspondence(&map_LH_RH_small, data_file);
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

    inline int parse_dist_prob_line1(Vector2d<float>* pdist_vec,
                                     Vector2d<float>* pprob_vec,
                                     const char*& c_ptr) {  // NOLINT
        char* end;
        int i = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        int j = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        float v1 = std::strtof(c_ptr, &end);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        float v2 = std::strtof(c_ptr, &end);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        pdist_vec->set(--i, --j, v1);
        pprob_vec->set(i, j, v2);
        return 0;
    }

    void load_3D_dist_prob(const char* file) {
        timer load_timer;
        std::cerr << "Load 3D_dist_prob : " << file << std::endl;
        for (int ix = 0; ix < NUM_PARAM; ix++) {
            Dist_Prob[ix] = Vector2d<float>(MRI_POINTS, MRI_POINTS, -1.0);
        }
        std::error_code err_code;
        mio::mmap_source in_mmap = mio::make_mmap_source(file, err_code);
        if (err_code) {
            handle_mio_error<int>(err_code);
        } else {
            const char* c_ptr = in_mmap.data();
            std::size_t nrecords = 0;
            do {
                if (parse_dist_prob_line1(&Dist_Prob[0], &Dist_Prob[1],
                                          c_ptr) != 0)
                    break;
                nrecords++;
            } while (1);  // while
        }
        maxDistance = *std::max_element(Dist_Prob[0].data_vec().begin(),
                                        Dist_Prob[0].data_vec().end());
        PRINT_RUNTIME_MEMUSED(load_timer, "Load 3D Dist : ", std::cout);
    }

    void load_intraThalamus_dist(const char* file) {
        Dist_intraThalamus =
            Vector2d<float>(thalPopulations, thalPopulations, -1.0);
        std::cerr << "Load intraThalamus_dist : " << file << std::endl;
        mio_load_vector(&Dist_intraThalamus, file);
    }

    void par_load_weight_factors(const char* file) {
        std::cerr << "load_weight_factors : " << file << std::endl;
        timer run_timer;
        weight_factor = Vector4d<float>(MRI_POINTS_FULL, MRI_POINTS_FULL,
                                        NUM_LAYERS, NUM_LAYERS, -1.0);
        PRINT_RUNTIME_MEMUSED(run_timer, "Init weight_factor : ", std::cout);
        //
        run_timer.reset();
        std::vector<std::size_t> vt_nrecords =
            par_mio_load_vector(&weight_factor, file);
        fmt::print("Read weight_factor Records {} \n",
                   std::accumulate(vt_nrecords.begin(), vt_nrecords.end(), 0));
        PRINT_RUNTIME_MEMUSED(run_timer, "Loaded weight_factor : ", std::cout);
    }

    void load_weight_factors(const char* file) {
        std::cerr << "seq load_weight_factors : " << file << std::endl;
        timer run_timer;
        weight_factor = Vector4d<float>(MRI_POINTS_FULL, MRI_POINTS_FULL,
                                        NUM_LAYERS, NUM_LAYERS, -1.0);
        PRINT_RUNTIME_MEMUSED(run_timer,
                              "Seq Init weight_factor : ", std::cout);
        //
        run_timer.reset();
        mio_load_vector(&weight_factor, file);
        PRINT_RUNTIME_MEMUSED(run_timer,
                              "Seq Loaded weight_factor : ", std::cout);
    }

    void par_load_weight_factors_INPY(const char* file) {
        weight_factor_INPY =
            Vector2d<float>(MRI_POINTS_FULL, MRI_POINTS_FULL, -1);
        std::cerr << "load_weightFactors_INPY" << std::endl;
        std::vector<std::size_t> vt_nrecords =
            par_mio_load_vector(&weight_factor_INPY, file);
        fmt::print("Read weight_factor_INPY Records {} \n",
                   std::accumulate(vt_nrecords.begin(), vt_nrecords.end(), 0));
    }

    void load_weight_factors_INPY(const char* file) {
        weight_factor_INPY =
            Vector2d<float>(MRI_POINTS_FULL, MRI_POINTS_FULL, -1);
        std::cerr << "load_weightFactors_INPY" << std::endl;
        mio_load_vector(&weight_factor_INPY, file);
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
    int load_source_data(const DataConfig& data_cfg, bool parallel_run) {
        timer load_timer;
        //
        useDelays = true;
        load_3D_subnet(data_cfg.left_subnet.c_str(), 0);   // left
        load_3D_subnet(data_cfg.right_subnet.c_str(), 1);  // right
        //
        load_3D_LH_RH_correspondence(data_cfg.left_right_map.c_str(), true);
        load_3D_LH_RH_correspondence(data_cfg.left_right_map_small.c_str(),
                                     false);
        PRINT_RUNTIME_MEMUSED(load_timer,
                              "Load 3D Subnet & LH-RH corres. : ", std::cout);

        mapLayerNameId_IN = fillDict_IN(mapLayerNameId_IN);
        mapLayerNameId_PY = fillDict_PY(mapLayerNameId_PY);

        //
        load_timer.reset();
        load_3D_dist_prob(data_cfg.dist3d_probability.c_str());
        PRINT_RUNTIME_MEMUSED(load_timer,
                              "Load dist3d_probability : ", std::cout);
        //
        load_timer.reset();
        if (parallel_run) {
            par_load_weight_factors_INPY(data_cfg.weight_factor_inpy.c_str());
        } else {
            load_weight_factors_INPY(data_cfg.weight_factor_inpy.c_str());
        }
        load_timer.measure_accumulate_print("Load IN PY : ", std::cout, false);
        print_memory_usage<uint64_t>("; ", std::cout);
        //
        //    load_weightFactors_binary(argv[3]);
        load_timer.reset();
        if (parallel_run) {
            par_load_weight_factors(data_cfg.weight_factor.c_str());
        } else {
            load_weight_factors(data_cfg.weight_factor.c_str());
        }
        PRINT_RUNTIME_MEMUSED(load_timer, "Load Weight Factor : ", std::cout);

        load_timer.reset();
        load_thalCort_dist(data_cfg.thalamus_cortex_distance.c_str());
        load_intraThalamus_dist(data_cfg.intra_thalamus_distance.c_str());
        PRINT_RUNTIME_MEMUSED(
            load_timer,
            "Load Thalamus-Cortex & Intra Thalamic Dist : ", std::cout);
        return 0;
    }
};
#endif  // SOURCE_DATA_HPP
