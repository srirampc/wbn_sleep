#include "fmt/core.h"
#include "fmt/format.h"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <random>
#include <vector>

#include <catch_amalgamated.hpp>
#include <source_data.hpp>
#include <vector_nd.hpp>

const char* DATA_CFG_PATH = "./config/data_config.json";

int handle_error(const std::error_code& error) {
    const auto& errmsg = error.message();
    std::cerr << "error mapping file:  " << errmsg << " exiting..."
              << std::endl;
    return error.value();
}

void validate_wfactor(const Vector4d<float>& weight_factor, const char* file) {
    timer run_timer;
    std::error_code err_code;
    mio::mmap_source in_mmap;
    in_mmap.map(file, err_code);
    if (err_code) {
        handle_error(err_code);
        return;
    }
    int nrecord = 0;

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
        if (weight_factor(fromNeuron, toNeuron, fromLayer, toLayer) != factor)
            std::cout << " -> Failed at record : " << nrecord << std::endl;
        REQUIRE(weight_factor(fromNeuron, toNeuron, fromLayer, toLayer) ==
                factor);

    } while (1);  // while
    PRINT_RUNTIME_MEMUSED(run_timer, "Validate weight_factor : ", std::cout);
}

void par_load_weight_factors(Vector4d<float>* pweight_factor,
                             const char* file) {
    timer run_timer;
    std::cerr << "par_load_weight_factors : " << file << std::endl;
    std::vector<std::size_t> vt_nrecords =
        par_mio_load_vector(pweight_factor, file);
    fmt::print("Par. Records {}\n", fmt::join(vt_nrecords, ", "));
    fmt::print("Par. Read weight_factor Records {}\n",
               std::accumulate(vt_nrecords.begin(), vt_nrecords.end(),
                               std::size_t(0)));
    PRINT_RUNTIME_MEMUSED(run_timer, "Par. Loaded weight_factor : ", std::cout);
}

void load_weight_factors(Vector4d<float>* pweight_factor, const char* file) {
    timer run_timer;
    std::cerr << "seq load_weight_factors : " << file << std::endl;
    mio_load_vector(pweight_factor, file);
    PRINT_RUNTIME_MEMUSED(run_timer, "Seq. Loaded weight_factor : ", std::cout);
}

TEST_CASE("Load Weight Factors", "[load_weight_factors]") {
    DataConfig data_cfg = DataConfig::parse(DATA_CFG_PATH);
    Vector4d<float> weight_factor(
        SourceData::MRI_POINTS_FULL, SourceData::MRI_POINTS_FULL,
        SourceData::NUM_LAYERS, SourceData::NUM_LAYERS, -1.0);
    load_weight_factors(&weight_factor, data_cfg.weight_factor.c_str());
    validate_wfactor(weight_factor, data_cfg.weight_factor.c_str());
}

TEST_CASE("Parallel Load Weight Factors", "[par_load_weight_factors]") {
    DataConfig data_cfg = DataConfig::parse(DATA_CFG_PATH);
    Vector4d<float> weight_factor(
        SourceData::MRI_POINTS_FULL, SourceData::MRI_POINTS_FULL,
        SourceData::NUM_LAYERS, SourceData::NUM_LAYERS, -1.0);
    par_load_weight_factors(&weight_factor, data_cfg.weight_factor.c_str());
    validate_wfactor(weight_factor, data_cfg.weight_factor.c_str());
}

void validate_intraThalamus_dist(const Vector2d<float>& Dist_intraThalamus,
                                 const char* file) {
    std::cerr << "Validate intraThalamus_dist : " << file << std::endl;
    std::size_t nrecord = 0;
    FILE* f = fopen(file, "r");
    assert(f);
    while (1) {
        int from, to;
        float dist;
        if (fscanf(f, "%d %d %f\n", &from, &to, &dist) != 3)
            break;
        from--;
        to--;
        assert(from < SourceData::MRI_POINTS);
        assert(to < SourceData::MRI_POINTS);
        nrecord++;
        if (Dist_intraThalamus(from, to) != dist) {
            std::cout << " -> Failed at Record : " << std::endl;
        }
        REQUIRE(Dist_intraThalamus(from, to) == dist);
    }
    fclose(f);
}

void load_intraThalamus_dist(Vector2d<float>* pDist_intraThalamus,
                             const char* file) {
    timer load_timer;
    std::cerr << "Load intraThalamus_dist : " << file << std::endl;
    mio_load_vector(pDist_intraThalamus, file);
    PRINT_RUNTIME_MEMUSED(load_timer,
                          "Load Thalamic Cortex Dist : ", std::cout);
}

TEST_CASE("Intra-Thalamus Distance", "[intra_thalamus_distance]") {
    DataConfig data_cfg = DataConfig::parse(DATA_CFG_PATH);
    Vector2d<float> Dist_intraThalamus(SourceData::thalPopulations,
                                       SourceData::thalPopulations, -1.0);
    load_intraThalamus_dist(&Dist_intraThalamus,
                            data_cfg.intra_thalamus_distance.c_str());
    //
    timer load_timer;
    validate_intraThalamus_dist(Dist_intraThalamus,
                                data_cfg.intra_thalamus_distance.c_str());
    PRINT_RUNTIME_MEMUSED(load_timer,
                          "Validate Thalamic Cortex Dist : ", std::cout);
}

float validate_3D_dist_prob(const std::array<Vector2d<float>, 2>& Dist_Prob,
                            const char* file) {
    std::cerr << "Validate 3D_dist_prob : " << file << std::endl;
    FILE* f = fopen(file, "r");
    assert(f);
    float maxDistance = -1.0;
    std::size_t nrecord = 0;
    // int i = 0;
    while (1) {
        int from, to;
        float dist, conn_prob;
        if (fscanf(f, "%d %d %f %f\n", &from, &to, &dist, &conn_prob) != 4)
            break;
        from--;
        to--;  // MATLAB vs C indexing
        assert(from < SourceData::MRI_POINTS);
        assert(to < SourceData::MRI_POINTS);
        if (Dist_Prob[0](from, to) != dist)
            std::cout << "  -> Failed dist at Records " << nrecord << std::endl;
        REQUIRE(Dist_Prob[0](from, to) == dist);
        //
        if (Dist_Prob[1](from, to) != conn_prob)
            std::cout << "  -> Failed conn prob at Records " << nrecord
                      << std::endl;
        REQUIRE(Dist_Prob[1](from, to) == conn_prob);
        if (dist > maxDistance)  // get max distance for scaling later
            maxDistance = dist;
    }  // while

    fclose(f);
    return maxDistance;
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

float load_3D_dist_prob(
    std::array<Vector2d<float>, SourceData::NUM_PARAM>& Dist_Prob,  // NOLINT
    const char* file) {
    timer load_timer;
    std::cerr << "Load 3D_dist_prob : " << file << std::endl;
    for (int ix = 0; ix < SourceData::NUM_PARAM; ix++) {
        Dist_Prob[ix] = Vector2d<float>(SourceData::MRI_POINTS,
                                        SourceData::MRI_POINTS, -1.0);
    }
    std::error_code err_code;
    mio::mmap_source in_mmap = mio::make_mmap_source(file, err_code);
    if (err_code) {
        handle_mio_error<int>(err_code);
    } else {
        const char* c_ptr = in_mmap.data();
        std::size_t nrecords = 0;
        do {
            if (parse_dist_prob_line1(&Dist_Prob[0], &Dist_Prob[1], c_ptr) != 0)
                break;
            nrecords++;
        } while (1);  // while
    }
    float maxDistance = *std::max_element(Dist_Prob[0].data_vec().begin(),
                                          Dist_Prob[0].data_vec().end());
    PRINT_RUNTIME_MEMUSED(load_timer, "Load 3D Dist : ", std::cout);
    return maxDistance;
}

TEST_CASE("Load 3D dist prob", "[load_3D_dist_prob]") {
    timer load_timer;
    DataConfig data_cfg = DataConfig::parse(DATA_CFG_PATH);
    std::array<Vector2d<float>, SourceData::NUM_PARAM> Dist_Prob;
    float maxDistance =
        load_3D_dist_prob(Dist_Prob, data_cfg.dist3d_probability.c_str());
    load_timer.reset();
    float testMax =
        validate_3D_dist_prob(Dist_Prob, data_cfg.dist3d_probability.c_str());
    PRINT_RUNTIME_MEMUSED(load_timer, "Validate 3D Dist : ", std::cout);
    if (testMax != maxDistance)
        std::cout << "  -> Failed conn prob at max dist " << maxDistance << " "
                  << testMax << std::endl;
    REQUIRE(testMax == maxDistance);
}

void validate_weight_factors_INPY(const Vector2d<float>& weight_factor_INPY,
                                  const char* file) {
    std::cerr << "load_weightFactors_INPY" << std::endl;
    std::size_t nrecord = 0;
    FILE* f = fopen(file, "r");
    assert(f);
    while (1) {
        int fromNeuron, toNeuron;  //,factor;
        float factor;
        if (fscanf(f, "%d %d %f \n", &fromNeuron, &toNeuron, &factor) != 3)
            break;
        fromNeuron--;
        toNeuron--;  // MATLAB vs C indexing
        assert(fromNeuron < SourceData::MRI_POINTS_FULL);
        assert(toNeuron < SourceData::MRI_POINTS_FULL);
        nrecord++;
        if (weight_factor_INPY(fromNeuron, toNeuron) != factor)
            std::cout << "Failed at Records " << nrecord << std::endl;
        REQUIRE(weight_factor_INPY(fromNeuron, toNeuron) == factor);
    }  // while
    fclose(f);
}

void par_load_weight_factors_INPY(Vector2d<float>* pweight_factor_INPY,
                                  const char* file) {
    timer load_timer;
    std::cerr << "par_load_weightFactors_INPY" << std::endl;
    std::vector<std::size_t> vt_nrecords =
        par_mio_load_vector(pweight_factor_INPY, file);
    //
    fmt::print("Parallel Read weight_factor_INPY Records {} \n",
               std::accumulate(vt_nrecords.begin(), vt_nrecords.end(),
                               std::size_t(0)));
}

void load_weight_factors_INPY(Vector2d<float>* pweight_factor_INPY,
                              const char* file) {
    timer load_timer;
    std::cerr << "load_weightFactors_INPY" << std::endl;
    mio_load_vector(pweight_factor_INPY, file);
    PRINT_RUNTIME_MEMUSED(load_timer,
                          "Load Weight Factors IN PY : ", std::cout);
    load_timer.reset();
}

TEST_CASE("Load Weight Factors INPY", "[load_weight_factors_inpy]") {
    DataConfig data_cfg = DataConfig::parse(DATA_CFG_PATH);
    Vector2d<float> weight_factor_INPY(SourceData::MRI_POINTS_FULL,
                                       SourceData::MRI_POINTS_FULL, -1);
    load_weight_factors_INPY(&weight_factor_INPY,
                             data_cfg.weight_factor_inpy.c_str());
    validate_weight_factors_INPY(weight_factor_INPY,
                                 data_cfg.weight_factor_inpy.c_str());
}

TEST_CASE("Parallel Load Weight Factors INPY",
          "[par_load_weight_factors_inpy]") {
    DataConfig data_cfg = DataConfig::parse(DATA_CFG_PATH);
    Vector2d<float> weight_factor_INPY(SourceData::MRI_POINTS_FULL,
                                       SourceData::MRI_POINTS_FULL, -1);
    par_load_weight_factors_INPY(&weight_factor_INPY,
                                 data_cfg.weight_factor_inpy.c_str());
    validate_weight_factors_INPY(weight_factor_INPY,
                                 data_cfg.weight_factor_inpy.c_str());
}
//
// int random_test() {
//     std::vector<std::mt19937*> vrGenPtr;
//     std::vector<std::mt19937::result_type> results;
// #pragma omp parallel default(none) shared(vrGenPtr, results)
//     {
//         int n_threads = omp_get_num_threads();
//         int thread_id = omp_get_thread_num();
// #pragma omp single
//         {
//             vrGenPtr.resize(n_threads);
//             results.resize(n_threads);
//         }
//         std::random_device rd;
//         results[thread_id] = rd();
//         vrGenPtr[thread_id] = new std::mt19937(results[thread_id]);
//     }
//     fmt::print("Rand Values {} \n", fmt::join(results, ", "));
//     for (std::mt19937* rptr : vrGenPtr) {
//         delete rptr;
//     }
//     return 0;
// }
//
// int main(int argc, char* argv[]) {
//     AppArguments in_args;
//     CLI11_PARSE(in_args.app, argc, argv);  // Parse command line arguments
//     in_args.print();
//     if (in_args.data_cfg_path.length() > 0) {
//         try {
//             in_args.data_cfg = DataConfig::parse(in_args.data_cfg_path);
//         } catch (njson::exception& e) {
//             std::cerr << "Failed to parse Data Config File "
//                       << in_args.data_cfg_path << ":" << e.what() <<
//                       std::endl;
//             return -1;
//         }
//         in_args.data_cfg.print();
//     }
//
//     // load_3D_dist_prob(in_args.data_cfg.dist3d_probability.c_str());
//     //
//     load_intraThalamus_dist(in_args.data_cfg.intra_thalamus_distance.c_str());
//     // load_weight_factors_INPY(in_args.data_cfg.weight_factor_inpy.c_str());
//     //
//     par_load_weight_factors_INPY(in_args.data_cfg.weight_factor_inpy.c_str());
//     // load_weight_factors(in_args.data_cfg.weight_factor.c_str());
//     par_load_weight_factors(in_args.data_cfg.weight_factor.c_str());
//     return 0;
// }
