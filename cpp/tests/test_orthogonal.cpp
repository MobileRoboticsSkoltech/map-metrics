#include <gtest/gtest.h>

#include "parse_utils.h"

#include "max_clique.h"
#include "metrics.h"
#include "orth_extract.h"
#include "clustering.h"

#include <iostream>
#include <chrono>

#ifdef TEST_ORTHOGONAL_METRICS

TEST(Orthogonal, MOM){
    std::vector<cilantro::VectorSet3d> PCs = parse_utils::depth_parse::GetDepthPCs("data/depth/pcs");
    std::vector<Eigen::Matrix4d> poses = parse_utils::depth_parse::GetDepthPoses("data/depth/poses");
    auto pc = cilantro::PointCloud3d(PCs[2]);  // PC-0091.ply

    auto orth_subset = orth_extract::ExtractOrthogonalSubset(pc, config::CustomConfig(5, 0.2, 30, 5));

    auto start_time = std::chrono::system_clock::now();
    double result_mom = metrics::ComputeOrthogonalMetrics(
                            PCs,
                            poses,
                            config::CustomConfig(5, 0.2, 30, 5),
                            &metrics_utils::metrics_algorithm::ComputeEigenvalues,
                            orth_subset
                        );  
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    
    std::cout << "MOM: " << result_mom << std::endl;
    std::cout << "Elapsed time (ms): " << elapsed_milliseconds.count() << std::endl;

    ASSERT_LE(std::fabs(result_mom - 0.00203637), 1e-8);
}

#endif // TEST_ORTHOGONAL_METRICS