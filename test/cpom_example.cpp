
#include <bvh/v2/bvh.h>
#include <bvh/v2/vec.h>
#include <bvh/v2/ray.h>
#include <bvh/v2/node.h>
#include <bvh/v2/default_builder.h>
#include <bvh/v2/thread_pool.h>
#include <bvh/v2/executor.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/tri.h>
#include <bvh/v2/dist_point_triangle.h>

#include "load_obj.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

using Scalar  = double;
using Index   = size_t;
using Vec3    = bvh::v2::Vec<Scalar, 3>;
using BBox    = bvh::v2::BBox<Scalar, 3>;
using Tri     = bvh::v2::Tri<Scalar, 3>;
using Node    = bvh::v2::Node<Scalar, 3>;
using Bvh     = bvh::v2::Bvh<Node>;
using Ray     = bvh::v2::Ray<Scalar, 3>;


auto read_queries(const std::string& path){

    std::ifstream infile(path);
    std::string line;
    std::vector<Vec3> queries, results;
    while (std::getline(infile, line)) {
        Scalar x, y, z, a, b, c;
        std::istringstream iss(line);
        if (!(iss >> x >> y >> z >> a >> b >> c)) { break; }
        
        queries.emplace_back(x, y, z);
        results.emplace_back(a, b, c);

    }

    return std::make_tuple(queries, results);
}

bool isclose(Scalar a, Scalar b) {
    return std::fabs(a - b) < 1.0e-5;
}


int main() {
    // Obj file format
    std::string mandlebulbPath = "C:\\Users\\tyler\\src\\GitHub\\Libraries\\bvh\\test\\mandleBulb.obj";
    auto tris = load_obj<Scalar>(mandlebulbPath);

    // Text. 6 comma separated floats per line
    // Each line has 3 floats that are the query point, and 3 that are the result point
    // as determined by Maya's closest point algorithm
    std::string queryPath = "C:\\Users\\tyler\\src\\GitHub\\Libraries\\bvh\\test\\mandleBulbQueries.txt";
    auto [queries, ground_truth_results] = read_queries(queryPath);

    bvh::v2::ThreadPool thread_pool;
    bvh::v2::ParallelExecutor executor(thread_pool);

    // Get triangle centers and bounding boxes (required for BVH builder)
    std::vector<BBox> bboxes(tris.size());
    std::vector<Vec3> centers(tris.size());
    std::vector<Vec3> normals(tris.size());
    executor.for_each(0, tris.size(), [&] (size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            bboxes[i]  = tris[i].get_bbox();
            centers[i] = tris[i].get_center();
            normals[i] = tris[i].get_normal();
        }
    });

    typename bvh::v2::DefaultBuilder<Node>::Config config;
    config.quality = bvh::v2::DefaultBuilder<Node>::Quality::High;
    auto bvh = bvh::v2::DefaultBuilder<Node>::build(thread_pool, bboxes, centers, config);

    static constexpr size_t invalid_id = std::numeric_limits<size_t>::max();

    auto best_prim_idx = invalid_id;
    Vec3 best_point(0), best_bary(0);
    Scalar best_dist2 = std::numeric_limits<Scalar>::max();

    auto leafFunc = [&](Vec3& p, size_t begin, size_t end) {
        for (Index i = begin; i < end; ++i) {
            auto [prim_point, prim_bary] = bvh::v2::closest_point_tri(p, tris[bvh.prim_ids[i]]);
            auto prim_dist2 = bvh::v2::length_squared<Scalar, 3>(prim_point - p);
            if (prim_dist2 < best_dist2) {

                best_prim_idx = i;
                best_point = prim_point;
                best_bary = prim_bary;

                best_dist2 = prim_dist2;
            }
        }
        return best_dist2;
    };

    static constexpr size_t stack_size = 64;
    size_t badCount = 0;

    for (size_t zz = 0; zz < 10000; ++zz){
        for (size_t j = 0; j < queries.size(); ++j){
            auto& qp = queries[j];
            bvh::v2::SmallStack<Bvh::Index, stack_size> nodeStack;

            best_dist2 = std::numeric_limits<Scalar>::max();
            bvh.closest_point(qp, bvh.get_root(), nodeStack, leafFunc, best_dist2);
            //std::cout << "Closest TriIdx " << best_prim_idx << std::endl;
            //std::cout << "Closest Bary" << best_bary[0] << ", " << best_bary[1] << ", " << best_bary[2] << std::endl;
            //std::cout << "Closest Point" << best_point[0] << ", " << best_point[1] << ", " << best_point[2] << std::endl;

            // allclose(best_point, ground_truth_results[i])
            auto &gt = ground_truth_results[j];
            for (size_t c = 0; c < 3; ++c) {
                if (!isclose(gt[c], best_point[c])) {
                    //std::cout << "BAD!!!" << std::endl;
                    //std::cout << "Closest Point: " << best_point[0] << ", " << best_point[1] << ", " << best_point[2] << std::endl;
                    //std::cout << "Ground Truth : " << gt[0] << ", " << gt[1] << ", " << gt[2] << std::endl;
                    badCount++;
                    break;
                }
            }
        }
    }

    std::cout << "BadCount: " << badCount << std::endl;

    return 0;
}
