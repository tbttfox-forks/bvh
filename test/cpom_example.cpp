
#include <bvh/v2/bvh.h>
#include <bvh/v2/vec.h>
#include <bvh/v2/ray.h>
#include <bvh/v2/node.h>
#include <bvh/v2/default_builder.h>
#include <bvh/v2/thread_pool.h>
#include <bvh/v2/executor.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/tri.h>

#include "load_obj.h"

#include <iostream>
#include <string>
#include <vector>
#include <tuple>

using Scalar  = float;
using Vec3    = bvh::v2::Vec<Scalar, 3>;
using BBox    = bvh::v2::BBox<Scalar, 3>;
using Tri     = bvh::v2::Tri<Scalar, 3>;
using Node    = bvh::v2::Node<Scalar, 3>;
using Bvh     = bvh::v2::Bvh<Node>;
using Ray     = bvh::v2::Ray<Scalar, 3>;


template <typename T>
auto read_queries(const std::string& path){
    std::ifstream is(file);
    static constexpr size_t max_line = 1024;
    char line[max_line];

    std::vector<Vec<T, 3>> queries, results;
    while (is.getline(ptr, max_line)) {
        if (*ptr == '\0' || *ptr == '#')
            continue;
        remove_eol(ptr);
        auto x = std::strtof(ptr, &ptr);
        auto y = std::strtof(ptr, &ptr);
        auto z = std::strtof(ptr, &ptr);
        queries.emplace_back(x, y, z);

        auto x = std::strtof(ptr, &ptr);
        auto y = std::strtof(ptr, &ptr);
        auto z = std::strtof(ptr, &ptr);
        results.emplace_back(x, y, z);
    }

    return std::make_tuple(queries, results);
}




int main() {
    // Obj file format
    std::string mandlebulbPath = "C:\\Users\\tyler\\src\\GitHub\\Libraries\\bvh\\test\\mandlebulb.obj";
    auto tris = load_obj<Scalar>(mandlebulbPath);

    // Text. 6 comma separated floats per line
    // Each line has 3 floats that are the query point, and 3 that are the result point
    // as determined by Maya's closest point algorithm
    std::string queryPath = "C:\\Users\\tyler\\src\\GitHub\\Libraries\\bvh\\test\\queries.txt";
    auto [queries, ground_truth_results] = read_queries<Scalar>(queryPath);

    bvh::v2::ThreadPool thread_pool;
    bvh::v2::ParallelExecutor executor(thread_pool);

    // Get triangle centers and bounding boxes (required for BVH builder)
    std::vector<BBox> bboxes(tris.size());
    std::vector<Vec3> centers(tris.size());
    executor.for_each(0, tris.size(), [&] (size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            bboxes[i]  = tris[i].get_bbox();
            centers[i] = tris[i].get_center();
        }
    });

    typename bvh::v2::DefaultBuilder<Node>::Config config;
    config.quality = bvh::v2::DefaultBuilder<Node>::Quality::High;
    auto bvh = bvh::v2::DefaultBuilder<Node>::build(thread_pool, bboxes, centers, config);

    // Permuting the primitive data allows to remove indirections during traversal, which makes it faster.
    static constexpr size_t invalid_id = std::numeric_limits<size_t>::max();

    auto best_prim_idx = invalid_id;
    Vec3 best_point(0), best_bary(0);


    auto leafFunc = [&](Scalar best_dist2, size_t begin, size_t end) {
        for (Index i = begin; i < end; ++i) {
            auto [prim_point, prim_bary] = closeest_point_tri(p, tris[i]);
            auto prim_dist2 = length_squared(prim_point - p);
            if (prim_dist2 < best_dist2) {

                best_prim_idx = i;
                best_point = prim_point;
                best_bary = prim_bary;

                best_dist2 = prim_dist2;
            }
        }
        return best_dist2;
    };

    for (size_t j = 0; j < queries.size(); ++j){
        auto& qp = queries[j];
        bvh.closest_point(qp, bvh.get_root().index, stack, leafFunc);
        //best_prim_idx
        //best_bary
        //best_point
        // allclose(best_point, ground_truth_results[i])
    }
    return 0;

}
