#ifndef BVH_V2_BVH_H
#define BVH_V2_BVH_H

#include "bvh/v2/node.h"

#include <cstddef>
#include <iterator>
#include <vector>
#include <stack>
#include <utility>

namespace bvh::v2 {

template <typename Node>
struct Bvh {
    using Index = typename Node::Index;
    using Scalar = typename Node::Scalar;

    std::vector<Node> nodes;
    std::vector<size_t> prim_ids;

    Bvh() = default;
    Bvh(Bvh&&) = default;

    Bvh& operator = (Bvh&&) = default;

    bool operator == (const Bvh& other) const = default;
    bool operator != (const Bvh& other) const = default;

    /// Returns the root node of this BVH.
    const Node& get_root() const { return nodes[0]; }

    /// Extracts the BVH rooted at the given node index.
    inline Bvh extract_bvh(size_t root_id) const;

    Octant get_vec_octant(Vec<Scalar, Node::dimension>& p) const;

    template <typename Stack>
    void closest_point(Vec<Scalar, Node::dimension>& p, Index, Stack&) const;

    /// Intersects the BVH with a single ray, using the given function to intersect the contents
    /// of a leaf. The algorithm starts at the node index `top` and uses the given stack object.
    /// When `IsAnyHit` is true, the function stops at the first intersection (useful for shadow
    /// rays), otherwise it finds the closest intersection. When `IsRobust` is true, a slower but
    /// numerically robust ray-box test is used, otherwise a fast, but less precise test is used.
    template <bool IsAnyHit, bool IsRobust, typename Stack, typename LeafFn, typename InnerFn = IgnoreArgs>
    inline void intersect(Ray<Scalar, Node::dimension>& ray, Index top, Stack&, LeafFn&&, InnerFn&& = {}) const;

    inline void serialize(OutputStream&) const;
    static inline Bvh deserialize(InputStream&);
};

template <typename Node>
auto Bvh<Node>::extract_bvh(size_t root_id) const -> Bvh {
    assert(root_id != 0);

    Bvh bvh;
    bvh.nodes.emplace_back();

    std::stack<std::pair<size_t, size_t>> stack;
    stack.emplace(root_id, 0);
    while (!stack.empty()) {
        auto [src_id, dst_id] = stack.top();
        stack.pop();
        auto& src_node = nodes[src_id];
        auto& dst_node = bvh.nodes[dst_id];
        dst_node = src_node;
        if (src_node.is_leaf()) {
            dst_node.index.first_id = static_cast<typename Index::Type>(bvh.prim_ids.size());
            std::copy_n(
                prim_ids.begin() + src_node.index.first_id,
                src_node.index.prim_count,
                std::back_inserter(bvh.prim_ids));
        } else {
            size_t first_id = bvh.nodes.size();
            dst_node.index.first_id = static_cast<typename Index::Type>(first_id);
            bvh.nodes.emplace_back();
            bvh.nodes.emplace_back();
            stack.emplace(src_node.index.first_id + 0, first_id + 0);
            stack.emplace(src_node.index.first_id + 1, first_id + 1);
        }
    }
    return bvh;
}

template <typename Node>
BVH_ALWAYS_INLINE Octant Bvh<Node>::get_vec_octant(Vec<Scalar, Node::dimension> &p) const {
    static_assert(Node::dimension <= Octant::max_dim);
    Octant octant;
    static_for<0, Node::dimension>([&] (size_t i) {
        octant.value |= std::signbit(p[i]) * (uint32_t{1} << i);
    });
    return octant;
}

template <typename Node>
template <typename Stack>
void Bvh<Node>::closest_point(Vec<Scalar, Node::dimension>& p, Index start, Stack& stack) const{

    // Initialize the return values
    Scalar best_dist2 = std::numeric_limits<Scalar>::max();
    Vec<Scalar, Node::dimension> best_point, best_bary;
    Index best_prim_idx = std::numeric_limits<Index>::max();

    // Add root to the stack
    stack.push(start);

    // keep going until the stack is empty
restart:
    while (!stack.is_empty()) {
        // Take closest node out of the stack
        auto top = stack.pop();

        // Start drilling down the hierarchy, picking the closer of the two children
        // each time as the "top", and putting the other on the stack
        while (top.prim_count == 0) {
            auto& left  = nodes[top.first_id];
            auto& right = nodes[top.first_id + 1];

            // TODO: I don't know the time/memory tradeoff for storing these values
            // or computing them. Perhaps that's something I could template
            Scalar near_dist2 = length_squared(left.get_bbox().closest_point(p));
            Scalar far_dist2 = length_squared(right.get_bbox().closest_point(p));

            // Pick the closer one to keep checking
            // and store the other one on the stack
            auto near = left.index;
            auto far = right.index;
            if (far_dist2 < near_dist2) {
                std::swap(near, far);
                std::swap(near_dist2, far_dist2);
            }

            // If the near node is farther than my best dist
            // then I can just skip everything, and restart the while loop
            if (near_dist2 > best_dist2) goto restart;
            top = near;
            if (far_dist2 > best_dist2) continue;
            stack.push(far)
        }

        // Here we are at a leaf node
        // To keep the same structure, perhaps I should do this as a closure
        for (Index i = top.first_id, i < top.first_id + top.prim_count; ++i) {
            Vec<Scalar, Node::dimension> prim_point, prim_bary;
            auto prim_dist2 = point_tri_dist(precomputed_tris[i], p, prim_point, prim_bary);
            if (prim_dist2 < best_dist2) {
                best_prim_idx = i;
                best_point = prim_point;
                best_bary = prim_bary;
            }
        }
    }




}



template <typename Node>
template <bool IsAnyHit, bool IsRobust, typename Stack, typename LeafFn, typename InnerFn>
void Bvh<Node>::intersect(Ray<Scalar, Node::dimension>& ray, Index start, Stack& stack, LeafFn&& leaf_fn, InnerFn&& inner_fn) const {
    auto inv_dir = ray.get_inv_dir();
    auto inv_org = -inv_dir * ray.org;
    auto inv_dir_pad = Ray<Scalar, Node::dimension>::pad_inv_dir(inv_dir);
    auto octant = ray.get_octant();

    auto intersect_node = [&] (const Node& node) {
        return IsRobust
            ? node.intersect_robust(ray, inv_dir, inv_dir_pad, octant)
            : node.intersect_fast(ray, inv_dir, inv_org, octant);
    };

    stack.push(start);
restart:
    while (!stack.is_empty()) {
        auto top = stack.pop();
        while (top.prim_count == 0) {
            auto& left  = nodes[top.first_id];
            auto& right = nodes[top.first_id + 1];

            inner_fn(left, right);

            auto intr_left  = intersect_node(left);
            auto intr_right = intersect_node(right);

            bool hit_left  = intr_left.first <= intr_left.second;
            bool hit_right = intr_right.first <= intr_right.second;

            if (hit_left) {
                auto near = left.index;
                if (hit_right) {
                    auto far = right.index;
                    if (!IsAnyHit && intr_left.first > intr_right.first)
                        std::swap(near, far);
                    stack.push(far);
                }
                top = near;
            } else if (hit_right)
                top = right.index;
            else [[unlikely]]
                goto restart;
        }

        [[maybe_unused]] auto was_hit = leaf_fn(top.first_id, top.first_id + top.prim_count);
        if constexpr (IsAnyHit) {
            if (was_hit) return;
        }
    }
}

template <typename Node>
void Bvh<Node>::serialize(OutputStream& stream) const {
    stream.write(nodes.size());
    stream.write(prim_ids.size());
    for (auto&& node : nodes)
        node.serialize(stream);
    for (auto&& prim_id : prim_ids)
        stream.write(prim_id);
}

template <typename Node>
Bvh<Node> Bvh<Node>::deserialize(InputStream& stream) {
    Bvh bvh;
    bvh.nodes.resize(stream.read<size_t>());
    bvh.prim_ids.resize(stream.read<size_t>());
    for (auto& node : bvh.nodes)
        node = Node::deserialize(stream);
    for (auto& prim_id : bvh.prim_ids)
        prim_id = stream.read<size_t>();
    return bvh;
}

} // namespace bvh::v2

#endif
