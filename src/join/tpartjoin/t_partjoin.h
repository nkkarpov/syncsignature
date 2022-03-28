//
// Created by nkarpov on 3/24/21.
//
#pragma once

#include "binary_node.h"
#include "string_label.h"
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace join {
using std::map;
using std::unordered_map;
using std::unordered_set;
using std::vector;

/*template <typename Label, typename VerificationAlgorithm>
class TPartJoinTI {
public:
    void
    execute_join(std::vector<node::Node<Label>>& trees_collection,
                 std::vector<node::BinaryNode<Label>>& binary_trees_collection,
                 std::vector<std::pair<int, int>>& candidates,
                 std::vector<join::JoinResultElement>& join_result, int T,
                 int t2bin, double sim, double distance_threshold);

    void verify_candidates(std::vector<node::Node<Label>>& trees_collection,
                           std::vector<std::pair<int, int>>& candidates,
                           std::vector<join::JoinResultElement>& join_result,
                           double distance_threshold);

    void convert_trees_to_binary_trees(
        std::vector<node::Node<Label>>& trees_collection,
        std::vector<node::BinaryNode<Label>>& binary_trees_collection);

    void
    retrieve_candidates(std::vector<std::vector<size_t>>& tours_collection,
                        std::vector<std::vector<size_t>>& signatures_collection,
                        std::vector<std::pair<int, int>>& candidates,
                        double sim, double distance_threshold);

    void lowerbound(const std::vector<std::vector<size_t>>& tours_collection,
                    std::vector<std::pair<int, int>>& candidates,
                    double distance_threshold);
    void
    convert_trees_to_tours(std::vector<node::Node<Label>>& trees_collection,
                           std::vector<std::vector<size_t>>& tours_collection);
    void partition_binary_tree(
        std::vector<node::BinaryNode<Label>>& binary_trees_collection,
        std::vector<std::vector<size_t>>& signature_collection, int T);

    void partition_tree(std::vector<node::Node<Label>>& trees_collection,
                        std::vector<std::vector<size_t>>& signature_collection,
                        int T);
    void upperbound(std::vector<node::Node<Label>>& trees_collection,
                    std::vector<std::pair<int, int>>& candidates,
                    std::vector<join::JoinResultElement>& join_result,
                    double distance_threshold);
};
#include "t_partjoin_impl.h"*/
} // namespace join