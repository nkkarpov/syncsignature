//
// Created by nkarpov on 3/24/21.
//
#pragma once
#include "binary_node.h"
#include "string_label.h"
#include <join_result_element.h>
#include <map>
#include <t_join_ti.h>
#include <t_minjoin.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>
#include <random>

namespace join {
using std::map;
using std::unordered_map;
using std::unordered_set;
using std::vector;

template <typename Label, typename VerificationAlgorithm>
class TBallJoinTI {
public:
    void execute_parallel_join(
        std::vector<node::Node<Label>>& trees_collection,
        std::vector<node::BinaryNode<Label>>& binary_trees_collection,
        std::vector<std::pair<int, int>>& candidates,
        std::vector<join::JoinResultElement>& join_result, double fraction,
        double tau, double distance_threshold, int p, int number_of_threads);

    /*void
    execute_join(std::vector<node::Node<Label>>& trees_collection,
                 std::vector<node::BinaryNode<Label>>& binary_trees_collection,
                 std::vector<std::pair<int, int>>& candidates,
                 std::vector<join::JoinResultElement>& join_result,
                 double fraction, double tau, double distance_threshold, int p);*/

    void verify_candidates(std::vector<node::Node<Label>>& trees_collection,
                           std::vector<std::pair<int, int>>& candidates,
                           std::vector<join::JoinResultElement>& join_result,
                           double distance_threshold);

    void convert_trees_to_binary_trees(
        std::vector<node::Node<Label>>& trees_collection,
        std::vector<node::BinaryNode<Label>>& binary_trees_collection);

    void lowerbound(const std::vector<std::vector<size_t>>& preorder,
                    std::vector<std::pair<int, int>>& candidates,
                    double distance_threshold);
    void convert_trees_to_tours(
        std::vector<node::Node<Label>>& trees_collection,
        std::vector<std::vector<std::pair<size_t, int>>>& tours_collection,
        std::vector<std::vector<size_t>>& preorder, int p);
    void partition_binary_tree(
        std::vector<node::BinaryNode<Label>>& binary_trees_collection,
        std::vector<std::vector<std::pair<size_t, int>>>& tours_collection,
        std::vector<signature_t>& signature_collection, int first, int last,
        int z, double tau, int p, int shift = 0, int step = 1);

    void upperbound(std::vector<node::Node<Label>>& trees_collection,
                    std::vector<std::pair<int, int>>& candidates,
                    std::vector<join::JoinResultElement>& join_result,
                    double distance_threshold);
};
#include "t_balljoin_impl.h"
} // namespace join