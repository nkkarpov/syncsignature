//
// Created by nkarpov on 3/10/21.
//

#pragma once
#include <cassert>
#include <future>
#include <thread>
#include <tang_join_ti.h>
#include "touzet_baseline_tree_index.h"
#include <set>
#include <random>

namespace join {
template <typename Label, typename VerificationAlgorithm>
class TMinJoinTI {
public:
    TMinJoinTI();
    /*void execute_join(std::vector<node::Node<Label>>& trees_collection,
                      std::vector<std::pair<int, int>>& candidates,
                      std::vector<join::JoinResultElement>& join_result,
                      double fraction, const double tau,
                      const double distance_threshold, int p);*/
    void
    execute_parallel_join(std::vector<node::Node<Label>>& trees_collection,
                          std::vector<std::pair<int, int>>& candidates,
                          std::vector<join::JoinResultElement>& join_result,
                          double fraction, double tau,
                          double distance_threshold, int p, int number_of_threads, int W);

    void verify_candidates(std::vector<node::Node<Label>>& trees_collection,
                           std::vector<std::pair<int, int>>& candidates,
                           std::vector<join::JoinResultElement>& join_result,
                           const double distance_threshold);

    void convert_trees_to_tours(
        std::vector<node::Node<Label>>& trees_collection,
        std::vector<std::vector<std::pair<size_t, int>>>& tours_collection,
        std::vector<std::vector<size_t>>& preorder, int p);
    void upperbound(std::vector<node::Node<Label>>& trees_collection,
                    std::vector<std::pair<int, int>>& candidates,
                    std::vector<join::JoinResultElement>& join_result,
                    const double distance_threshold);
    void lowerbound(const std::vector<std::vector<std::pair<size_t, int>>>&
                        tours_collection,
                    std::vector<std::pair<int, int>>& candidates,
                    const double distance_threshold);

    /*bool filter(const std::vector<size_t>& x, const std::vector<size_t>& y,
                const double distance_threshold);*/

private:
    label::LabelDictionary<Label> ld_;
};
#include "t_minjoin_impl.h"
} // namespace join
