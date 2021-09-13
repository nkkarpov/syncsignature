//
// Created by nkarpov on 3/9/21.
//

#include "bracket_notation_parser.h"
#include "node.h"
#include "string_label.h"
#include "t_balljoin.h"
#include "t_join_ti.h"
#include "t_minjoin.h"
#include "touzet_baseline_tree_index.h"
#include "unit_cost_model.h"
#include <apted_tree_index.h>
#include <binary_node.h>
#include <cassert>
#include <tang_join_ti.h>
#include <vector>
#include <random>

int main(int argc, char** argv) {
    assert(argc == 8 || argc == 9);
    std::string input_file_name = std::string(argv[1]);
    int distance_threshold = atoi(argv[3]);
    int id = atoi(argv[2]);
    double fraction = atof(argv[4]);
    fprintf(stderr, "fraction = %lf\n", fraction);
    int R = atoi(argv[5]);
    double sim = atof(argv[6]);
    int seed = atoi(argv[7]);
    int number_of_threads = 0;
    if (argc == 9) {
        number_of_threads = atoi(argv[8]);
    }
//    [[maybe_unused]] std::mt19937 generator(seed);
    std::cerr << std::mt19937::max() << std::endl;
    using Label = label::StringLabel;
    using CostModel = cost_model::UnitCostModelLD<Label>;
    using LabelDictionary = label::LabelDictionary<Label>;
    LabelDictionary ld;
    std::vector<node::Node<Label>> trees_collection;
    parser::BracketNotationParser bnp;
    bnp.parse_collection(trees_collection, input_file_name, 1000, false);
    srand(seed);
    if (id == 0) {
        join::TMinJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>> alg;
        std::vector<std::pair<int, int>> alg_candidates;
        std::vector<join::JoinResultElement> alg_join_result;
        for (auto iter = 0; iter < R; ++iter) {
            if (number_of_threads == 0) {
                /*alg.execute_join(trees_collection, alg_candidates,
                                 alg_join_result, fraction, sim,
                                 (double)distance_threshold, rand());*/
            } else {
                alg.execute_parallel_join(
                    trees_collection, alg_candidates, alg_join_result, fraction,
                    sim, (double)distance_threshold, rand(), number_of_threads);
            }
            alg_candidates.clear();
            alg_join_result.clear();
        }
    } else if (id == 1) {
        std::vector<std::pair<int, int>> candidates;
        std::vector<join::JoinResultElement> join_result;
        std::vector<node::BinaryNode<Label>> binary_trees_collection;
        join::TBallJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
            algpart;
        for (int iter = 0; iter < R; iter++) {
            if (number_of_threads == 0) {
                /*algpart.execute_join(trees_collection, binary_trees_collection,
                                     candidates, join_result, fraction, sim,
                                     distance_threshold, rand());*/
            } else {
                algpart.execute_parallel_join(
                    trees_collection, binary_trees_collection, candidates,
                    join_result, fraction, sim, distance_threshold, rand(),
                    number_of_threads);
            }
            binary_trees_collection.clear();
            candidates.clear();
            join_result.clear();
        }
    } else if (id == 2) {
        std::vector<
            std::pair<int, std::vector<label_set_converter::LabelSetElement>>>
            sets_collection;
        std::vector<std::pair<int, int>> candidates;
        std::vector<join::JoinResultElement> join_result;

        join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
            ted_join_algorithm;
        for (int iter = 0; iter < R; iter++) {

            ted_join_algorithm.execute_parallel_join(
                trees_collection, sets_collection, candidates, join_result,
                (double)distance_threshold, number_of_threads);
            sets_collection.clear();
            candidates.clear();
            join_result.clear();
        }
    } else if (id == 3) {
        join::TangJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>> alg;
        std::vector<node::BinaryNode<Label>> binary_trees_collection;
        std::unordered_set<std::pair<int, int>, join::hashintegerpair>
            candidates;
        std::vector<join::JoinResultElement> join_result;
        for (int iter = 0; iter < R; iter++) {
            alg.execute_join(trees_collection, binary_trees_collection,
                             candidates, join_result,
                             (double)distance_threshold);
            binary_trees_collection.clear();
            candidates.clear();
            join_result.clear();
        }
    }
    // swiss : 122772
    // python: 35958
    // json: 40795
    return 0;
}