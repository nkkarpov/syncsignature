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
#include <random>
#include <tang_join_ti.h>
#include <vector>

int main(int argc, char** argv) {
    //    assert(argc == 8 || argc == 9);
    //    assert(argc == 10);
    std::string input_file_name = std::string(argv[1]);
    int distance_threshold = atoi(argv[3]);
    int id = atoi(argv[2]);
    double fraction = atof(argv[4]);
    fprintf(stderr, "fraction = %lf\n", fraction);
    int R = atoi(argv[5]);
    double sim = atof(argv[6]);
    int seed = atoi(argv[7]);
    int number_of_threads = atoi(argv[8]);
    int down = atoi(argv[9]);

    int W = 1;

    //    int up = atoi(argv[10]);
    /*if (argc == 9) {
        number_of_threads = atoi(argv[8]);
    }*/
    //    assert(argc == 11);
    fprintf(stderr, "number of threads = %d\n", number_of_threads);
    std::cerr << std::mt19937::max() << std::endl;
    using Label = label::StringLabel;
    using CostModel = cost_model::UnitCostModelLD<Label>;
    using LabelDictionary = label::LabelDictionary<Label>;
    LabelDictionary ld;
    std::vector<node::Node<Label>> trees_collection;
    parser::BracketNotationParser bnp;
    bnp.parse_collection(trees_collection, input_file_name, down, (int)1e9);
    std::vector<int> borders{ 100, 200, 400, 600, 800, 1000, 2000, 3000 };
    fprintf(stderr, "read is done (%lu)\n", trees_collection.size());
    fprintf(stderr, "read is done (%lu)\n", trees_collection.size());
    std::vector<int> sizes(trees_collection.size());
    for (auto i = 0; i < trees_collection.size(); ++i)
        sizes[i] = trees_collection[i].get_tree_size();
    if (id == 8) {
        printf("size = %d\n", sizes.size());
        printf("minsize = %u\n", sizes[0]);
        double ave = 0;
        for (int i = 0; i < sizes.size(); ++i) {
            ave += sizes[i];
        }
        ave /= sizes.size();
        printf("ave = %lf\n", ave);
    }
    fprintf(stderr, "id = %d\n", id);
    srand(time(nullptr));
    if (id == 0) {
        join::TMinJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>> alg;
        std::vector<std::pair<int, int>> alg_candidates;
        std::vector<join::JoinResultElement> alg_join_result;
        for (auto i = 0u; i < borders.size(); i++) {
            auto s = std::lower_bound(sizes.begin(), sizes.end(), borders[i]) -
                     sizes.begin();
            fprintf(stderr, "s = %ld\n", s);
            std::vector<node::Node<Label>> tmp(trees_collection.begin() + s,
                                               trees_collection.end());
            for (auto iter = 0; iter < R; ++iter) {
                printf("%d,%d,", borders[i], (int)1e5);
                alg.execute_parallel_join(
                    tmp, alg_candidates, alg_join_result, fraction, sim,
                    (double)distance_threshold, rand(), number_of_threads, W);
                alg_candidates.clear();
                alg_join_result.clear();
            }
        }
    } else if (id == 1) {
        std::vector<std::pair<int, int>> candidates;
        std::vector<join::JoinResultElement> join_result;
        std::vector<node::BinaryNode<Label>> binary_trees_collection;
        join::TBallJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
            algpart;
        for (auto i = 0u; i < borders.size(); i++) {
            auto s = std::lower_bound(sizes.begin(), sizes.end(), borders[i]) -
                     sizes.begin();
            std::vector<node::Node<Label>> tmp(trees_collection.begin() + s,
                                               trees_collection.end());
            for (int iter = 0; iter < R; iter++) {
                printf("%d,%d,", borders[i], (int)1e5);
                algpart.execute_parallel_join(tmp, binary_trees_collection,
                                              candidates, join_result, fraction,
                                              sim, distance_threshold, rand(),
                                              number_of_threads);

                binary_trees_collection.clear();
                candidates.clear();
                join_result.clear();
            }
        }
    } else if (id == 2) {
        std::vector<
            std::pair<int, std::vector<label_set_converter::LabelSetElement>>>
            sets_collection;
        std::vector<std::pair<int, int>> candidates;
        std::vector<join::JoinResultElement> join_result;

        join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
            ted_join_algorithm;
        for (auto i = 0u; i < borders.size(); i++) {
            auto s = std::lower_bound(sizes.begin(), sizes.end(), borders[i]) -
                     sizes.begin();
            std::vector<node::Node<Label>> tmp(trees_collection.begin() + s,
                                               trees_collection.end());
            for (int iter = 0; iter < R; iter++) {
                printf("%d,%d,", borders[i], (int)1e5);
                ted_join_algorithm.execute_parallel_join(
                    tmp, sets_collection, candidates, join_result,
                    (double)distance_threshold, number_of_threads);
                sets_collection.clear();
                candidates.clear();
                join_result.clear();
            }
        }
    } else if (id == 3) {

        //        auto s =
        //            std::lower_bound(sizes.begin(), sizes.end(), 0) -
        //            sizes.begin();
        //        auto t = std::upper_bound(sizes.begin(), sizes.end(),
        //                                  1000 + distance_threshold) -
        //                 sizes.begin();
        fprintf(stderr, "s = %d\n", (int)0);
        {
            join::TMinJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                alg;
            std::vector<std::pair<int, int>> alg_candidates;
            std::vector<join::JoinResultElement> alg_join_result;
            for (auto iter = 0; iter < R; ++iter) {
                printf("%d,%d,", down, (int)1e5);
                alg.execute_parallel_join(trees_collection, alg_candidates,
                                          alg_join_result, fraction, sim,
                                          (double)distance_threshold, rand(),
                                          number_of_threads, W);
                alg_candidates.clear();
                alg_join_result.clear();
            }
        }
    } else if (id == 4) {

        {
            std::vector<std::pair<int, int>> candidates;
            std::vector<join::JoinResultElement> join_result;
            std::vector<node::BinaryNode<Label>> binary_trees_collection;
            join::TBallJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                algpart;
            for (int iter = 0; iter < R; iter++) {
                printf("%d,%d,", down, (int)1e5);
                algpart.execute_parallel_join(
                    trees_collection, binary_trees_collection, candidates,
                    join_result, fraction, sim, distance_threshold, rand(),
                    number_of_threads);

                binary_trees_collection.clear();
                candidates.clear();
                join_result.clear();
            }
        }
        /*{
            std::vector<std::pair<
                int, std::vector<label_set_converter::LabelSetElement>>>
                sets_collection;
            std::vector<std::pair<int, int>> candidates;
            std::vector<join::JoinResultElement> join_result;

            join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                ted_join_algorithm;
            std::vector<node::Node<Label>> tmp(trees_collection.begin(),
                                               trees_collection.begin() + t);
            for (int iter = 0; iter < R; iter++) {
                printf("%d,%d,", 0, 1000 + distance_threshold);
                ted_join_algorithm.execute_parallel_join(
                    tmp, sets_collection, candidates, join_result,
                    (double)distance_threshold, number_of_threads);
                sets_collection.clear();
                candidates.clear();
                join_result.clear();
            }
        }*/
    } else if (id == 5) {
        std::vector<
            std::pair<int, std::vector<label_set_converter::LabelSetElement>>>
            sets_collection;
        std::vector<std::pair<int, int>> candidates;
        std::vector<join::JoinResultElement> join_result;

        join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
            ted_join_algorithm;
        for (int iter = 0; iter < R; iter++) {
            printf("%d,%d,", down, (int)1e5);
            ted_join_algorithm.execute_parallel_join(
                trees_collection, sets_collection, candidates, join_result,
                (double)distance_threshold, number_of_threads);
            sets_collection.clear();
            candidates.clear();
            join_result.clear();
        }
    } else if (id == 6) {
        {
            std::vector<std::pair<
                int, std::vector<label_set_converter::LabelSetElement>>>
                sets_collection;
            std::vector<std::pair<int, int>> candidates;
            std::vector<join::JoinResultElement> join_result;

            join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                ted_join_algorithm;
            for (int iter = 0; iter < R; iter++) {
                printf("%d,%d,", 0, (int)(1e5));
                ted_join_algorithm.execute_parallel_join(
                    trees_collection, sets_collection, candidates, join_result,
                    (double)distance_threshold, number_of_threads);
                sets_collection.clear();
                candidates.clear();
                join_result.clear();
            }
        }
        {
            auto t = std::upper_bound(sizes.begin(), sizes.end(),
                                      100 + distance_threshold) -
                     sizes.begin();
            std::vector<std::pair<
                int, std::vector<label_set_converter::LabelSetElement>>>
                sets_collection;
            std::vector<std::pair<int, int>> candidates;
            std::vector<join::JoinResultElement> join_result;

            join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                ted_join_algorithm;
            std::vector<node::Node<Label>> tmp(trees_collection.begin(),
                                               trees_collection.begin() + t);
            for (int iter = 0; iter < R; iter++) {
                printf("%d,%d,", 0, 100);
                ted_join_algorithm.execute_parallel_join(
                    tmp, sets_collection, candidates, join_result,
                    (double)distance_threshold, number_of_threads);
                sets_collection.clear();
                candidates.clear();
                join_result.clear();
            }
        }
    } else if (id == 7) {
        fprintf(stderr, "id = %d\n", id);
        std::map<std::string, double> cnt;
        for (auto& x : trees_collection) {
            std::function<void(const node::Node<Label>&)> dfs;
            std::map<std::string, int> q;
            int t = 0;
            dfs = [&](const node::Node<Label>& s) {
                q[s.label().to_string()]++;
                t++;
                for (const auto& c : s.get_children()) {
                    dfs(c);
                }
            };
            dfs(x);
            for (const auto& s : q) {
                cnt[s.first] = std::max(cnt[s.first], 1.0 * s.second / t);
            }
        }
        std::vector<std::pair<double, std::string>> arr;
        arr.reserve(cnt.size());
        for (const auto& t : cnt) {
            arr.emplace_back(t.second, t.first);
        }
        sort(begin(arr), end(arr));
        //        reverse(begin(arr), end(arr));
        for (auto x : arr) {
            fprintf(stderr, "%s -> %lf\n", x.second.c_str(), x.first);
        }
    } else if (id == 10) {
        {
            join::TMinJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                alg;
            std::vector<std::pair<int, int>> alg_candidates;
            std::vector<join::JoinResultElement> alg_join_result;
            for (auto iter = 0; iter < 1; ++iter) {
                //                printf("%d,%d,", down, (int)1e5);
                alg.execute_parallel_join(trees_collection, alg_candidates,
                                          alg_join_result, fraction, sim,
                                          (double)distance_threshold, rand(),
                                          number_of_threads, W);
                printf("%lu %lu\n", trees_collection.size(),
                       alg_join_result.size());
                for (int i = 0; i < alg_join_result.size(); i++) {
                    printf("%d %d\n", alg_join_result[i].tree_id_1,
                           alg_join_result[i].tree_id_2);
                }
                alg_candidates.clear();
                alg_join_result.clear();
            }
        }
    } else if (id == 11) {
        {
            std::vector<std::pair<int, int>> candidates;
            std::vector<join::JoinResultElement> join_result;
            std::vector<node::BinaryNode<Label>> binary_trees_collection;
            join::TBallJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                algpart;
            for (int iter = 0; iter < 1; iter++) {
                //                printf("%d,%d,", down, (int)1e5);
                algpart.execute_parallel_join(
                    trees_collection, binary_trees_collection, candidates,
                    join_result, fraction, sim, distance_threshold, rand(),
                    number_of_threads);

                printf("%lu %lu\n", trees_collection.size(),
                       join_result.size());
                for (int i = 0; i < join_result.size(); i++) {
                    printf("%d %d\n", join_result[i].tree_id_1,
                           join_result[i].tree_id_2);
                }
                binary_trees_collection.clear();
                candidates.clear();
                join_result.clear();
            }
        }
    } else if (id == 12) {
        std::vector<
            std::pair<int, std::vector<label_set_converter::LabelSetElement>>>
            sets_collection;
        std::vector<std::pair<int, int>> candidates;
        std::vector<join::JoinResultElement> join_result;

        join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
            ted_join_algorithm;
        for (int iter = 0; iter < R; iter++) {
            //            printf("%d,%d,", down, (int)1e5);
            ted_join_algorithm.execute_parallel_join(
                trees_collection, sets_collection, candidates, join_result,
                (double)distance_threshold, number_of_threads);
            printf("%lu %lu\n", trees_collection.size(), join_result.size());
            for (int i = 0; i < join_result.size(); i++) {
                printf("%d %d\n", join_result[i].tree_id_1,
                       join_result[i].tree_id_2);
            }
            sets_collection.clear();
            candidates.clear();
            join_result.clear();
        }
    } else if (id == 13) {
        for (int iter = 0; iter < R; iter++) {
            join::TMinJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                alg;
            std::vector<std::pair<int, int>> alg_candidates;
            std::vector<join::JoinResultElement> alg_join_result;
            std::vector<std::pair<int, int>> tt;
            std::vector<int> result;
            for (auto q = 0; q < 10; ++q) {
                printf("%d,%d,", down, (int)1e5);
                alg.execute_parallel_join(trees_collection, alg_candidates,
                                          alg_join_result, fraction, sim,
                                          (double)distance_threshold, rand(),
                                          number_of_threads, W);
                alg_candidates.clear();
                for (auto x : alg_join_result) {
                    auto a = std::min(x.tree_id_1, x.tree_id_2);
                    auto b = std::max(x.tree_id_1, x.tree_id_2);
                    tt.emplace_back(a, b);
                }
                std::sort(tt.begin(), tt.end());
                tt.resize(std::unique(tt.begin(), tt.end()) - tt.begin());
                result.push_back((int)tt.size());
                alg_join_result.clear();
            }
            printf("EJoin,%d", distance_threshold);
            for (int i = 0; i < (int)result.size(); i++) {
                printf(",%d", result[i]);
            }
            printf("\n");
        }
        for (int iter = 0; iter < R; iter++) {
            std::vector<std::pair<int, int>> candidates;
            std::vector<join::JoinResultElement> join_result;
            std::vector<node::BinaryNode<Label>> binary_trees_collection;
            join::TBallJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                algpart;

            std::vector<std::pair<int, int>> tt;
            std::vector<int> result;
            for (auto q = 0; q < 10; ++q) {
                printf("%d,%d,", down, (int)1e5);
                algpart.execute_parallel_join(
                    trees_collection, binary_trees_collection, candidates,
                    join_result, fraction, sim, distance_threshold, rand(),
                    number_of_threads);

                binary_trees_collection.clear();
                candidates.clear();

                for (auto x : join_result) {
                    auto a = std::min(x.tree_id_1, x.tree_id_2);
                    auto b = std::max(x.tree_id_1, x.tree_id_2);
                    tt.emplace_back(a, b);
                }
                std::sort(tt.begin(), tt.end());
                tt.resize(std::unique(tt.begin(), tt.end()) - tt.begin());
                result.push_back((int)tt.size());
                join_result.clear();
            }
            printf("BJoin,%d", distance_threshold);
            for (int i = 0; i < (int)result.size(); i++) {
                printf(",%d", result[i]);
            }
            printf("\n");
        }
    }

    return 0;
}
/*
 * 5   : 32307902 717.355247
 * 10  : 36454496 758.652922
 * alg : 32307931
 * 5   :  2672167
 * 5   :  2668545
 * cand: 32393360
 * */
