//
// Created by nkarpov on 6/12/21.
//

#pragma once
#include "binary_tree_converter.h"

template <typename Label, typename VerificationAlgorithm>
void TBallJoinTI<Label, VerificationAlgorithm>::convert_trees_to_binary_trees(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<node::BinaryNode<Label>>& binary_trees_collection) {
    binary_tree_converter::Converter<Label> btsc;
    btsc.convert(trees_collection, binary_trees_collection);
}

template <typename Label, typename VerificationAlgorithm>
void TBallJoinTI<Label, VerificationAlgorithm>::execute_parallel_join(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<node::BinaryNode<Label>>& binary_trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result, double fraction,
    double tau, double distance_threshold, int p, int number_of_threads) {
    std::vector<std::vector<std::pair<size_t, int>>> tours_collection(
        trees_collection.size());
    std::vector<std::vector<size_t>> preorder(trees_collection.size());
    //    binary_trees_collection.reserve(trees_collection.size(),
    //    node::BinaryNode<Label>(""));
    binary_trees_collection.clear();
    for (size_t i = 0; i < trees_collection.size(); ++i) {
        binary_trees_collection.push_back(
            std::move(node::BinaryNode(trees_collection[i].label())));
    }
    if (p % 2 == 0)
        p ^= 1;
    fprintf(stderr, "BJoin parallel p = %d\n", p);
    using candidate_t = std::pair<int, int>;

    auto a = current_time();
    {

        auto partition = [&](int shift) {
            for (auto i = shift; i < (int)trees_collection.size();
                 i += number_of_threads) {

                detour(trees_collection[i], tours_collection[i], preorder[i],
                       p);
            }
            binary_tree_converter::Converter<Label> btsc;

            btsc.convert(trees_collection, binary_trees_collection, shift,
                         number_of_threads);
        };

        std::vector<std::future<void>> tasks(number_of_threads);
        for (int index = 0; index < number_of_threads; index++) {
            tasks[index] = std::async(std::launch::async, partition, index);
        }
        for (auto index = 0; index < number_of_threads; index++) {
            tasks[index].get();
        }
    }
    fprintf(stderr, "BJoin partition done\n");
    auto b = current_time();
    auto partition_time = get_time(a, b);
    auto join_time = 0.0;
    fprintf(stderr, "partition in %lf\n", get_time(a, b));

    std::vector<int> sizes(tours_collection.size());
    fprintf(stderr, "sz %d %d\n", (int)preorder.size(), (int)sizes.size());
    for (auto i = 0; i < (int)sizes.size(); ++i) {
        sizes[i] = (int)preorder[i].size();
    }

    std::vector<std::pair<int, int>> segments;
    std::vector<int> start;
    int K = (int)distance_threshold;
    int _shift = int(ceil(K / fraction));
    int val = std::max(10 * K, 70);
    //    _shift = 10000;
    for (int current = 0; !sizes.empty() && current < sizes.back(); current += _shift) {
        auto lower_key = current;
        auto upper_key = current + _shift + K;

        auto lower = std::lower_bound(begin(sizes), end(sizes), lower_key) -
                     begin(sizes);
        auto upper = std::upper_bound(begin(sizes), end(sizes), upper_key) -
                     begin(sizes);
        segments.emplace_back(lower, upper);
        start.push_back(lower_key);
    }
    for (int id = 0; id < (int)segments.size(); id++) {
        fprintf(stderr, "id = %d\n", id);
        if (segments[id].first + 1 >= segments[id].second)
            continue;

        auto f1 = current_time();
        int current = start[id];
        auto z =
            std::max((int)(ceil(fraction * current / distance_threshold)), 1);
        std::vector<signature_t> signatures_collection;
        fprintf(stderr, "id = %d, cur = %d, z = %d, [%d, %d)\n", id, current, z,
                segments[id].first, segments[id].second);
        signatures_collection.resize(tours_collection.size());
        {
            auto split = [&](int shift) {
                partition_binary_tree(binary_trees_collection, tours_collection,
                                      signatures_collection, segments[id].first,
                                      segments[id].second, z, tau, p, shift,
                                      number_of_threads);
            };
            fprintf(stderr, "split start\n");
            std::vector<std::future<void>> tasks(number_of_threads);
            for (int index = 0; index < number_of_threads; index++) {
                tasks[index] = std::async(std::launch::async, split, index);
            }
            for (auto index = 0; index < number_of_threads; index++) {
                tasks[index].get();
            }
        }
        fprintf(stderr, "split done\n");
        auto f2 = current_time();
        partition_time += get_time(f1, f2);
        fprintf(stderr, "signature done in %lf s\n", get_time(f1, f2));
         {
            unsigned long sz = 0;
            for (auto index = segments[id].first; index < segments[id].second;
                 ++index) {
                sz += signatures_collection[index].size();
            }
            fprintf(stderr, "sz = %lu\n", sz);
            std::vector<std::pair<size_t, std::pair<int, int>>> d;
            d.reserve(sz);
            for (auto index = segments[id].first; index < segments[id].second;
                 ++index) {
                auto& signature = signatures_collection[index];
                for (auto& part : signature) {
                    d.push_back(std::make_pair(
                        part.hash, std::make_pair(part.pos, index)));
                }
            }
            sort(begin(d), end(d));
            std::vector<std::pair<size_t, size_t>> intervals;
            for (size_t i = 0; i < d.size();) {
                size_t j = i;
                while (j < d.size() && d[i].first == d[j].first)
                    j++;
                if (j - i > 1) {
                    intervals.push_back(std::make_pair(i, j));
                }
                i = j;
            }
            std::vector<std::future<std::vector<std::pair<int, int>>>> tasks(
                number_of_threads);
            auto f = [&](int shift) {
                std::vector<std::pair<int, int>> output;
                for (size_t index = shift; index < intervals.size();
                     index += number_of_threads) {
                    auto len = intervals[index].second - intervals[index].first;
                    if (len * 10 > segments[id].second - segments[id].first && len > 100) {
                        //                            fprintf(stderr, "skip len = %d\n", len);
                        continue;
                    }
                    for (size_t i = intervals[index].first;
                         i < intervals[index].second; ++i) {
                        std::set<int> v;

                        for (size_t j = i + 1; j < intervals[index].second;
                             ++j) {
                            auto& o1 = d[i];
                            auto& o2 = d[j];
                            if (o2.second.second == o1.second.second)
                                continue;
                            if (abs(o2.second.first - o1.second.first) >
                                distance_threshold) {
                                break;
                            }
                            if (abs((int)tours_collection[o2.second.second]
                                        .size() -
                                    (int)tours_collection[o1.second.second]
                                        .size()) > 2 * distance_threshold) {
                                continue;
                            }
                            int ind1 = o1.second.second;
                            int ind2 = o2.second.second;
                            if (v.find(ind2) != v.end())
                                continue;
                            v.insert(ind2);
                            if (ind1 > ind2)
                                std::swap(ind1, ind2);
                            output.push_back(std::make_pair(ind1, ind2));
                        }
                    }
                }
                return output;
            };
            for (int i = 0; i < number_of_threads; ++i) {
                tasks[i] = std::async(std::launch::async, f, i);
            }

            std::vector<std::pair<int, int>> output;
            for (int i = 0; i < number_of_threads; ++i) {
                auto tmp = tasks[i].get();
                for (auto& x : tmp) {
                    output.push_back(x);
                }
            }
            sort(begin(output), end(output));
            for (size_t i = 0; i < output.size();) {
                auto j = i;
                while (j < output.size() && output[i] == output[j])
                    j++;
                if (j - i >= tau) {
                    candidates.push_back(output[i]);
                }
                i = j;
            }
        }
        auto f3 = current_time();
        fprintf(stderr, "process buckets done in %lf s\n", get_time(f2, f3));
        join_time += get_time(f2, f3);

    }
    sort(begin(candidates), end(candidates));
    candidates.resize(unique(begin(candidates), end(candidates)) -
                      begin(candidates));
    auto c = current_time();
    fprintf(stderr, "whole join done in %lf seconds\n", get_time(a, c));
    {
        unsigned seed =
            std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(begin(candidates), end(candidates),
                     std::default_random_engine(seed));
        std::vector<std::future<std::vector<candidate_t>>> tasks(
            number_of_threads);
        auto f = [&](int shift) {
            std::vector<std::pair<int, int>> output;
            label::LabelDictionary<Label> ld;
            typename VerificationAlgorithm::AlgsCostModel cm(ld);
            ted_ub::LGMTreeIndex<typename VerificationAlgorithm::AlgsCostModel>
                lgm_algorithm(cm);
            VerificationAlgorithm ted_algorithm(cm);

            for (size_t i = shift; i < candidates.size();
                 i += number_of_threads) {
                if (preorder[candidates[i].first].size() + preorder[candidates[i].second].size() <= distance_threshold) {
                    output.push_back(
                        std::make_pair(candidates[i].first, candidates[i].second));
                    continue;
                }
                if (edit(preorder[candidates[i].first],
                         preorder[candidates[i].second],
                         distance_threshold) > distance_threshold) {
                    continue;
                }

                node::TreeIndexLGM ti_1;
                node::TreeIndexLGM ti_2;
                node::index_tree(ti_1, trees_collection[candidates[i].first],
                                 ld, cm);
                node::index_tree(ti_2, trees_collection[candidates[i].second],
                                 ld, cm);

                if (lgm_algorithm.ted_k(ti_1, ti_2, (int)distance_threshold) <=
                    distance_threshold) {
                    output.push_back(std::make_pair(candidates[i].first,
                                                    candidates[i].second));
                    continue;
                }
                typename VerificationAlgorithm::AlgsTreeIndex ti_1v;
                typename VerificationAlgorithm::AlgsTreeIndex ti_2v;

                node::index_tree(ti_1v, trees_collection[candidates[i].first],
                                 ld, cm);
                node::index_tree(ti_2v, trees_collection[candidates[i].second],
                                 ld, cm);
                if (ted_algorithm.ted_k(ti_1v, ti_2v,
                                        (int)distance_threshold) <=
                    distance_threshold) {
                    output.push_back(std::make_pair(candidates[i].first,
                                                    candidates[i].second));
                    continue;
                }
            }
            return output;
        };
        for (int i = 0; i < number_of_threads; i++) {
            tasks[i] = std::async(std::launch::async, f, i);
        }
        for (int i = 0; i < number_of_threads; ++i) {
            auto tmp = tasks[i].get();
            for (auto& x : tmp) {
                join_result.push_back(JoinResultElement(x.first, x.second, 0));
            }
        }
    }

    auto d = current_time();
    fprintf(stderr, "alg done in %lf seconds\n", get_time(a, d));
//    auto join_time = get_time(a, c);
    auto verification_time = get_time(c, d);
    auto total_time = get_time(a, d);
    auto output = (int)join_result.size();
    printf("%s,%lf,%lf,%d,%lf,%lf,%lf,%lf,%d,%d\n", "BJoin", tau, fraction,
           (int)distance_threshold, partition_time, join_time, verification_time, total_time,
           output, number_of_threads);
    fflush(stdout);
}

template <typename Label, typename VerificationAlgorithm>
void TBallJoinTI<Label, VerificationAlgorithm>::lowerbound(
    const std::vector<std::vector<size_t>>& preorder,
    std::vector<std::pair<int, int>>& candidates,
    const double distance_threshold) {
    for (auto i = 0u; i < candidates.size(); i++) {
        while (i < candidates.size() &&
               edit(preorder[candidates[i].first],
                    preorder[candidates[i].second],
                    2 * (int)distance_threshold) > 2 * distance_threshold) {
            std::swap(candidates[i], candidates.back());
            candidates.pop_back();
        }
    }
}
template <typename Label, typename VerificationAlgorithm>
void TBallJoinTI<Label, VerificationAlgorithm>::convert_trees_to_tours(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::vector<std::pair<size_t, int>>>& tours_collection,
    std::vector<std::vector<size_t>>& preorder, int p) {
    tours_collection.resize(trees_collection.size());
    preorder.resize(trees_collection.size());
    for (auto i = 0u; i < trees_collection.size(); ++i) {
        detour(trees_collection[i], tours_collection[i], preorder[i], p);
    }
}

template <typename Label, typename VerificationAlgorithm>
void TBallJoinTI<Label, VerificationAlgorithm>::verify_candidates(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    label::LabelDictionary<Label> ld;
    typename VerificationAlgorithm::AlgsCostModel cm(ld);
    VerificationAlgorithm ted_algorithm(cm);

    for (const auto& pair : candidates) {
        ld.clear();
        typename VerificationAlgorithm::AlgsTreeIndex ti_1;
        typename VerificationAlgorithm::AlgsTreeIndex ti_2;
        node::index_tree(ti_1, trees_collection[pair.first], ld, cm);
        node::index_tree(ti_2, trees_collection[pair.second], ld, cm);
        double ted_value = ted_algorithm.ted_k(ti_1, ti_2, distance_threshold);
        if (ted_value <= distance_threshold)
            join_result.emplace_back(pair.first, pair.second, ted_value);
    }
}
size_t hash(const vector<size_t>& vec, size_t p) {
    size_t result = 0;
    for (auto it = vec.begin(); it != vec.end(); it++) {
        result = mult(result, p);
        result += (size_t)*it;
        //        if (result >= BIG_PRIME)
        //            result -= BIG_PRIME;
    }
    return result;
}
template <typename Label, typename VerificationAlgorithm>
void TBallJoinTI<Label, VerificationAlgorithm>::partition_binary_tree(
    std::vector<node::BinaryNode<Label>>& binary_trees_collection,
    std::vector<std::vector<std::pair<size_t, int>>>& tours_collection,
    std::vector<signature_t>& signature_collection, int first, int last, int z,
    double tau, int p, int shift, int step) {
    fprintf(stderr, "partition shift = %d step = %d\n", shift, step);
    fprintf(stderr, "partition first = %d last = %d\n", first, last);
    fprintf(stderr, "partition z = %d tau = %lf\n", z, tau);
    int number = 0;
    std::vector<std::vector<int>> graph;
    std::vector<size_t> rank;
    //    z *= 2;
    int tt;
    std::function<void(node::BinaryNode<Label>*)> dfs;
    dfs = [&](node::BinaryNode<Label>* node) {
        assert(node != nullptr);
        int id = node->get_id();
        /*if (tt == 34250) {
            fprintf(stderr, "id = %d\n", id);
        }*/

        rank[id] = hash(node->label().to_string(), hash_bool[0], p);
        if (node->has_left_child()) {
            //            assert(node->get_left_child() != nullptr);
            auto left_id = node->get_left_child()->get_id();
            //            assert(left_id != 0);
            graph[id].push_back(left_id);
            graph[left_id].push_back(id);
            dfs(node->get_left_child());
        }
        if (node->has_right_child()) {
            //            assert(node->get_right_child() != nullptr);
            auto right_id = node->get_right_child()->get_id();
            //            assert(right_id != 0);
            graph[id].push_back(right_id);
            graph[right_id].push_back(id);
            dfs(node->get_right_child());
        }
    };

    std::function<int(int)> check;
    vector<int> v;
    part_t part;
    int version = 1;
    vector<int> d;
    check = [&](int s) {
        std::queue<int> q;
        vector<std::pair<int, size_t>> t;
        int distance = 0;
        int cnt = 0;
        q.push(s);
        d[s] = 0;
        v[s] = version;
        double value = (double)s; //????
        while (!q.empty()) {
            auto x = q.front();
            q.pop();
            if (distance != d[x] && cnt > z) {
                break;
            }
            if (x != s && rank[x] < rank[s]) {
                return 0;
            }
            t.push_back(std::make_pair(x, rank[x]));
            cnt++;
            distance = d[x];
            for (auto& y : graph[x]) {
                if (v[y] != version) {
                    v[y] = version;
                    d[y] = d[x] + 1;
                    q.push(y);
                }
            }
        }
        sort(begin(t), end(t));
        vector<size_t> tmp(t.size());
        for (int i = 0; i < (int)t.size(); ++i) {
            tmp[i] = t[i].second;
        }
        part.hash = hash(tmp, p);
        part.pos = value;
        return 1;
    };
    assert(last <= signature_collection.size());
//    assert(last <= tours_collection.size() / 2);
    for (int i = first + shift; i < last; i += step) {
        tt = i;
        signature_collection[i].clear();
        version = 1;
        number++;
        auto n = (int)tours_collection[i].size() / 2;

        rank.resize(n);
        v.resize(n);
        std::fill(begin(v), end(v), 0);
        d.resize(n);
        graph.resize(n);
        dfs(&binary_trees_collection[i]);
        vector<int> q(n);
        for (int j = 0; j < (int)q.size(); ++j) {
            q[j] = j;
        }
        sort(begin(q), end(q), [&](const auto& o1, const auto& o2) {
            return rank[o1] < rank[o2];
        });
        for (auto& s : q) {
            version++;
            if (check(s)) {
                signature_collection[i].push_back(part);
                if ((int)signature_collection[i].size() > 5 * tau)
                    break;
            }
        }
        for (int j = 0; j < n; j++) {
            graph[j].clear();
        }
    }
}

template <typename Label, typename VerificationAlgorithm>
void TBallJoinTI<Label, VerificationAlgorithm>::upperbound(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {
    label::LabelDictionary<Label> ld;
    typename VerificationAlgorithm::AlgsCostModel cm(ld);
    ted_ub::LGMTreeIndex<typename VerificationAlgorithm::AlgsCostModel>
        lgm_algorithm(cm);
    // TODO: Index trees only once for LGM And Verification using a TreeIndex
    //       that is a superset of TreeIndexLGM and
    //       VerificationAlgorithm::AlgsTreeIndex.

    std::vector<std::pair<int, int>>::iterator it = candidates.begin();
    while (it != candidates.end()) {
        node::TreeIndexLGM ti_1;
        node::TreeIndexLGM ti_2;
        ld.clear();
        node::index_tree(ti_1, trees_collection[it->first], ld, cm);
        node::index_tree(ti_2, trees_collection[it->second], ld, cm);
        double ub_value = lgm_algorithm.ted_k(ti_1, ti_2, distance_threshold);
        if (ub_value <= distance_threshold) {
            join_result.emplace_back(it->first, it->second, ub_value);
            *it = candidates.back();
            candidates.pop_back();
        } else {
            ++it;
        }
    }
}