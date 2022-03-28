

//
// Created by nkarpov on 3/10/21.
//

#pragma once

template <typename Label, typename VerificationAlgorithm>
TMinJoinTI<Label, VerificationAlgorithm>::TMinJoinTI() {
    ld_ = label::LabelDictionary<Label>();
}

static const size_t BIG_PRIME = 1000000000000000003;

size_t mult(size_t a, size_t b) { return a * b; }
size_t hash(const std::string& str, size_t init, int p) {
    size_t result = init;
    for (int i = 0; i < (int)str.size(); i++) {
        result = mult(result, p);
        result += (size_t)str[str.size() - i - 1];
    }
    return result;
}

static size_t hash_bool[] = { 11162313925820027003u, 5808496959448994821u };

template <class T>
size_t hash(T first, T last, int p) {
    size_t result = 0;
    for (auto it = first; it != last; it++) {
        result = mult(result, p);
        result += *it;
    }
    return result;
}

struct part_t {
    double pos{};
    size_t hash{};
};

template <typename Label>
void detour(const node::Node<Label>& s,
            std::vector<std::pair<size_t, int>>& tour,
            std::vector<size_t>& preorder, int p) {
    int cnt = preorder.size();
    tour.push_back(
        std::make_pair(hash(s.label().to_string(), hash_bool[0], p), cnt));
    preorder.push_back(tour.back().first);
    for (auto& x : s.get_children()) {
        detour(x, tour, preorder, p);
    }
    tour.push_back(
        std::make_pair(hash(s.label().to_string(), hash_bool[1], p), cnt));
}
using signature_t = std::vector<part_t>;

void create_signature(const std::vector<size_t>& preorder,
                      const std::vector<std::pair<size_t, int>>& tour,
                      std::vector<part_t>& signature, const int z, const int W,
                      int p) {
    assert((int)tour.size() >= W);
//    assert(W == 1);
    assert(z >= 1);
    signature.clear();
    std::vector<size_t> q(tour.size() - W + 1);
    for (auto i = 0u; i < tour.size() - W + 1; ++i) {
//        q[i] = tour[i].first;
        q[i] = 0;
        for (int j = 0; j < W; j++) {
            q[i] *= int(1e9) + 7;
            q[i] += tour[i + j].first;
        }
    }
    std::vector<int> anchors;
//    anchors.push_back(0);
    auto len = (int)tour.size()- W + 1;
    for (auto i = z; i + z < len;) {
        int d = 1;
        for (; d <= z && q[i] < q[i - d] && q[i] < q[i + d]; ++d)
            ;
        if (d == z + 1 && i != 0)
            anchors.push_back(i);
        i += d;
    }
    /*if (anchors.back() != (int)len) {
        anchors.push_back((int)len);
    }*/
//    assert(!anchors.empty());
    double U = 1; // was 1.0
    for (auto i = 0u; i + 1 < anchors.size(); i++) {
        if ((double)anchors[i + 1] - (double)anchors[i] > U) {
            std::vector<std::pair<int, size_t>> tmp;
            for (int j = anchors[i]; j < anchors[i + 1]; ++j) {
                assert(tour[j].second < preorder.size());
                tmp.push_back(
                    std::make_pair(tour[j].second, preorder[tour[j].second]));
            }
            /*sort(begin(tmp), end(tmp));
            tmp.resize(
                unique(begin(tmp), end(tmp),
                       [](auto& a, auto& b) { return a.first == b.first; }) -
                begin(tmp));*/
            assert(tmp.size() > 0);
            std::vector<size_t> tt(tmp.size());
            for (int j = 0; j < (int)tmp.size(); ++j) {
                tt[j] = tmp[j].second;
            }
            signature.push_back(
                part_t{ anchors[i] / 2.0, hash(tt.begin(), tt.end(), p) });
        }
    }

//    assert(!signature.empty());
}

template <class T>
inline int slide(T first1, T last1, T first2, T last2) {
    return std::distance(first1,
                         std::mismatch(first1, last1, first2, last2).first);
}

template <class T>
int edit(T first1, T last1, T first2, T last2, int k) {
    auto len1 = std::distance(first1, last1);
    auto len2 = std::distance(first2, last2);
    if (abs(len1 - len2) > k) {
        return (k + 1);
    }
    if (len1 > len2) {
        return edit(first2, last2, first1, last1, k);
    }
    int* fp;
    int* fc;
    std::vector<int> fc_buf(2 * k + 1);
    std::vector<int> fp_buf(2 * k + 1);
    int h, dl, du, d;
    fc_buf[k] = fp_buf[k] = -1;

    for (h = 0; h <= k; h++) {

        if ((h & 1) == 0) {
            fc = &fc_buf[0] + k;
            fp = &fp_buf[0] + k;
        } else {
            fc = &fp_buf[0] + k;
            fp = &fc_buf[0] + k;
        }
        dl = -std::min(1 + ((k - (len2 - len1)) / 2), (long)h);
        du = std::min(1 + k / 2 + (len2 - len1), (long)h);

        fp[dl - 1] = fp[dl] = fp[du] = fp[du + 1] = -1;

        for (d = dl; d <= du; d++) {
            int r = std::max(std::max(fp[d - 1], fp[d] + 1), fp[d + 1] + 1);

            if ((r >= len1) || (r + d >= len2)) {
                fc[d] = r;
            } else {
                fc[d] = slide(first1 + r, last1, first2 + r + d, last2) + r;
            }

            if ((d == len2 - len1) && (fc[d] >= len1)) {
                return (h == -1) ? (k + 1) : h;
            }
        }
    }
    return (k + 1);
}
template <class T>
int edit(const T& a, const T& b, int k) {
    return edit(begin(a), end(a), begin(b), end(b), k);
}

template <typename Label, typename VerificationAlgorithm>
void TMinJoinTI<Label, VerificationAlgorithm>::execute_parallel_join(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result, double fraction,
    double tau, double distance_threshold, int p, int number_of_threads,
    int W) {
    std::vector<std::vector<std::pair<size_t, int>>> tours_collection(
        trees_collection.size());
    std::vector<std::vector<size_t>> preorder(trees_collection.size());
    if (p % 2 == 0)
        p ^= 1;
    fprintf(stderr, "EJoin parallel p = %d\n", p);
    std::cerr << "prime = " << BIG_PRIME << std::endl;
    fprintf(stderr, "W = %d\n", W);
    auto a = current_time();
    {
        auto partition = [&](int shift) {
            for (auto i = shift; i < (int)trees_collection.size();
                 i += number_of_threads) {
                detour(trees_collection[i], tours_collection[i], preorder[i],
                       p);
            }
            return;
        };

        std::vector<std::future<void>> tasks(number_of_threads);
        for (int index = 0; index < number_of_threads; index++) {
            tasks[index] = std::async(std::launch::async, partition, index);
        }
        for (auto index = 0; index < number_of_threads; index++) {
            tasks[index].get();
        }
    }

    std::vector<int> sizes(tours_collection.size());
    for (auto i = 0; i < (int)sizes.size(); ++i) {
        sizes[i] = (int)preorder[i].size();
        if (i) {
            assert(sizes[i] >= sizes[i - 1]);
        }
    }
    auto b = current_time();
    auto partition_time = get_time(a, b);
    auto join_time = 0.0;
    //    fprintf(stderr, "partition in %lf\n", get_time(a, b));
    std::vector<std::pair<int, int>> segments;
    std::vector<int> start;
    int K = (int)distance_threshold;
    int _shift = int(ceil(K / fraction));
    //    int tt = 0;
//    int val = factor * K;

    for (int current = trees_collection.empty() ? 0 : sizes[0];
         !sizes.empty() && current < sizes.back(); current += _shift) {
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

        if (segments[id].first + 1 >= segments[id].second)
            continue;
        auto f1 = current_time();
        int current = start[id];
        auto z =
            std::max((int)(ceil(fraction * current / distance_threshold)), 1);

        std::vector<signature_t> signatures_collection;
        signatures_collection.resize(tours_collection.size());
        fprintf(stderr, "id = %d, cur = %d, z = %d, [%d, %d)\n", id, current, z,
                segments[id].first, segments[id].second);
        {
            auto split = [&](int shift) {
                for (auto i = segments[id].first + shift;
                     i < segments[id].second; i += number_of_threads) {
                    create_signature(preorder[i], tours_collection[i],
                                     signatures_collection[i], z, W, p);
                }
            };

            std::vector<std::future<void>> tasks(number_of_threads);
            for (int index = 0; index < number_of_threads; index++) {
                tasks[index] = std::async(std::launch::async, split, index);
            }
            for (auto index = 0; index < number_of_threads; index++) {
                tasks[index].get();
            }
        }
        auto f2 = current_time();
        partition_time += get_time(f1, f2);
        fprintf(stderr, "signature done in %lf s\n", get_time(f1, f2));

        if (sizes[segments[id].first] <= -1) {
        /*    std::vector<std::pair<
                int, std::vector<label_set_converter::LabelSetElement>>>
                sets_collection;
            std::vector<std::pair<int, int>> tcandidates;
            std::vector<join::JoinResultElement> tjoin_result;
            using CostModel = cost_model::UnitCostModelLD<Label>;
            join::TJoinTI<Label, ted::TouzetBaselineTreeIndex<CostModel>>
                ted_join_algorithm;
            std::vector<node::Node<Label>> tmp;
            for (int i = segments[id].first; i < segments[id].second; i++) {
                tmp.push_back(trees_collection[i]);
            }
            ted_join_algorithm.execute_parallel_join(
                tmp, sets_collection, tcandidates, tjoin_result,
                (double)distance_threshold, number_of_threads, 1);
            sets_collection.clear();
            assert(tjoin_result.empty());
            for (auto x : tcandidates) {
                candidates.emplace_back(
                    std::min(x.first, x.second) + segments[id].first,
                    std::max(x.first, x.second) + segments[id].first);
            }
            fprintf(stderr, "number of candidates = %d\n",
                    (int)tcandidates.size());
            tcandidates.clear();
            tjoin_result.clear();
*/
        } else {
            unsigned long sz = 0;
            for (auto index = segments[id].first; index < segments[id].second;
                 ++index) {
                sz += signatures_collection[index].size();
            }
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
//                        fprintf(stderr, "skip len = %d\n", len);
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
        join_time += get_time(f2, f3);
        fprintf(stderr, "process buckets done in %lf s\n", get_time(f2, f3));
    }

    for (auto tt : candidates) {
        if (tt.first >= tt.second) {
            fprintf(stderr, "fuck!\n");
        }
        assert(tt.first < tt.second);
    }
    fprintf(stderr, "check done\n");
    sort(begin(candidates), end(candidates));
    candidates.resize(unique(begin(candidates), end(candidates)) -
                      begin(candidates));

    auto c = current_time();
    fprintf(stderr, "whole join done in %lf seconds\n", get_time(a, c));
    fprintf(stderr, "number of candidates = %d\n", (int)candidates.size());
    for (auto& tour : tours_collection) {
        for (auto& x : tour) {
            x.second = 0;
        }
    }
    {
        unsigned seed =
            std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(begin(candidates), end(candidates),
                     std::default_random_engine(seed));

        using candidate_t = std::pair<int, int>;
        ;
        std::vector<std::future<std::vector<candidate_t>>> tasks(
            number_of_threads);
        auto f = [&](int shift) {
            std::vector<std::pair<int, int>> output;
            auto ld = ld_;
            typename VerificationAlgorithm::AlgsCostModel cm(ld);
            ted_ub::LGMTreeIndex<typename VerificationAlgorithm::AlgsCostModel>
                lgm_algorithm(cm);
            VerificationAlgorithm ted_algorithm(cm);
            int small_tree = 0;
            for (size_t i = shift; i < candidates.size();
                 i += number_of_threads) {
                if (abs(sizes[candidates[i].first] -
                        sizes[candidates[i].second]) > distance_threshold) {
                    continue;
                }
                if (sizes[candidates[i].first] + sizes[candidates[i].second] <=
                    (int)distance_threshold) {
                    output.push_back(std::make_pair(candidates[i].first,
                                                    candidates[i].second));
                    small_tree++;
                    //                    fprintf(stderr, "(%d, %d)", (int)i,
                    //                    (int)output.size());
                    continue;
                }
                if (edit(preorder[candidates[i].first],
                         preorder[candidates[i].second],
                         (int)distance_threshold) > distance_threshold) {
                    //                    fprintf(stderr, "break %d\n", (int)i);
                    continue;
                }
                node::TreeIndexLGM ti_1;
                node::TreeIndexLGM ti_2;
                ld.clear();
                node::index_tree(ti_1, trees_collection[candidates[i].first],
                                 ld, cm);
                node::index_tree(ti_2, trees_collection[candidates[i].second],
                                 ld, cm);

                double ubted =
                    lgm_algorithm.ted_k(ti_1, ti_2, (int)distance_threshold);
                if (ubted <= distance_threshold) {
                    output.push_back(std::make_pair(candidates[i].first,
                                                    candidates[i].second));
                    continue;
                }
                typename VerificationAlgorithm::AlgsTreeIndex ti_1v;
                typename VerificationAlgorithm::AlgsTreeIndex ti_2v;
                ld.clear();
                node::index_tree(ti_1v, trees_collection[candidates[i].first],
                                 ld, cm);
                node::index_tree(ti_2v, trees_collection[candidates[i].second],
                                 ld, cm);
                double ted =
                    ted_algorithm.ted_k(ti_1v, ti_2v, (int)distance_threshold);
                if (ted > distance_threshold) {
                    continue;
                }
                output.push_back(
                    std::make_pair(candidates[i].first, candidates[i].second));
            }
            fprintf(stderr, "number of small trees %d in %d\n", small_tree,
                    shift);
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
    auto verification_time = get_time(c, d);
    auto total_time = get_time(a, d);
    auto output = (int)join_result.size();
    printf("%s,%lf,%lf,%d,%lf,%lf,%lf,%lf,%d,%d\n", "EJoin", tau, fraction,
           (int)distance_threshold, partition_time, join_time,
           verification_time, total_time, output, number_of_threads);
    fflush(stdout);
}

template <typename Label, typename VerificationAlgorithm>
void TMinJoinTI<Label, VerificationAlgorithm>::lowerbound(
    const std::vector<std::vector<std::pair<size_t, int>>>& tours_collection,
    std::vector<std::pair<int, int>>& candidates,
    const double distance_threshold) {
    for (auto i = 0u; i < candidates.size(); i++) {
        while (i < candidates.size() &&
               edit(tours_collection[candidates[i].first],
                    tours_collection[candidates[i].second],
                    2 * (int)distance_threshold) > 2 * distance_threshold) {
            std::swap(candidates[i], candidates.back());
            candidates.pop_back();
        }
    }
}

template <typename Label, typename VerificationAlgorithm>
void TMinJoinTI<Label, VerificationAlgorithm>::verify_candidates(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    typename VerificationAlgorithm::AlgsCostModel cm(ld_);
    VerificationAlgorithm ted_algorithm(cm);

    for (const auto& pair : candidates) {
        ld_.clear();
        typename VerificationAlgorithm::AlgsTreeIndex ti_1;
        typename VerificationAlgorithm::AlgsTreeIndex ti_2;
        node::index_tree(ti_1, trees_collection[pair.first], ld_, cm);
        node::index_tree(ti_2, trees_collection[pair.second], ld_, cm);
        double ted_value = ted_algorithm.ted_k(ti_1, ti_2, distance_threshold);
        if (ted_value <= distance_threshold)
            join_result.emplace_back(pair.first, pair.second, ted_value);
    }
}

template <typename Label, typename VerificationAlgorithm>
void TMinJoinTI<Label, VerificationAlgorithm>::convert_trees_to_tours(
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
void TMinJoinTI<Label, VerificationAlgorithm>::upperbound(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    typename VerificationAlgorithm::AlgsCostModel cm(ld_);
    ted_ub::LGMTreeIndex<typename VerificationAlgorithm::AlgsCostModel>
        lgm_algorithm(cm);
    // TODO: Index trees only once for LGM And Verification using a TreeIndex
    //       that is a superset of TreeIndexLGM and
    //       VerificationAlgorithm::AlgsTreeIndex.

    std::vector<std::pair<int, int>>::iterator it = candidates.begin();
    while (it != candidates.end()) {
        node::TreeIndexLGM ti_1;
        node::TreeIndexLGM ti_2;
        ld_.clear();
        node::index_tree(ti_1, trees_collection[it->first], ld_, cm);
        node::index_tree(ti_2, trees_collection[it->second], ld_, cm);
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