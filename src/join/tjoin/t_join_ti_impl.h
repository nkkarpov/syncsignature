// The MIT License (MIT)
// Copyright (c) 2017 Thomas Huetter
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// Implements the TJoin tree similarity join.

#pragma once

template <typename Label, typename VerificationAlgorithm>
TJoinTI<Label, VerificationAlgorithm>::TJoinTI() {
    ld_ = label::LabelDictionary<Label>();
    pre_candidates_ = 0;
    sum_subproblem_counter_ = 0;
    number_of_labels_ = 0;
    il_lookups_ = 0;
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::execute_parallel_join(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<
        std::pair<int, std::vector<label_set_converter::LabelSetElement>>>&
        sets_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    double distance_threshold, int number_of_threads, int stop) {
    fprintf(stderr, "start TJoin\n");
    auto z = current_time();
    convert_trees_to_sets(trees_collection, sets_collection);
    retrieve_candidates(sets_collection, candidates, distance_threshold);
    fprintf(stderr, "number of threads = %d\n", number_of_threads);
    auto stamp = current_time();
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(begin(candidates), end(candidates),
                 std::default_random_engine(seed));
    assert(number_of_threads >= 1);
    if (stop){
        return;
    }
    {
        fprintf(stderr, "number of candidates %d\n", (int)candidates.size());
        fprintf(stderr, "join time = %lfs\n", get_time(z, stamp));
        auto f = [&](int shift) {
            auto ld = label::LabelDictionary<Label>();
            typename VerificationAlgorithm::AlgsCostModel cm(ld);
            VerificationAlgorithm ted_algorithm(cm);
            ted_ub::LGMTreeIndex<typename VerificationAlgorithm::AlgsCostModel>
                lgm_algorithm(cm);

            fprintf(stderr, "start %d\n", shift);

            std::vector<std::pair<int, int>> output;
            for (size_t i = shift; i < candidates.size();
                 i += number_of_threads) {
                ld.clear();
                node::TreeIndexLGM ti_1;
                node::TreeIndexLGM ti_2;
                auto id1 = candidates[i].first;
                auto id2 = candidates[i].second;
                node::index_tree(ti_1, trees_collection[id1], ld, cm);
                node::index_tree(ti_2, trees_collection[id2], ld, cm);
                double ub_value =
                    lgm_algorithm.ted_k(ti_1, ti_2, distance_threshold);
                if (ub_value <= distance_threshold) {
                    output.emplace_back(id1, id2);
                    continue;
                }
                typename VerificationAlgorithm::AlgsTreeIndex ti_1v;
                typename VerificationAlgorithm::AlgsTreeIndex ti_2v;
                ld.clear();
                node::index_tree(ti_1v, trees_collection[id1], ld, cm);
                node::index_tree(ti_2v, trees_collection[id2], ld, cm);
                double ted_value =
                    ted_algorithm.ted_k(ti_1v, ti_2v, distance_threshold);
                if (ted_value <= distance_threshold)
                    output.emplace_back(id1, id2);
            }
            return output;
        };
        std::vector<std::future<std::vector<std::pair<int, int>>>> tasks(
            number_of_threads);
        for (int i = 0; i < number_of_threads; ++i) {
            tasks[i] = std::async(std::launch::async, f, i);
        }

        for (int i = 0; i < number_of_threads; ++i) {
            auto output = tasks[i].get();
            for (auto& x : output) {
                join_result.emplace_back(x.first, x.second, 0);
            }
        }
    }
    auto njoin = join_result.size();
    auto d = current_time();
    printf("%s,%lf,%lf,%d,%lf,%lf,%lf,%lf,%lu,%d\n", "TJoin", 0.0, 0.4,
           (int)distance_threshold, 0.0, 0.0, 0.0, get_time(z, d), njoin,
           number_of_threads);
    fflush(stdout);
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::convert_trees_to_sets(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<
        std::pair<int, std::vector<label_set_converter::LabelSetElement>>>&
        sets_collection) {

    // Convert trees to sets and get the result.
    label_set_converter::Converter<Label> lsc;
    lsc.assignFrequencyIdentifiers(trees_collection, sets_collection);
    number_of_labels_ = lsc.get_number_of_labels();
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::retrieve_candidates(
    std::vector<
        std::pair<int, std::vector<label_set_converter::LabelSetElement>>>&
        sets_collection,
    std::vector<std::pair<int, int>>& candidates,
    const double distance_threshold) {

    // Initialize candidate index.
    candidate_index::CandidateIndex c_index;

    // Retrieve candidates from the candidate index.
    c_index.lookup(sets_collection, candidates, number_of_labels_,
                   distance_threshold);

    // Copy the number of pre-candidates.
    pre_candidates_ = c_index.get_number_of_pre_candidates();
    // Copy the number of inverted list lookups.
    il_lookups_ = c_index.get_number_of_il_lookups();
}

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::upperbound(
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
        ld_.clear();
        node::TreeIndexLGM ti_1;
        node::TreeIndexLGM ti_2;
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

template <typename Label, typename VerificationAlgorithm>
void TJoinTI<Label, VerificationAlgorithm>::verify_candidates(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    typename VerificationAlgorithm::AlgsCostModel cm(ld_);
    VerificationAlgorithm ted_algorithm(cm);
    //  printf("candidates for verify %lu\n", candidates.size());
    // Verify each pair in the candidate set
    for (const auto& pair : candidates) {
        typename VerificationAlgorithm::AlgsTreeIndex ti_1;
        typename VerificationAlgorithm::AlgsTreeIndex ti_2;
        ld_.clear();
        node::index_tree(ti_1, trees_collection[pair.first], ld_, cm);
        node::index_tree(ti_2, trees_collection[pair.second], ld_, cm);
        double ted_value = ted_algorithm.ted_k(ti_1, ti_2, distance_threshold);
        if (ted_value <= distance_threshold)
            join_result.emplace_back(pair.first, pair.second, ted_value);

        // Sum up all number of subproblems
        sum_subproblem_counter_ += ted_algorithm.get_subproblem_count();
    }
}

template <typename Label, typename VerificationAlgorithm>
long long int
TJoinTI<Label, VerificationAlgorithm>::get_number_of_pre_candidates() const {
    return pre_candidates_;
}

template <typename Label, typename VerificationAlgorithm>
long long int
TJoinTI<Label, VerificationAlgorithm>::get_subproblem_count() const {
    return sum_subproblem_counter_;
}

template <typename Label, typename VerificationAlgorithm>
long long int
TJoinTI<Label, VerificationAlgorithm>::get_number_of_il_lookups() const {
    return il_lookups_;
}
