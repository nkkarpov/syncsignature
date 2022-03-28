/*

//
// Created by nkarpov on 3/24/21.
//

#pragma once
#include "binary_tree_converter.h"
template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::convert_trees_to_binary_trees(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<node::BinaryNode<Label>>& binary_trees_collection) {

    // Convert trees to sets and get the result.
    binary_tree_converter::Converter<Label> btsc;
    btsc.convert(trees_collection, binary_trees_collection);
}

template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::execute_join(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<node::BinaryNode<Label>>& binary_trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result, int T, int t2bin,
    const double sim, const double distance_threshold) {

    auto a = current_time();

    printf("convert trees to binary trees start\n");
    if (t2bin) {
        convert_trees_to_binary_trees(trees_collection,
                                      binary_trees_collection);
    }
    printf("convert trees to binary trees end\n");
    std::vector<std::vector<size_t>> signature_collection;
    std::vector<std::vector<size_t>> tours_collection;
    convert_trees_to_tours(trees_collection, tours_collection);

    if (t2bin) {
        printf("partition_binary_tree start\n");
        partition_binary_tree(binary_trees_collection, signature_collection, T);
        printf("partition_binary_tree end\n");
    } else {
//        printf("partition_tree start\n");
        partition_tree(trees_collection, signature_collection, T);
//        printf("partition_tree end\n");
    }

    auto b = current_time();
    retrieve_candidates(tours_collection, signature_collection, candidates, sim,
                        distance_threshold);
    auto ncandidates = candidates.size();
    auto c = current_time();
    lowerbound(tours_collection, candidates, distance_threshold);
    upperbound(trees_collection, candidates, join_result, distance_threshold);
    auto d = current_time();
    auto nverificatons = candidates.size();
    verify_candidates(trees_collection, candidates, join_result,
                      distance_threshold);
    auto e = current_time();
    auto njoin = join_result.size();
    // name, sim, T / K, K, W, join time, number of candidates, filtering time,
    // number of verifications, verification time, output, total verification
    // time, total time
    printf(
        "TPartJoin, %lf, %lf, %d, %d, %lf, %lu, %lf, %lu, %lf, %lu, %lf, %lf\n",
        sim, T / distance_threshold, (int)(distance_threshold), 1,
        get_time(a, c), ncandidates, get_time(c, d), nverificatons,
        get_time(d, e), njoin, get_time(c, e), get_time(a, e));
    fflush(stdout);
}
template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::lowerbound(
    const std::vector<std::vector<size_t>>& tours_collection,
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
void TPartJoinTI<Label, VerificationAlgorithm>::convert_trees_to_tours(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::vector<size_t>>& tours_collection) {
    tours_collection.resize(trees_collection.size());
    for (auto i = 0u; i < trees_collection.size(); ++i) {
        tours_collection[i].reserve(trees_collection[i].get_tree_size() * 2);
    }
    for (auto i = 0u; i < trees_collection.size(); ++i) {
        int counter = 0;
        detour(trees_collection[i], tours_collection[i], counter);
    }
}

template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::retrieve_candidates(
    std::vector<std::vector<size_t>>& tours_collection,
    std::vector<std::vector<size_t>>& signatures_collection,
    std::vector<std::pair<int, int>>& candidates, double sim,
    double distance_threshold) {
    using candidate_t = std::pair<int, int>;
    struct frequency_t {
        candidate_t c;
        double f;
    };
//    printf("start join...\n");
    candidates.clear();
    std::unordered_map<size_t, std::vector<std::pair<int, size_t>>> table;
    auto start = current_time();
    std::vector<int> sz(signatures_collection.size());
    for (int i = 0; i < (int)tours_collection.size(); i++) {
        sz[i] = tours_collection[i].size() / 2;
    }
    std::vector<frequency_t> result;
    for (auto index = 0; index < (int)signatures_collection.size(); ++index) {
        std::vector<std::pair<int, int>> output;
        */
/*if (index % 10000 == 0) {
            auto t = current_time();
            printf("join %d in %lf\n", index, get_time(start, t));
            printf("size = %lu\n", result.size());
        }*//*

        auto& signature = signatures_collection[index];
        for (auto& part : signature) {
            std::vector<std::pair<int, int>> tmp;
            auto it = table.find(part);
            if (it != table.end()) {
                auto& vec = it->second;
                for (auto i = 0; i < (int)vec.size(); ++i) {
                    while (i < (int)vec.size() &&
                           sz[vec[i].first] + distance_threshold < sz[index]) {
                        swap(vec[i], vec.back());
                        vec.pop_back();
                    }
                    if (i >= (int)vec.size()) {
                        break;
                    }
                    auto candidate = std::pair(vec[i].first, index);
                    if (candidate.first > candidate.second) {
                        std::swap(candidate.first, candidate.second);
                    }
                    tmp.push_back(candidate);
                }
            }
            sort(tmp.begin(), tmp.end());
            tmp.resize(unique(tmp.begin(), tmp.end()) - tmp.begin());
            for (auto& item : tmp) {
                output.push_back(item);
            }
        }

        for (auto& part : signature) {
            table[part].emplace_back(index, part);
        }
        sort(output.begin(), output.end());
        for (size_t i = 0; i < output.size();) {
            auto j = i;
            while (j < output.size() && output[i] == output[j])
                j++;
            double den = (double)signatures_collection[output[i].first].size();
            den += (double)signatures_collection[output[i].second].size();
            result.push_back(frequency_t{ output[i], (double)(j - i) / den });
            i = j;
        }
    }
    for (size_t i = 0; i < result.size(); i++) {
        while (i <= result.size() && result[i].f < sim) {
            std::swap(result[i], result.back());
            result.pop_back();
        }
    }
    for (const auto& item : result) {
        candidates.push_back(item.c);
    }
}
template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::verify_candidates(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::pair<int, int>>& candidates,
    std::vector<join::JoinResultElement>& join_result,
    const double distance_threshold) {

    label::LabelDictionary<Label> ld;
    typename VerificationAlgorithm::AlgsCostModel cm(ld);
    VerificationAlgorithm ted_algorithm(cm);
    typename VerificationAlgorithm::AlgsTreeIndex ti_1;
    typename VerificationAlgorithm::AlgsTreeIndex ti_2;

    // Verify each pair in the candidate set
    for (const auto& pair : candidates) {
        node::index_tree(ti_1, trees_collection[pair.first], ld, cm);
        node::index_tree(ti_2, trees_collection[pair.second], ld, cm);
        double ted_value = ted_algorithm.ted_k(ti_1, ti_2, distance_threshold);
        if (ted_value <= distance_threshold)
            join_result.emplace_back(pair.first, pair.second, ted_value);
        // Sum up all number of subproblems
    }
    //    printf("\n");
}

template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::partition_tree(
    std::vector<node::Node<Label>>& trees_collection,
    std::vector<std::vector<size_t>>& signature_collection, int T, int p) {
    assert(0);
    signature_collection.resize(trees_collection.size());
    int number = 0;
    int z;
    auto start = current_time();
    std::vector<std::vector<int>> graph;
    std::vector<size_t> rank;
    std::vector<int> f;
    std::vector<int> g;
    std::function<void(node::Node<Label>&)> dfs;
    //    std::function<void(node::Node<Label>*)> dfs_id;

    dfs = [&](node::Node<Label>& node) {
        int id = node.get_id();
        rank[id] = hash(node.label().to_string(), hash_bool[0]);
        for (auto& x : node.children_) {
            auto child_id = x.get_id();
            graph[id].push_back(child_id);
            graph[child_id].push_back(id);
            dfs(x);
        }
    };
    std::function<int(int, int)> count;
    std::map<std::pair<int, int>, int> mcnt;
    count = [&](int x, int p) {
        using std::pair;
        int ret = 0;
        if (mcnt.count(pair(x, p))) {
            return mcnt[pair(x, p)];
        }
        if (f[x]) {
            mcnt[pair(x, p)] = ret;
            return ret;
        }
        ret += 1;
        for (auto& y : graph[x]) {
            if (y == p)
                continue;
            ret += count(y, x);
        }
        mcnt[pair(x, p)] = ret;
        return ret;
    };
    auto check = [&](int x) {
        if (graph[x].size() <= 1) {
            return 0;
        }
        int t = count(x, x);
        for (auto& y : graph[x]) {
            int q = count(y, x);
            t = std::min(q, t);
        }
        return (int)(t >= z);
    };
    std::function<void(int, int)> mark;
    mark = [&](int x, int p) {
        if (f[x]) {
            return;
        }
        g[x] = 0;
        for (auto& y : graph[x]) {
            if (y == p) {
                continue;
            }
            mark(y, x);
        }
    };
    std::vector<int> st;
    std::vector<int> stn;
    std::vector<size_t> parts;
    std::function<void(node::Node<Label>&)> split;
    split = [&](node::Node<Label>& node) {
        int id = node.get_id();
        st.push_back(rank[id]);
        //        str.push_back(id);
        if (f[id]) {
            stn.push_back((int)st.size() - 1);
        }
        for (auto &x : node.children_) {
            split(x);
        }
        st.push_back(rank[id]);
        //        str.push_back(id);
        if (f[id]) {
            parts.push_back(hash(begin(st) + stn.back(), end(st), p));
            while ((int)st.size() > (int)stn.back()) {
                st.pop_back();
                //                str.pop_back();
            }
            stn.pop_back();
        }
    };
    for (int i = 0; i < (int)trees_collection.size(); ++i) {
        mcnt.clear();
        parts.clear();
        auto& tree = trees_collection[i];
        number++;
        auto n = tree.get_tree_size();
        z = n / (T);
        rank.resize(n);
        graph.resize(n);
        f.resize(n);
        g.resize(n);
        for (int j = 0; j < n; j++) {
            graph[j].clear();
            f[j] = 0;
            g[j] = 1;
        }
        dfs(tree);
        std::vector<int> p(n);
        for (int j = 0; j < n; j++) {
//            printf("%d: ", j);
//            for (int k = 0; k < graph[j].size(); k++) printf("%d%c", graph[j][k], ",\n"[k + 1 == graph[j].size()]);
            p[j] = j;
        }
        return;
        sort(begin(p), end(p), [&](const auto& o1, const auto& o2) {
            return rank[o1] < rank[o2];
        });

        int anch = 0;
        for (int j = 0; j < n; j++) {
            auto x = p[j];
            if (!g[x]) {
                continue;
            }

            if (check(x)) {
                anch++;
                f[x] = 1;
                mcnt.clear();
            }
        }
        f[tree.get_id()] = 1;
        split(tree);
        assert(st.size() == 0);
        assert(stn.size() == 0);
        signature_collection[i] = parts;
        */
/*if (number % 10'000 == 0) {
            auto t = current_time();
            printf("split %d in %lf\n", number, get_time(start, t));
            printf("size %d, z = %d\n", n, z);
            printf("produce %d parts\n", (int)parts.size());
            printf("with %d anchors\n", anch);
            for (auto& part : parts)
                printf("%lu ", part);
            printf("\n");
        }*//*

        for (int j = 0; j < n; j++) {
            graph[j].clear();
        }
    }
    printf("done %d\n", number);
}

template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::partition_binary_tree(
    std::vector<node::BinaryNode<Label>>& binary_trees_collection,
    std::vector<std::vector<size_t>>& signature_collection, int T) {
    //    int max_size = binary_trees_collection.back().get_tree_size();
    //    printf("start split %d\n", T);
    signature_collection.resize(binary_trees_collection.size());
    int number = 0;
    int z = 0;
    auto start = current_time();
    std::vector<std::vector<int>> graph;
    std::vector<size_t> rank;
    std::vector<int> f;
    std::vector<int> g;
    std::function<void(node::BinaryNode<Label>*)> dfs;
    dfs = [&](node::BinaryNode<Label>* node) {
        int id = node->get_id();
        rank[id] = hash(node->label().to_string(), hash_bool[0]);
        if (node->has_left_child()) {
            auto left_id = node->get_left_child()->get_id();
            graph[id].push_back(left_id);
            graph[left_id].push_back(id);
            dfs(node->get_left_child());
        }
        if (node->has_right_child()) {
            auto right_id = node->get_right_child()->get_id();
            graph[id].push_back(right_id);
            graph[right_id].push_back(id);
            dfs(node->get_right_child());
        }
    };

    std::function<int(int, int)> count;

    std::map<std::pair<int, int>, int> mcnt;
    count = [&](int x, int p) {
        using std::pair;
        int ret = 0;
        if (mcnt.count(pair(x, p))) {
            return mcnt[pair(x, p)];
        }
        if (f[x]) {
            mcnt[pair(x, p)] = ret;
            return ret;
        }
        ret += 1;
        for (auto& y : graph[x]) {
            if (y == p)
                continue;
            ret += count(y, x);
        }
        mcnt[pair(x, p)] = ret;
        return ret;
    };

    auto check = [&](int x) {
        if (graph[x].size() <= 1) {
            return 0;
        }
//        int t = count(x, x);
        for (auto& y : graph[x]) {
            int q = count(y, x);
            if (q < z) return 0;
        }
        return 1;
    };
    std::function<void(int, int)> mark;
    mark = [&](int x, int p) {
        if (f[x]) {
            return;
        }
//        g[x] = 0;
        mcnt.erase(std::pair(x, p));
        mcnt.erase(std::pair(p, x));
        for (auto& y : graph[x]) {
            if (y == p) {
                continue;
            }
            mark(y, x);
        }
    };
    std::vector<int> st;
    std::vector<int> stn;
    std::vector<size_t> parts;
    std::function<void(node::BinaryNode<Label>*)> split;
    split = [&](node::BinaryNode<Label>* node) {
        int id = node->get_id();
        st.push_back(rank[id]);
        //        str.push_back(id);
        if (f[id]) {
            stn.push_back((int)st.size() - 1);
        }
        if (node->has_left_child())
            split(node->get_left_child());
        if (node->has_right_child())
            split(node->get_right_child());
        st.push_back(rank[id]);
        //        str.push_back(id);
        if (f[id]) {
            parts.push_back(hash(begin(st) + stn.back(), end(st)));
            while ((int)st.size() > (int)stn.back()) {
                st.pop_back();
                //                str.pop_back();
            }
            stn.pop_back();
        }
    };
    for (int i = 0; i < (int)binary_trees_collection.size(); ++i) {
        parts.clear();
        mcnt.clear();
        auto& tree = binary_trees_collection[i];
        number++;
        auto n = tree.get_tree_size();
        z = n / T;
//        printf("size %d, z = %d\n", n, z);
        rank.resize(n);
        graph.resize(n);
        f.resize(n);
        g.resize(n);
        for (int j = 0; j < n; j++) {
            graph[j].clear();
            f[j] = 0;
            g[j] = 1;
        }
        dfs(&tree);
        std::vector<int> p(n);
        for (int j = 0; j < n; j++) {
            p[j] = j;
        }
        sort(begin(p), end(p), [&](const auto& o1, const auto& o2) {
            return rank[o1] < rank[o2];
        });

        int anch = 0;
        for (int j = 0; j < n; j++) {
            auto x = p[j];
            if (!g[x]) {
                continue;
            }
            */
/*int temp = count(x, x);
            assert(temp >= T);
            if (temp < T) {
                mark(x, x);
                continue;
            }*//*

            if (check(x)) {
                anch++;
                mark(x, x);
                f[x] = 1;
//                mcnt.clear();
            }
        }
        f[tree.get_id()] = 1;
        split(&tree);
        assert(st.size() == 0);
        assert(stn.size() == 0);
        signature_collection[i] = parts;
        */
/*if (number % 1'000 == 0) {
            auto t = current_time();
            printf("split %d in %lf\n", number, get_time(start, t));
            printf("size %d, z = %d\n", n, z);
            printf("produce %d parts\n", (int)parts.size());
            printf("with %d anchors\n", anch);
            for (auto& part : parts)
                printf("%lu ", part);
            printf("\n");
        }*//*

        for (int j = 0; j < n; j++) {
            graph[j].clear();
        }
    }
//    printf("done %d\n", number);
}

template <typename Label, typename VerificationAlgorithm>
void TPartJoinTI<Label, VerificationAlgorithm>::upperbound(
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
    node::TreeIndexLGM ti_1;
    node::TreeIndexLGM ti_2;

    std::vector<std::pair<int, int>>::iterator it = candidates.begin();
    while (it != candidates.end()) {
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
}*/
