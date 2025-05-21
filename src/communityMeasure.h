#ifndef CommunityMeasure_hpp
#define CommunityMeasure_hpp

#include <cstdint>
#include <memory>
#include <vector>
#include <numeric>
#include <limits>
#include <set>
#include <algorithm>
#include "eccezioni.h"
#include "matriciCommunity.h"
#include "grafi.h"

// ***********************************************
// ***********************************************
// ***********************************************

class CommunityMeasure {
public:
    virtual std::vector<std::pair<double, double>> globalValue(std::vector<std::uint_fast64_t>&, std::shared_ptr<std::vector<double>>) const = 0;
    virtual std::pair<double, double> localValue(std::set<std::uint_fast64_t>& x, std::shared_ptr<std::vector<double>>) const = 0;
    virtual std::pair<double, double> localValue(std::vector<bool>& x, std::shared_ptr<std::vector<double>>) const = 0;
    virtual ~CommunityMeasure() {};
    // questi metodi sono utilizzati in LouvainMethod
    virtual void LMStart(const UGraph& grafo,
                         const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) = 0;
    virtual double LMAt(std::uint_fast64_t node, std::uint_fast64_t cluster_from, std::uint_fast64_t cluster_to,
                        std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) const = 0;
    virtual void LMUpdate(std::uint_fast64_t node,
                          std::uint_fast64_t cluster_from,
                          std::uint_fast64_t cluster_to,
                          std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) = 0;
    virtual double LMValue() const = 0;
};


// ***********************************************
// ***********************************************
// ***********************************************

class PersistenceMeasure : public CommunityMeasure {  ///ex PersistenceMeasureV2
private:
    UGraph& __grafo;
    std::shared_ptr<std::vector<std::uint_fast64_t>> __in_arc;
    std::shared_ptr<std::vector<std::uint_fast64_t>> __out_arc;
    std::shared_ptr<std::vector<double>> __prs;
    std::vector<double> __pi; // serve per la strategia di Stefano che usa il risultato del duale
                                            // start inizializza a 0 questo vettore se non si risolve il duale

    using LMNODE = std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, double>>;
    using LMCLUSTER = std::map<std::uint_fast64_t, double>;

    LMCLUSTER __lm_clusters;
    LMNODE __lm_node;
public:
    PersistenceMeasure(UGraph& grafo) : __grafo(grafo) {
        __in_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __out_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __prs = std::make_shared<std::vector<double>>(__grafo.Size());
    }

    std::vector<std::pair<double, double>> globalValue(std::vector<std::uint_fast64_t>& membership, std::shared_ptr<std::vector<double>> pi) const {
        if (membership.size() != __grafo.Size()) {
            std::string err_str = "error: globalValue";
            throw_line(err_str);
        }

        auto num_node = __grafo.Size();

        std::uint_fast64_t maximum = *std::max_element(membership.begin(), membership.end());
        std::vector<std::pair<double, double>> local_values(maximum, std::make_pair(0.0, 0.0));

        for (std::uint_fast64_t i = 1; i <= maximum; i++){
            std::set<std::uint_fast64_t> x_1;
            for (std::uint_fast64_t j = 0; j < num_node; j++) {
                if (membership.at(j) == i) {
                    x_1.insert(j);
                }
            }
            if (x_1.size() > 0) {
                auto local_value = localValue(x_1, pi);
                local_values[i - 1] = local_value;
            } else {
                local_values[i - 1] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
            }
        }
        return local_values;
    }

    std::pair<double, double> localValue(std::set<std::uint_fast64_t>& x, std::shared_ptr<std::vector<double>> pi) const {
        double somma_interni_interni = 0.0;
        double somma_strength = 0.0;
        double somma_pi = 0.0;
        for (auto v1_it = x.begin(); v1_it != x.end(); v1_it++) {
            auto v1 = *v1_it;
            somma_strength += __grafo.Strength(v1);
            for (auto v2_it = std::next(v1_it); v2_it != x.end(); v2_it++) {
                auto v2 = *v2_it;
                auto peso = __grafo.Weight(v1, v2).second;
                somma_interni_interni += peso;
            }
            somma_pi += (pi == nullptr ? 0 : pi->at(v1));
        }
        double risultato_persistence = (2.0 * somma_interni_interni) / ((double) somma_strength);
        return std::make_pair(risultato_persistence, - somma_pi);
    }

    std::pair<double, double> localValue(std::vector<bool>& x, std::shared_ptr<std::vector<double>> pi) const {
        std::set<std::uint_fast64_t> xset;
        for (std::uint_fast64_t v = 0; v < x.size(); ++v) {
            if (x.at(v)) {
                xset.insert(v);
            }
        }
        return localValue(xset, pi);
    }

    // per LouvainMethod
    void LMStart(const UGraph& grafo,
                 const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) {
        __lm_clusters.clear();
        __lm_node.clear();
        for (auto it_big_node_1 = clusters.begin(); it_big_node_1 != clusters.end(); ++it_big_node_1) {
            auto cm_c1_iter = __lm_node.insert({it_big_node_1->first, std::map<std::uint_fast64_t, double>()});

            auto c1_pair = localValue(*(it_big_node_1->second), nullptr);
            auto c1_value = c1_pair.first + c1_pair.second;

            cm_c1_iter.first->second.insert({it_big_node_1->first, c1_value});
            __lm_clusters[it_big_node_1->first] = c1_value;

            for (auto it_big_node_2 = std::next(it_big_node_1); it_big_node_2 != clusters.end(); ++it_big_node_2) {

                auto c2_pair = localValue(*(it_big_node_2->second), nullptr);
                auto c2_value = c2_pair.first + c2_pair.second;

                std::set<std::uint_fast64_t> c3_set(it_big_node_1->second->begin(), it_big_node_1->second->end());
                c3_set.insert(it_big_node_2->second->begin(), it_big_node_2->second->end());

                auto c3_pair = localValue(c3_set, nullptr);
                auto c3_value = c3_pair.first + c3_pair.second;

                cm_c1_iter.first->second.insert({it_big_node_2->first, c3_value - c1_value - c2_value});
            }
        }
    }

    double LMAt(std::uint_fast64_t node,
                std::uint_fast64_t cluster_from,
                std::uint_fast64_t cluster_to,
                std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) const {

        //double at_end = it_modularity_from->second + it_modularity_to->second;
        double at_end = std::accumulate(std::begin(__lm_clusters),
                                        std::end(__lm_clusters),
                                        0.0,
                                        [] (double value, const LMCLUSTER::value_type& p){
                                            return value + p.second;
                                        });


        {
            auto it_clusters_from = clusters.find(cluster_from);
            for (auto it_node_in_cluster = it_clusters_from->second->begin() ; it_node_in_cluster != it_clusters_from->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                at_end -= valore;
            }
        }
        {
            auto it_clusters_to = clusters.find(cluster_to);
            at_end += __lm_node.at(node).at(node);
            for (auto it_node_in_cluster = it_clusters_to->second->begin() ; it_node_in_cluster != it_clusters_to->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                at_end += valore;
            }
        }
        return at_end;
    };

    void LMUpdate(std::uint_fast64_t node,
                  std::uint_fast64_t cluster_from,
                  std::uint_fast64_t cluster_to,
                  std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) {

        // da utilizzare prima di modificare clusters
        {
            auto it_modularity_from = __lm_clusters.find(cluster_from);
            auto it_clusters_from = clusters.find(cluster_from);
            for (auto it_node_in_cluster = it_clusters_from->second->begin() ; it_node_in_cluster != it_clusters_from->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                it_modularity_from->second -= valore;
            }
        }
        {
            auto it_modularity_to = __lm_clusters.find(cluster_to);
            auto it_clusters_to = clusters.find(cluster_to);
            it_modularity_to->second += __lm_node.at(node).at(node);
            for (auto it_node_in_cluster = it_clusters_to->second->begin() ; it_node_in_cluster != it_clusters_to->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                it_modularity_to->second += valore;
            }
        }
    };

    double LMValue() const {
        double result = std::accumulate(std::begin(__lm_clusters),
                                        std::end(__lm_clusters),
                                        0.0,
                                        [] (double value, const LMCLUSTER::value_type& p){
                                            return value + p.second;
                                        });
        return result;
    };


private:
    void BuildFrom(std::map<std::uint_fast64_t, std::set<std::uint_fast64_t>>& communities) {
        std::fill(__in_arc->begin(), __in_arc->end(), 0);
        std::fill(__out_arc->begin(), __out_arc->end(), 0);
        std::fill(__prs->begin(), __prs->end(), 0.0);

        for (std::uint_fast64_t row = 0; row  < __grafo.Size(); ++row) {
            if (communities.find(row) != communities.end()) {
                auto& s = communities.at(row);
                for (auto v1 : s) {
                    for (auto v2 : s) {
                        __in_arc->at(row) += __grafo.Weight(v1, v2).first;
                    }
                    for (std::uint_fast64_t v2 = 0; v2  < __grafo.Size(); ++v2) {
                        if (s.find(v2) == s.end()) {
                            __out_arc->at(row) += __grafo.Weight(v1, v2).first;
                        }
                    }
                }
                __in_arc->at(row) /= 2;
                __prs->at(row) = ((double) __in_arc->at(row)) / ((double) (__in_arc->at(row) + __out_arc->at(row)));
            }
        }
        return;
    }
};


// ***********************************************
// ***********************************************
// ***********************************************

class PersistenceModularityMeasure : public CommunityMeasure { // ex: PersistenceModularityMeasureV2
    // P_{\cal C}
protected:
    UGraph& __grafo;

    std::shared_ptr<std::vector<std::uint_fast64_t>> __in_arc;
    std::shared_ptr<std::vector<std::uint_fast64_t>> __out_arc;
    std::shared_ptr<std::vector<double>> __prs;
    std::vector<double> __pi;

    using LMNODE = std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, double>>;
    using LMCLUSTER = std::map<std::uint_fast64_t, double>;

    LMCLUSTER __lm_clusters;
    LMNODE __lm_node;
public:
    PersistenceModularityMeasure(UGraph& grafo) : __grafo(grafo) {
        __in_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __out_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __prs = std::make_shared<std::vector<double>>(__grafo.Size());
    }

    std::vector<std::pair<double, double>> globalValue(std::vector<std::uint_fast64_t>& membership, std::shared_ptr<std::vector<double>> pi) const {
        if (membership.size() != __grafo.Size()) {
            std::string err_str = "error: globalValue";
            throw_line(err_str);
        }

        auto num_node = __grafo.Size();

        std::uint_fast64_t maximum = *std::max_element(membership.begin(), membership.end());
        std::vector<std::pair<double, double>> local_values(maximum, std::make_pair(0.0, 0.0));

        for (std::uint_fast64_t i = 1; i <= maximum; i++){
            std::set<std::uint_fast64_t> x_1;
            //std::cout << "globalValue:" << i << std::endl;

            for (std::uint_fast64_t j = 0; j < num_node; j++) {
                if (membership.at(j) == i) {
                    //std::cout << j << ", ";
                    x_1.insert(j);
                }
            }
            //std::cout << "." << std::endl;

            if (x_1.size() > 0) {
                auto local_value = localValue(x_1, pi);
                local_values[i - 1] = local_value;
            } else {
                local_values[i - 1] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
            }
        }
        return local_values;
    }

    std::pair<double, double> localValue(std::set<std::uint_fast64_t>& x, std::shared_ptr<std::vector<double>> pi) const {
        double somma_interni_interni = 0.0;
        double somma_strength = 0.0;
        double somma_pi = 0.0;
        for (auto v1_it = x.begin(); v1_it != x.end(); v1_it++) {
            auto v1 = *v1_it;
            somma_strength += __grafo.Strength(v1);
            for (auto v2_it = std::next(v1_it); v2_it != x.end(); v2_it++) {
                auto v2 = *v2_it;
                auto peso = __grafo.Weight(v1, v2).second;
                somma_interni_interni += peso;
            }
            somma_pi += (pi == nullptr ? 0 : pi->at(v1));
        }

        //double risultato = ((2.0 * somma_interni_interni) / ((double) somma_strength)) - (((double) somma_strength) /  ((double) __grafo.Degree()));
        double risultato_persistence = (2.0 * somma_interni_interni) /
                                                ((double) somma_strength);

        double risultato_null = ((double) somma_strength) /  (__grafo.Degree());

        double risultato = risultato_persistence - risultato_null;



        //std::cout << "somma_interni_interni:" << somma_interni_interni << std::endl;
        //std::cout << "somma_strength:" << somma_strength << std::endl;
        //std::cout << "__grafo.Degree():" << __grafo.Degree() << std::endl;
        //std::cout << "risultato_persistence:" << risultato_persistence << std::endl;
        //std::cout << "risultato_null:" << risultato_null << std::endl;

        return std::make_pair(risultato, - somma_pi);
    }

    std::pair<double, double> localValue(std::vector<bool>& x, std::shared_ptr<std::vector<double>> pi) const {
        std::set<std::uint_fast64_t> xset;
        for (std::uint_fast64_t v = 0; v < x.size(); ++v) {
            if (x.at(v)) {
                xset.insert(v);
            }
        }
        return localValue(xset, pi);
    }

    // per LouvainMethod
    void LMStart(const UGraph& grafo,
                 const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) {
        __lm_clusters.clear();
        __lm_node.clear();
        for (auto it_big_node_1 = clusters.begin(); it_big_node_1 != clusters.end(); ++it_big_node_1) {
            auto cm_c1_iter = __lm_node.insert({it_big_node_1->first, std::map<std::uint_fast64_t, double>()});

            auto c1_pair = localValue(*(it_big_node_1->second), nullptr);
            auto c1_value = c1_pair.first + c1_pair.second;

            cm_c1_iter.first->second.insert({it_big_node_1->first, c1_value});
            __lm_clusters[it_big_node_1->first] = c1_value;

            for (auto it_big_node_2 = std::next(it_big_node_1); it_big_node_2 != clusters.end(); ++it_big_node_2) {

                auto c2_pair = localValue(*(it_big_node_2->second), nullptr);
                auto c2_value = c2_pair.first + c2_pair.second;

                std::set<std::uint_fast64_t> c3_set(it_big_node_1->second->begin(), it_big_node_1->second->end());
                c3_set.insert(it_big_node_2->second->begin(), it_big_node_2->second->end());

                auto c3_pair = localValue(c3_set, nullptr);
                auto c3_value = c3_pair.first + c3_pair.second;

                cm_c1_iter.first->second.insert({it_big_node_2->first, c3_value - c1_value - c2_value});
            }
        }
    }

    double LMAt(std::uint_fast64_t node,
                std::uint_fast64_t cluster_from,
                std::uint_fast64_t cluster_to,
                std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) const {

        //double at_end = it_modularity_from->second + it_modularity_to->second;
        double at_end = std::accumulate(std::begin(__lm_clusters),
                                        std::end(__lm_clusters),
                                        0.0,
                                        [] (double value, const LMCLUSTER::value_type& p){
                                            return value + p.second;
                                        });


        {
            auto it_clusters_from = clusters.find(cluster_from);
            for (auto it_node_in_cluster = it_clusters_from->second->begin() ; it_node_in_cluster != it_clusters_from->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                at_end -= valore;
            }
        }
        {
            auto it_clusters_to = clusters.find(cluster_to);
            at_end += __lm_node.at(node).at(node);
            for (auto it_node_in_cluster = it_clusters_to->second->begin() ; it_node_in_cluster != it_clusters_to->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                at_end += valore;
            }
        }
        return at_end;
    };

    void LMUpdate(std::uint_fast64_t node,
                  std::uint_fast64_t cluster_from,
                  std::uint_fast64_t cluster_to,
                  std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) {

        // da utilizzare prima di modificare clusters
        {
            auto it_modularity_from = __lm_clusters.find(cluster_from);
            auto it_clusters_from = clusters.find(cluster_from);
            for (auto it_node_in_cluster = it_clusters_from->second->begin() ; it_node_in_cluster != it_clusters_from->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                it_modularity_from->second -= valore;
            }
        }
        {
            auto it_modularity_to = __lm_clusters.find(cluster_to);
            auto it_clusters_to = clusters.find(cluster_to);
            it_modularity_to->second += __lm_node.at(node).at(node);
            for (auto it_node_in_cluster = it_clusters_to->second->begin() ; it_node_in_cluster != it_clusters_to->second->end(); ++it_node_in_cluster) {
                auto node_id = *it_node_in_cluster;
                auto valore = (node <= node_id ? __lm_node.at(node).at(node_id) : __lm_node.at(node_id).at(node));
                it_modularity_to->second += valore;
            }
        }
    };

    double LMValue() const {
        double result = std::accumulate(std::begin(__lm_clusters),
                                        std::end(__lm_clusters),
                                        0.0,
                                        [] (double value, const LMCLUSTER::value_type& p){
                                            return value + p.second;
                                        });
        return result;
    };
private:
    void BuildFrom(std::map<std::uint_fast64_t, std::set<std::uint_fast64_t>>& communities) {
        std::fill(__in_arc->begin(), __in_arc->end(), 0);
        std::fill(__out_arc->begin(), __out_arc->end(), 0);
        std::fill(__prs->begin(), __prs->end(), 0.0);

        for (std::uint_fast64_t row = 0; row  < __grafo.Size(); ++row) {
            if (communities.find(row) != communities.end()) {
                auto& s = communities.at(row);
                for (auto v1 : s) {
                    for (auto v2 : s) {
                        __in_arc->at(row) += __grafo.Weight(v1, v2).first;
                    }
                    for (std::uint_fast64_t v2 = 0; v2  < __grafo.Size(); ++v2) {
                        if (s.find(v2) == s.end()) {
                            __out_arc->at(row) += __grafo.Weight(v1, v2).first;
                        }
                    }
                }
                __in_arc->at(row) /= 2;
                __prs->at(row) = ((double) __in_arc->at(row)) / ((double) (__in_arc->at(row) + __out_arc->at(row)));
            }
        }
        return;
    }
};


#endif /* CommunityMeasure_hpp */
