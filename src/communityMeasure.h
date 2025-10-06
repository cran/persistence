#ifndef CommunityMeasure_hpp
#define CommunityMeasure_hpp

#include <cstdint>
#include <memory>
#include <vector>
#include <numeric>
#include <limits>
#include <set>
#include <algorithm>
#include <fstream>
#include <filesystem>
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
    virtual void LMInit(std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>) = 0;
    virtual void LMStart(std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters) = 0;
    virtual std::pair<double, double> LMAt(std::uint_fast64_t node, std::uint_fast64_t cluster_from, std::uint_fast64_t cluster_to,
                        std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) const = 0;
    virtual void LMUpdate(std::uint_fast64_t node,
                          std::uint_fast64_t cluster_from,
                          std::uint_fast64_t cluster_to,
                          std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) = 0;
    virtual double LMValue() const = 0;
    virtual void LMRestart(std::uint_fast64_t, std::uint_fast64_t) = 0;
    virtual void LMRestart() = 0;

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

    std::vector<std::uint_fast64_t> __lm_active_node;
    std::uint_fast64_t __lm_active_node_number;
    
    std::vector<double> __lm_clusters_value;
    std::vector<double> __lm_cluster_strength;
    std::vector<double> __lm_cluster_in_arc;
    std::vector<double> __lm_cluster_start_in_arc;
    std::vector<std::vector<double>> __lm_node_to_cluster_in_arc_start;
    std::vector<std::vector<double>> __lm_node_to_cluster_in_arc_running;
    std::vector<std::vector<double>> __lm_cluster_to_cluster_in_arc;
    std::vector<double> __lm_node_strength;
    double __lm_value;
    
    std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> __lm_start_clusters;
public:
    PersistenceMeasure(UGraph& grafo) : __grafo(grafo) {
        __in_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __out_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __prs = std::make_shared<std::vector<double>>(__grafo.Size());
        
        __lm_active_node_number = 0;
        __lm_active_node.assign(__grafo.Size(), 0);
    }
    
    std::vector<std::pair<double, double>> globalValue(std::vector<std::uint_fast64_t>& membership, std::shared_ptr<std::vector<double>> pi) const {
        if (membership.size() != __grafo.Size()) {
            std::string err_str = "error: globalValue";
            throw_line(err_str);
        }

        auto num_node = __grafo.Size();

        std::uint_fast64_t maximum = *std::max_element(membership.begin(), membership.end());
        std::vector<std::pair<double, double>> local_values(maximum, std::make_pair(0.0, 0.0));

        for (std::uint_fast64_t i = 1; i <= maximum; i++) {
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
        auto lvs = localValueSplit(x, pi);
        auto somma_interni_interni = std::get<0>(lvs);
        auto somma_strength = std::get<1>(lvs);
        auto somma_pi = std::get<2>(lvs);

        double risultato_persistence = (2.0 * somma_interni_interni) / ((double) somma_strength);
        double risultato = risultato_persistence;
        return std::make_pair((std::isnan(risultato) || std::isinf(risultato) ? 0 : risultato), - somma_pi);
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
    
    std::tuple<double, double, double> localValueSplit(std::set<std::uint_fast64_t>& x, std::shared_ptr<std::vector<double>> pi) const {
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
        return std::make_tuple(somma_interni_interni, somma_strength, - somma_pi);
    }
    
    // per LouvainMethod
    void LMInit(std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters) {
        __lm_cluster_to_cluster_in_arc.assign(__grafo.Size(), std::vector<double>(__grafo.Size(), 0));

        std::uint_fast64_t k = 0;
        for (auto it_big_node = clusters->begin(); it_big_node != clusters->end(); ++it_big_node, ++k) {
            __lm_cluster_to_cluster_in_arc.at(it_big_node->first).at(it_big_node->first) = 0;
        }
        
        for (auto it_big_node_1 = clusters->begin(); it_big_node_1 != clusters->end(); ++it_big_node_1) {
            for (auto it_big_node_2 = std::next(it_big_node_1); it_big_node_2 != clusters->end(); ++it_big_node_2) {
                auto ctc = ClusterToClusterInArc(*(it_big_node_1->second), *(it_big_node_2->second));
                __lm_cluster_to_cluster_in_arc.at(it_big_node_1->first).at(it_big_node_2->first) = ctc;
                __lm_cluster_to_cluster_in_arc.at(it_big_node_2->first).at(it_big_node_1->first) = ctc;
            }
        }
        return;
    }
    
    void LMStart(std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters) {
        
        __lm_start_clusters = clusters;

        __lm_value = 0.0;
        __lm_clusters_value.assign(__grafo.Size(), std::numeric_limits<double>::quiet_NaN());
        __lm_cluster_strength.assign(__grafo.Size(), std::numeric_limits<double>::quiet_NaN());
        __lm_cluster_in_arc.assign(__grafo.Size(), 0);
        __lm_cluster_start_in_arc.assign(__grafo.Size(), 0);
        __lm_node_to_cluster_in_arc_start.assign(__grafo.Size(), std::vector<double>(__grafo.Size(), 0));
        __lm_node_to_cluster_in_arc_running.assign(__grafo.Size(), std::vector<double>(__grafo.Size(), 0));
        
        __lm_node_strength.assign(__grafo.Size(), std::numeric_limits<double>::quiet_NaN());
        __lm_active_node_number = clusters->size();

        std::uint_fast64_t k = 0;
        for (auto it_big_node = clusters->begin(); it_big_node != clusters->end(); ++it_big_node, ++k) {
            __lm_active_node.at(k) = it_big_node->first;
            auto big_node_values = localValueSplit(*(it_big_node->second), nullptr);
            auto big_node_in_arc = std::get<0>(big_node_values);
            auto big_node_strength = std::get<1>(big_node_values);
            double big_node_persistence = (2.0 * big_node_in_arc) / ((double) big_node_strength);
            __lm_value += big_node_persistence;

            __lm_node_to_cluster_in_arc_start.at(it_big_node->first).at(it_big_node->first) = 0;
            __lm_node_to_cluster_in_arc_running.at(it_big_node->first).at(it_big_node->first) = 0;
            __lm_cluster_in_arc.at(it_big_node->first) = big_node_in_arc;
            __lm_cluster_start_in_arc.at(it_big_node->first) = big_node_in_arc;
            __lm_node_strength.at(it_big_node->first) = big_node_strength;
            __lm_cluster_strength.at(it_big_node->first) = big_node_strength;
            __lm_clusters_value.at(it_big_node->first) = big_node_persistence;
        }
        
        for (auto it_big_node_1 = clusters->begin(); it_big_node_1 != clusters->end(); ++it_big_node_1) {
            for (auto it_big_node_2 = std::next(it_big_node_1); it_big_node_2 != clusters->end(); ++it_big_node_2) {
                auto ctc = ClusterToClusterInArc(*(it_big_node_1->second), *(it_big_node_2->second));
                __lm_node_to_cluster_in_arc_start.at(it_big_node_1->first).at(it_big_node_2->first) = ctc;
                __lm_node_to_cluster_in_arc_running.at(it_big_node_1->first).at(it_big_node_2->first) = ctc;
                __lm_node_to_cluster_in_arc_start.at(it_big_node_2->first).at(it_big_node_1->first) = ctc;
                __lm_node_to_cluster_in_arc_running.at(it_big_node_2->first).at(it_big_node_1->first) = ctc;
            }
        }
        return;
    }

    std::pair<double, double> LMAt(std::uint_fast64_t node,
                std::uint_fast64_t cluster_from,
                std::uint_fast64_t cluster_to,
                std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) const {

        auto cluster_from_in = __lm_cluster_in_arc.at(cluster_from);
        auto cluster_from_strength = __lm_cluster_strength.at(cluster_from);
        auto cluster_to_in = __lm_cluster_in_arc.at(cluster_to);
        auto cluster_to_strength = __lm_cluster_strength.at(cluster_to);
        auto node_to_cluster_from = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_from);
        auto node_to_cluster_to = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_to);
        auto node_strength = __lm_node_strength.at(node);
        auto node_start_in = __lm_cluster_start_in_arc.at(node);

        std::uint_fast64_t new_in_from = cluster_from_in - node_to_cluster_from - node_start_in;
        double r_from = ((2.0 * new_in_from) / ((double) (cluster_from_strength - node_strength)));
        auto r_old_from = __lm_clusters_value.at(cluster_from);

        std::uint_fast64_t new_in_to = cluster_to_in + node_to_cluster_to + node_start_in;
        double r_to = ((2.0 * new_in_to) / ((double) (cluster_to_strength + node_strength)));
        auto r_old_to = __lm_clusters_value.at(cluster_to);

        double new_value = __lm_value - r_old_from - r_old_to + (std::isnan(r_from) || std::isinf(r_from) ? 0 : r_from) + (std::isnan(r_to) || std::isinf(r_to) ? 0 : r_to);
        
        double var = (std::isnan(r_from) || std::isinf(r_from) ? 0 : r_from) - r_old_from + (std::isnan(r_to) || std::isinf(r_to) ? 0 : r_to) - r_old_to;
        return {new_value, var};
    };

    void LMUpdate(std::uint_fast64_t node,
                  std::uint_fast64_t cluster_from,
                  std::uint_fast64_t cluster_to,
                  std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) {
      
        //print(__lm_node_to_cluster_in_arc_running);
        
        auto cluster_from_in = __lm_cluster_in_arc.at(cluster_from);
        auto cluster_from_strength = __lm_cluster_strength.at(cluster_from);
        auto cluster_to_in = __lm_cluster_in_arc.at(cluster_to);
        auto cluster_to_strength = __lm_cluster_strength.at(cluster_to);
        auto node_to_cluster_from = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_from);
        auto node_to_cluster_to = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_to);
        auto node_strength = __lm_node_strength.at(node);
        auto node_start_in = __lm_cluster_start_in_arc.at(node);

        std::uint_fast64_t new_in_from = cluster_from_in - node_to_cluster_from - node_start_in;
        double per_from = (2.0 * new_in_from) / (1.0 * cluster_from_strength - node_strength);
        double r_from = (std::isnan(per_from) || std::isinf(per_from)  ? 0 : per_from);
        
        std::uint_fast64_t new_in_to = cluster_to_in + node_to_cluster_to + node_start_in;
        double per_to = (2.0 * new_in_to) / (1.0 * cluster_to_strength + node_strength);
        double r_to = (std::isnan(per_to) || std::isinf(per_to)  ? 0 : per_to);
        double old__lm_value = __lm_value;
        
        __lm_value = old__lm_value - __lm_clusters_value.at(cluster_from) - __lm_clusters_value.at(cluster_to) + (std::isnan(r_from) || std::isinf(r_from)  ? 0 : r_from) + (std::isnan(r_to) || std::isinf(r_to) ? 0 : r_to);

        __lm_cluster_in_arc.at(cluster_from) -= (__lm_node_to_cluster_in_arc_running.at(node).at(cluster_from) + node_start_in);
        for (std::uint_fast64_t k = 0; k < __lm_active_node_number; ++k) {
            auto n = __lm_active_node.at(k);
            if (n != node) {
                __lm_node_to_cluster_in_arc_running.at(n).at(cluster_from) -= __lm_node_to_cluster_in_arc_start.at(n).at(node);
            }
        }
        
        __lm_cluster_in_arc.at(cluster_to) += (__lm_node_to_cluster_in_arc_running.at(node).at(cluster_to) + node_start_in);
        for (std::uint_fast64_t k = 0; k < __lm_active_node_number; ++k) {
            auto n = __lm_active_node.at(k);
            if (n != node) {
                __lm_node_to_cluster_in_arc_running.at(n).at(cluster_to) += __lm_node_to_cluster_in_arc_start.at(node).at(n);
            }
        }
        for (std::uint_fast64_t k = 0; k < __lm_active_node_number; ++k) {
            auto n = __lm_active_node.at(k);
            auto value = __lm_node_to_cluster_in_arc_running.at(node).at(n);
            __lm_cluster_to_cluster_in_arc.at(cluster_from).at(n) -= value;
            __lm_cluster_to_cluster_in_arc.at(n).at(cluster_from) -= value;
            __lm_cluster_to_cluster_in_arc.at(cluster_to).at(n) += value;
            __lm_cluster_to_cluster_in_arc.at(n).at(cluster_to) += value;
            
            __lm_cluster_to_cluster_in_arc.at(n).at(n) = 0;
        }
        
        __lm_cluster_strength.at(cluster_to) += __lm_node_strength.at(node);
        __lm_cluster_strength.at(cluster_from) -= __lm_node_strength.at(node);

        __lm_clusters_value.at(cluster_from) = r_from;
        __lm_clusters_value.at(cluster_to) = r_to;
        
        return;
    };
    
    double LMValue() const {
        double result = std::accumulate(std::begin(__lm_clusters_value),
                                        std::end(__lm_clusters_value),
                                        0.0,
                                        [] (double value, double p){
                                            return value + (std::isnan(p) || std::isinf(p) ? 0 : p);
                                        });
        return result;
    };
    
    void LMRestart(std::uint_fast64_t from, std::uint_fast64_t to) {
        
        return;
    }
    
    void LMRestart() {

    }


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
    
    double ClusterToClusterInArc(std::set<std::uint_fast64_t>& cluster1, std::set<std::uint_fast64_t>& cluster2) const {
        double peso_totale = 0.0;
        for (auto v1_it = cluster1.begin(); v1_it != cluster1.end(); v1_it++) {
            auto v1 = *v1_it;
            for (auto v2_it = cluster2.begin(); v2_it != cluster2.end(); v2_it++) {
                auto v2 = *v2_it;
                auto peso = __grafo.Weight(v1, v2).second;
                peso_totale += peso;
            }
        }
        return peso_totale;
    }
    
    double ClusterStrength(std::set<std::uint_fast64_t>& cluster1, std::shared_ptr<std::vector<double>> pi) const {
        double strength = 0.0;// __grafo.Strength(nodo);
        for (auto v1_it = cluster1.begin(); v1_it != cluster1.end(); v1_it++) {
            auto v1 = *v1_it;
            strength += __grafo.Strength(v1);
        }
        return strength;
    }
};


// ***********************************************
// ***********************************************
// ***********************************************

class PersistenceModularityMeasure : public CommunityMeasure { // ex: PersistenceModularityMeasureV2
    // P_{\cal C}
//protected:
public:
    UGraph& __grafo;

    std::shared_ptr<std::vector<std::uint_fast64_t>> __in_arc;
    std::shared_ptr<std::vector<std::uint_fast64_t>> __out_arc;
    std::shared_ptr<std::vector<double>> __prs;
    std::vector<double> __pi;


    std::vector<std::uint_fast64_t> __lm_active_node;
    std::uint_fast64_t __lm_active_node_number;
    
    std::vector<double> __lm_clusters_value;
    std::vector<double> __lm_cluster_strength;
    std::vector<double> __lm_cluster_in_arc;
    std::vector<double> __lm_cluster_start_in_arc;
    std::vector<std::vector<double>> __lm_node_to_cluster_in_arc_start;
    std::vector<std::vector<double>> __lm_node_to_cluster_in_arc_running;
    std::vector<std::vector<double>> __lm_cluster_to_cluster_in_arc;
    std::vector<double> __lm_node_strength;
    double __lm_value;
    
    std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> __lm_start_clusters;
    
    std::fstream& __out_file;
public:
    PersistenceModularityMeasure(UGraph& grafo, std::fstream& out_file) : __grafo(grafo), __out_file(out_file) {
        __in_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __out_arc = std::make_shared<std::vector<std::uint_fast64_t>>(__grafo.Size());
        __prs = std::make_shared<std::vector<double>>(__grafo.Size());
        
        __lm_active_node_number = 0;
        __lm_active_node.assign(__grafo.Size(), 0);
    }

    std::vector<std::pair<double, double>> globalValue(std::vector<std::uint_fast64_t>& membership, std::shared_ptr<std::vector<double>> pi) const {
        if (membership.size() != __grafo.Size()) {
            std::string err_str = "error: globalValue";
            throw_line(err_str);
        }

        auto num_node = __grafo.Size();

        std::uint_fast64_t maximum = *std::max_element(membership.begin(), membership.end());
        std::vector<std::pair<double, double>> local_values(maximum, std::make_pair(0.0, 0.0));

        for (std::uint_fast64_t i = 1; i <= maximum; i++) {
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
        auto lvs = localValueSplit(x, pi);
        auto somma_interni_interni = std::get<0>(lvs);
        auto somma_strength = std::get<1>(lvs);
        auto somma_pi = std::get<2>(lvs);

        double risultato_persistence = (2.0 * somma_interni_interni) / ((double) somma_strength);
        double risultato_null = ((double) somma_strength) /  (__grafo.Strength());
        double risultato = risultato_persistence - risultato_null;
        return std::make_pair((std::isnan(risultato) || std::isinf(risultato) ? 0 : risultato), - somma_pi);
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
    void LMInit(std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters) {
        
        __lm_cluster_to_cluster_in_arc.assign(__grafo.Size(), std::vector<double>(__grafo.Size(), 0));

        std::uint_fast64_t k = 0;
        for (auto it_big_node = clusters->begin(); it_big_node != clusters->end(); ++it_big_node, ++k) {
            __lm_cluster_to_cluster_in_arc.at(it_big_node->first).at(it_big_node->first) = 0;
        }
        
        for (auto it_big_node_1 = clusters->begin(); it_big_node_1 != clusters->end(); ++it_big_node_1) {
            for (auto it_big_node_2 = std::next(it_big_node_1); it_big_node_2 != clusters->end(); ++it_big_node_2) {
                auto ctc = ClusterToClusterInArc(*(it_big_node_1->second), *(it_big_node_2->second));
                
                __lm_cluster_to_cluster_in_arc.at(it_big_node_1->first).at(it_big_node_2->first) = ctc;
                __lm_cluster_to_cluster_in_arc.at(it_big_node_2->first).at(it_big_node_1->first) = ctc;
            }
        }
        return;
    }
    
    void LMStart(std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters) {
        
        __lm_start_clusters = clusters;

        __lm_value = 0.0;
        __lm_clusters_value.assign(__grafo.Size(), std::numeric_limits<double>::quiet_NaN());
        __lm_cluster_strength.assign(__grafo.Size(), std::numeric_limits<double>::quiet_NaN());
        __lm_cluster_in_arc.assign(__grafo.Size(), 0);
        __lm_cluster_start_in_arc.assign(__grafo.Size(), 0);
        __lm_node_to_cluster_in_arc_start.assign(__grafo.Size(), std::vector<double>(__grafo.Size(), 0));
        __lm_node_to_cluster_in_arc_running.assign(__grafo.Size(), std::vector<double>(__grafo.Size(), 0));
        
        __lm_node_strength.assign(__grafo.Size(), std::numeric_limits<double>::quiet_NaN());
        __lm_active_node_number = clusters->size();

        std::uint_fast64_t k = 0;
        for (auto it_big_node = clusters->begin(); it_big_node != clusters->end(); ++it_big_node, ++k) {
            __lm_active_node.at(k) = it_big_node->first;
            auto big_node_values = localValueSplit(*(it_big_node->second), nullptr);
            auto big_node_in_arc = std::get<0>(big_node_values);
            auto big_node_strength = std::get<1>(big_node_values);
            auto big_node_pi = std::get<2>(big_node_values);
            double big_node_persistence = (2.0 * big_node_in_arc) / ((double) big_node_strength);
            double big_node_null = ((double) big_node_strength) /  (__grafo.Strength());
            double big_node_null_adj_persistence = big_node_persistence - big_node_null + big_node_pi;
            __lm_value += big_node_null_adj_persistence;

            __lm_node_to_cluster_in_arc_start.at(it_big_node->first).at(it_big_node->first) = 0;
            __lm_node_to_cluster_in_arc_running.at(it_big_node->first).at(it_big_node->first) = 0;
            __lm_cluster_in_arc.at(it_big_node->first) = big_node_in_arc;
            __lm_cluster_start_in_arc.at(it_big_node->first) = big_node_in_arc;
            __lm_node_strength.at(it_big_node->first) = big_node_strength;
            __lm_cluster_strength.at(it_big_node->first) = big_node_strength;
            __lm_clusters_value.at(it_big_node->first) = big_node_null_adj_persistence;
        }
        
        for (auto it_big_node_1 = clusters->begin(); it_big_node_1 != clusters->end(); ++it_big_node_1) {
            for (auto it_big_node_2 = std::next(it_big_node_1); it_big_node_2 != clusters->end(); ++it_big_node_2) {
                auto ctc = ClusterToClusterInArc(*(it_big_node_1->second), *(it_big_node_2->second));
                __lm_node_to_cluster_in_arc_start.at(it_big_node_1->first).at(it_big_node_2->first) = ctc;
                __lm_node_to_cluster_in_arc_running.at(it_big_node_1->first).at(it_big_node_2->first) = ctc;
                __lm_node_to_cluster_in_arc_start.at(it_big_node_2->first).at(it_big_node_1->first) = ctc;
                __lm_node_to_cluster_in_arc_running.at(it_big_node_2->first).at(it_big_node_1->first) = ctc;
            }
        }
        
        if (__out_file.is_open()) {
            auto membership = to_membership(*clusters, nullptr);
            __out_file << "LMStart;" << membership.first << __lm_value << "\n";
        }
        return;
    }

    std::pair<double, double> LMAt(std::uint_fast64_t node,
                std::uint_fast64_t cluster_from,
                std::uint_fast64_t cluster_to,
                std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) const {

        auto cluster_from_in = __lm_cluster_in_arc.at(cluster_from);
        auto cluster_from_strength = __lm_cluster_strength.at(cluster_from);
        auto cluster_to_in = __lm_cluster_in_arc.at(cluster_to);
        auto cluster_to_strength = __lm_cluster_strength.at(cluster_to);
        auto node_to_cluster_from = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_from);
        auto node_to_cluster_to = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_to);
        auto node_strength = __lm_node_strength.at(node);
        auto node_start_in = __lm_cluster_start_in_arc.at(node);

        //double risultato = ((2.0 * somma_interni_interni) / ((double) somma_strength)) - (((double) somma_strength) /  ((double) __grafo.Strength()));
        std::uint_fast64_t new_in_from = cluster_from_in - node_to_cluster_from - node_start_in;
        double r_from = ((2.0 * new_in_from) / ((double) (cluster_from_strength - node_strength))) -
                        (((double) (cluster_from_strength - node_strength)) /  ((double) __grafo.Strength()));
        auto r_old_from = __lm_clusters_value.at(cluster_from);

        std::uint_fast64_t new_in_to = cluster_to_in + node_to_cluster_to + node_start_in;
        double r_to = ((2.0 * new_in_to) / ((double) (cluster_to_strength + node_strength))) -
                (((double) (cluster_to_strength + node_strength)) /  ((double) __grafo.Strength()));
        auto r_old_to = __lm_clusters_value.at(cluster_to);

        double new_value = __lm_value - r_old_from - r_old_to + (std::isnan(r_from) || std::isinf(r_from) ? 0 : r_from) + (std::isnan(r_to) || std::isinf(r_to) ? 0 : r_to);
        
        
        
        /*
        {
            double value_check = 0;
            for (auto c : clusters) {
                std::set<std::uint_fast64_t> cl = *c.second;
                if (c.first == cluster_from) cl.erase(node);
                if (c.first == cluster_to) cl.insert(node);
                auto cluster_set = ClusterSet(cl);
                auto l_value = localValue(cluster_set, nullptr).first;
                value_check += l_value;
            }
            if (std::abs(value_check - new_value) > 0.001) {
                auto a = 3;
            }
        }*/
        /*{
            auto cluster_set_to = ClusterSet(*clusters.at(cluster_to));
            auto l_value_to_before = localValue(cluster_set_to, nullptr).first;
            cluster_set_to.insert(__lm_start_clusters->at(node)->begin(), __lm_start_clusters->at(node)->end());
            auto l_value_to = localValue(cluster_set_to, nullptr).first;
            
            auto cluster_set_from = ClusterSet(*clusters.at(cluster_from));
            auto l_value_from_before = localValue(cluster_set_from, nullptr).first;
            for (auto v : *__lm_start_clusters->at(node)) {
                cluster_set_from.erase(v);
            }
            auto l_value_from = localValue(cluster_set_from, nullptr).first;
            double misura_before_update = LMValue();
            double misura_new = misura_before_update - l_value_to_before - l_value_from_before + l_value_to + l_value_from;
            double misura_delta = new_value - misura_before_update;
            double l_value_delta = (l_value_to - l_value_to_before) + (l_value_from - l_value_from_before);
            if (std::abs(misura_delta - l_value_delta) > 1e-6) {
                auto aaaa = 3;
            }
        }*/
        
        double var = (std::isnan(r_from) || std::isinf(r_from) ? 0 : r_from) - r_old_from + (std::isnan(r_to) || std::isinf(r_to) ? 0 : r_to) - r_old_to;
        return {new_value, var};
    };

    void LMUpdate(std::uint_fast64_t node,
                  std::uint_fast64_t cluster_from,
                  std::uint_fast64_t cluster_to,
                  std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters) {
      
        //print(__lm_node_to_cluster_in_arc_running);
        
        auto cluster_from_in = __lm_cluster_in_arc.at(cluster_from);
        auto cluster_from_strength = __lm_cluster_strength.at(cluster_from);
        auto cluster_to_in = __lm_cluster_in_arc.at(cluster_to);
        auto cluster_to_strength = __lm_cluster_strength.at(cluster_to);
        auto node_to_cluster_from = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_from);
        auto node_to_cluster_to = __lm_node_to_cluster_in_arc_running.at(node).at(cluster_to);
        auto node_strength = __lm_node_strength.at(node);
        auto node_start_in = __lm_cluster_start_in_arc.at(node);

        std::uint_fast64_t new_in_from = cluster_from_in - node_to_cluster_from - node_start_in;
        double per_from = (2.0 * new_in_from) / (1.0 * cluster_from_strength - node_strength);
        double null_from = (cluster_from_strength - node_strength) / (1.0 * __grafo.Strength());
        double r_from = (std::isnan(per_from) || std::isinf(per_from)  ? 0 : per_from) - (std::isnan(null_from) || std::isinf(null_from)  ? 0 : null_from);
        
        std::uint_fast64_t new_in_to = cluster_to_in + node_to_cluster_to + node_start_in;
        double per_to = (2.0 * new_in_to) / (1.0 * cluster_to_strength + node_strength);
        double null_to = (cluster_to_strength + node_strength) / (1.0 * __grafo.Strength());
        double r_to = (std::isnan(per_to) || std::isinf(per_to)  ? 0 : per_to) - (std::isnan(null_to) || std::isinf(null_to)  ? 0 : null_to);
        double old__lm_value = __lm_value;
        
        __lm_value = old__lm_value - __lm_clusters_value.at(cluster_from) - __lm_clusters_value.at(cluster_to) + (std::isnan(r_from) || std::isinf(r_from)  ? 0 : r_from) + (std::isnan(r_to) || std::isinf(r_to) ? 0 : r_to);

        __lm_cluster_in_arc.at(cluster_from) -= (__lm_node_to_cluster_in_arc_running.at(node).at(cluster_from) + node_start_in);
        //for (std::uint_fast64_t n = 0; n < __lm_node_to_cluster_in_arc_start.size(); ++n) {
        for (std::uint_fast64_t k = 0; k < __lm_active_node_number; ++k) {
            auto n = __lm_active_node.at(k);
            if (n != node) {
                __lm_node_to_cluster_in_arc_running.at(n).at(cluster_from) -= __lm_node_to_cluster_in_arc_start.at(n).at(node);
            }
        }
        
//        print(__lm_node_to_cluster_in_arc_running);
        __lm_cluster_in_arc.at(cluster_to) += (__lm_node_to_cluster_in_arc_running.at(node).at(cluster_to) + node_start_in);
        //for (std::uint_fast64_t n = 0; n < __lm_node_to_cluster_in_arc_start.size(); ++n) {
        for (std::uint_fast64_t k = 0; k < __lm_active_node_number; ++k) {
            auto n = __lm_active_node.at(k);
            if (n != node) {
                __lm_node_to_cluster_in_arc_running.at(n).at(cluster_to) += __lm_node_to_cluster_in_arc_start.at(node).at(n);
            }
        }
        for (std::uint_fast64_t k = 0; k < __lm_active_node_number; ++k) {
            auto n = __lm_active_node.at(k);
            auto value = __lm_node_to_cluster_in_arc_running.at(node).at(n);
            __lm_cluster_to_cluster_in_arc.at(cluster_from).at(n) -= value;
            __lm_cluster_to_cluster_in_arc.at(n).at(cluster_from) -= value;
            __lm_cluster_to_cluster_in_arc.at(cluster_to).at(n) += value;
            __lm_cluster_to_cluster_in_arc.at(n).at(cluster_to) += value;
            
            __lm_cluster_to_cluster_in_arc.at(n).at(n) = 0;
        }
        
        __lm_cluster_strength.at(cluster_to) += __lm_node_strength.at(node);
        __lm_cluster_strength.at(cluster_from) -= __lm_node_strength.at(node);

        __lm_clusters_value.at(cluster_from) = r_from;
        __lm_clusters_value.at(cluster_to) = r_to;
        
       /* {
            for (auto& c1 : clusters) {
                auto ccheck_1 = ClusterSet(*c1.second);
                for (auto& c2 : clusters) {
                    if (c1.first == c2.first) continue;
                    auto ccheck_2 = ClusterSet(*c2.second);

                    auto v1 = ClusterToClusterInArc(ccheck_1, ccheck_2);
                    auto v2 = __lm_cluster_to_cluster_in_arc.at(c1.first).at(c2.first);
                    auto v3 = __lm_cluster_to_cluster_in_arc.at(c2.first).at(c1.first);
                    if (v1 != v2 || v1 != v3) {
                        auto a = 3;
                    }
                }
            }
            
        }*/
            
        //print(__lm_node_to_cluster_in_arc_running);
        if (__out_file.is_open()) {
            auto membership = to_membership(clusters, __lm_start_clusters);
            __out_file << "LMUpdate;" << membership.first << __lm_value << "\n";
        }
        
        return;
    };
    
    double LMValue() const {
        double result = std::accumulate(std::begin(__lm_clusters_value),
                                        std::end(__lm_clusters_value),
                                        0.0,
                                        [] (double value, double p){
                                            return value + (std::isnan(p) || std::isinf(p) ? 0 : p);
                                        });
        return result;
    };
    
    void LMRestart(std::uint_fast64_t from, std::uint_fast64_t to) {
        return;
    }
    
    void LMRestart() {
    }
    
    
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

    double ClusterToClusterInArc(std::set<std::uint_fast64_t>& cluster1, std::set<std::uint_fast64_t>& cluster2) const {
        double peso_totale = 0.0;
        for (auto v1_it = cluster1.begin(); v1_it != cluster1.end(); v1_it++) {
            auto v1 = *v1_it;
            for (auto v2_it = cluster2.begin(); v2_it != cluster2.end(); v2_it++) {
                auto v2 = *v2_it;
                auto peso = __grafo.Weight(v1, v2).second;
                peso_totale += peso;
            }
        }
        return peso_totale;
    }
    
    double ClusterStrength(std::set<std::uint_fast64_t>& cluster1, std::shared_ptr<std::vector<double>> pi) const {
        double strength = 0.0;// __grafo.Strength(nodo);
        for (auto v1_it = cluster1.begin(); v1_it != cluster1.end(); v1_it++) {
            auto v1 = *v1_it;
            strength += __grafo.Strength(v1);
        }
        return strength;
    }
    
    
    std::set<std::uint_fast64_t> ClusterSet(std::set<std::uint_fast64_t>& cluster) const {
        std::set<std::uint_fast64_t> r;// __grafo.Strength(nodo);
        for (auto v1 : cluster) {
            r.insert(__lm_start_clusters->at(v1)->begin(), __lm_start_clusters->at(v1)->end());
        }
        return r;
    }
    
    double ClusterPi(std::set<std::uint_fast64_t>& cluster1, std::shared_ptr<std::vector<double>> pi) const {
        double pi_value = 0.0; // (pi == nullptr ? 0 : pi->at(nodo));
        for (auto v1_it = cluster1.begin(); v1_it != cluster1.end(); v1_it++) {
            auto v1 = *v1_it;
            pi_value += (pi == nullptr ? 0 : pi->at(v1));
            
        }
        return -pi_value;
    }
    
    std::tuple<double, double, double> localValueSplit(std::set<std::uint_fast64_t>& x, std::shared_ptr<std::vector<double>> pi) const {
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
        //double risultato = ((2.0 * somma_interni_interni) / ((double) somma_strength)) - (((double) somma_strength) /  ((double) __grafo.Strength()));
        return std::make_tuple(somma_interni_interni, somma_strength, - somma_pi);
    }
        
    void print(std::vector<std::vector<double>>& m) {
        for (std::uint_fast64_t r = 0; r < m.size(); ++r) {
            for (std::uint_fast64_t c = 0; c < m.at(r).size(); ++c) {
                std::cout << (std::isnan(m.at(r).at(c)) ? "NNNNNNNN" : std::to_string(m.at(r).at(c)))  << " ";

            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    
    std::pair<std::string, std::vector<std::uint_fast64_t>> to_membership(const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& partizione,
                              std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> big_node) {
        std::string membership_string = "";
        std::vector<std::uint_fast64_t> membership_vector(__grafo.Size());
        std::uint_fast64_t m_index = 1;
        for (auto it = partizione.begin(); it != partizione.end(); ++it) {
            for (auto v : *it->second) {
                membership_vector.at(v) = m_index;
                if (big_node != nullptr) {
                    auto v_set_in_big_node = big_node->at(v);
                    for (auto k : *v_set_in_big_node) {
                        membership_vector.at(k) = m_index;
                    }
                }
            }
            ++m_index;
        }
        
        for (auto v : membership_vector) {
            membership_string += std::to_string(v) + ";";
        }
        return {membership_string, membership_vector};
    }
};




#endif /* CommunityMeasure_hpp */
