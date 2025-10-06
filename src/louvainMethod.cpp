#include <set>
#include <cstdint>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <limits>
#include <algorithm>

#include "eccezioni.h"
#include "random.h"
#include "communityMeasure.h"
#include "grafi.h"

// ***********************************************
// ***********************************************
// ***********************************************

const double TOLLERANZA = 1e-6;

// ***********************************************
// ***********************************************
// ***********************************************

std::shared_ptr<std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>>
RenameEdges(std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>& old_edges_list,
            std::vector<std::uint_fast64_t>& vertex_rename) {
    auto new_edges_list = std::make_shared<std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>>();
    for (auto& edge : old_edges_list) {
        auto v1 = std::get<0>(edge);
        auto v2 = std::get<1>(edge);
        auto weight = std::get<2>(edge);
        new_edges_list->push_back({vertex_rename.at(v1), vertex_rename.at(v2), weight});
    }
    return new_edges_list;
}

// ***********************************************
// ***********************************************
// ***********************************************

std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> RenamePartition(
                                                                                                             std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& old_partition,
                                                                                                             //std::map<std::uint_fast64_t, std::uint_fast64_t>& vertex_rename
                                                                                                             std::vector<std::uint_fast64_t>& vertex_rename
                                                                                                             ) {
    auto new_partition = std::make_shared<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>();
    for (auto& p_iter : old_partition) {
        auto p_new = std::make_shared<std::set<std::uint_fast64_t>>();
        for (auto& v : *p_iter.second) {
            p_new->insert(vertex_rename.at(v));
        }
        new_partition->insert({vertex_rename.at(p_iter.first), p_new});
    }
    return new_partition;
}

// ***********************************************
// ***********************************************
// ***********************************************

void BuildMembership(std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters_start,
                           std::vector<std::uint_fast64_t>& membership) {
    std::uint_fast64_t m_index = 1;
    for (auto it = clusters_start.begin(); it != clusters_start.end(); ++it) {
        for (auto v : *it->second) {
            membership.at(v) = m_index;
        }
        ++m_index;
    }
};

// ***********************************************
// ***********************************************
// ***********************************************

void BuildEdgeBetweenBigNode(const UGraph& grafo,
                             const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters,
                             std::vector<std::vector<bool>>& edge_between_start_cluster) {
    edge_between_start_cluster.assign(grafo.Size(), std::vector<bool>(grafo.Size(), false));
    for (auto it_c1 = clusters.begin(); it_c1 != clusters.end(); ++it_c1) {
        auto c1 = it_c1->first;
        for (auto it_c2 = std::next(it_c1); it_c2 != clusters.end(); ++it_c2) {
            auto c2 = it_c2->first;
            bool arco_trovato = false;
            for (auto& v1 : *it_c1->second) {
                for (auto& v2 : *it_c2->second) {
                    if (grafo.Weight(v1, v2).first) {
                        arco_trovato = true;
                        break;
                    }
                    if (arco_trovato) break;
                }
            }
            if (arco_trovato) {
                edge_between_start_cluster.at(c1).at(c2) = true;
                edge_between_start_cluster.at(c2).at(c1) = true;
            }
        }
    }
}

// ***********************************************
// ***********************************************
// ***********************************************

void EdgeVertexClusterStart(const std::vector<std::vector<bool>>& edge_between_start_cluster,
                            const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters,
                            std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>>& edge_vertex_clusters,
                            std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>>& edge_clusters_vertex) {
    edge_vertex_clusters.clear();
    edge_clusters_vertex.clear();
    for (auto it_c1 = clusters.begin(); it_c1 != clusters.end(); ++it_c1) {
        auto vc_v1_iter = edge_vertex_clusters.insert({it_c1->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
        auto cv_c1_iter = edge_clusters_vertex.insert({it_c1->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
        for (auto it_c2 = std::next(it_c1); it_c2 != clusters.end(); ++it_c2) {
            auto vc_v2_iter = edge_vertex_clusters.insert({it_c2->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
            auto cv_c2_iter = edge_clusters_vertex.insert({it_c2->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
            bool arco_trovato = edge_between_start_cluster.at(it_c1->first).at(it_c2->first);
            if (arco_trovato) {
                vc_v1_iter.first->second.insert({it_c2->first, 1});
                cv_c1_iter.first->second.insert({it_c2->first, 1});
                vc_v2_iter.first->second.insert({it_c1->first, 1});
                cv_c2_iter.first->second.insert({it_c1->first, 1});
            }
        }
    }
    return;
}

// ***********************************************
// ***********************************************
// ***********************************************

void EdgeVertexClusterUpdate(std::uint_fast64_t v,
                             std::uint_fast64_t cluster_from,
                             std::uint_fast64_t cluster_to,
                             std::vector<std::vector<bool>>& edge_between_start_cluster,
                             std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters,
                             std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>>& edge_vertex_clusters,
                             std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>>& edge_clusters_vertex)

{
    {
        auto it_cv_cluster_from = edge_clusters_vertex.find(cluster_from);
        for (auto it_cv_cluster_from_vert = it_cv_cluster_from->second.begin(); it_cv_cluster_from_vert != it_cv_cluster_from->second.end(); ) {
            auto v1 = it_cv_cluster_from_vert->first;

            if (edge_between_start_cluster.at(v).at(v1))
            {
                auto it_vc_v1 = edge_vertex_clusters.find(v1);
                auto it_vc_v1_cluster_from = it_vc_v1->second.find(cluster_from);
                --it_vc_v1_cluster_from->second;
                if (it_vc_v1_cluster_from->second == 0) {
                    it_vc_v1->second.erase(it_vc_v1_cluster_from);
                }

                --it_cv_cluster_from_vert->second;
                if (it_cv_cluster_from_vert->second == 0) {
                    it_cv_cluster_from->second.erase(it_cv_cluster_from_vert++);
                } else {
                    ++it_cv_cluster_from_vert;
                }
            } else {
                ++it_cv_cluster_from_vert;
            }
        }
    }
    {
        auto it_cv_cluster_to = edge_clusters_vertex.find(cluster_to);
        for (auto v1_iter : edge_vertex_clusters) {
            auto v1 = v1_iter.first;
            if (!edge_between_start_cluster.at(v).at(v1)) continue;

            auto it_cv_cluster_to_v1 = it_cv_cluster_to->second.insert({v1, 1});
            if (it_cv_cluster_to_v1.second == false) {
                ++it_cv_cluster_to_v1.first->second;
            }
            auto it_vc_v1 = edge_vertex_clusters.find(v1);
            auto it_vc_v1_cluster_to = it_vc_v1->second.insert({cluster_to, 1});
            if (it_vc_v1_cluster_to.second == false) {
                ++it_vc_v1_cluster_to.first->second;
            }
        }
    }
    return;
}

// ***********************************************
// ***********************************************
// ***********************************************

std::tuple<bool, std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>>
moveNodes(std::vector<std::uint_fast64_t>& where_i_am
            ,std::vector<std::vector<bool>>& edge_between_start_cluster
            ,std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>>& edge_vertex_clusters
            ,std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>>& edge_clusters_vertex
            ,Random& rnd
            ,UGraph& grafo
            ,CommunityMeasure& misura
            ,std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> big_node
            ,std::fstream& debug_file
          )
{
    //std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>> edge_vertex_clusters_OLD;
    //std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>> edge_clusters_vertex_OLD;
    std::vector<std::uint_fast64_t> vertex_order(big_node->size(), 0);

    auto clusters_local = std::make_shared<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>();
    {
        std::uint_fast64_t k = 0;
        for(auto& it : *big_node) {
            (*clusters_local)[it.first] = std::make_shared<std::set<std::uint_fast64_t>>();
            (*clusters_local)[it.first]->insert(it.first);
            where_i_am[it.first] = it.first;
            vertex_order.at(k) = it.first;
            ++k;
        }
    }

    BuildEdgeBetweenBigNode(grafo, *big_node, edge_between_start_cluster);

    misura.LMStart(big_node);

    bool changed_global = false;
    bool trovato_miglioramento = true;

    auto misura_best = misura.LMValue();

    while (trovato_miglioramento)
    {

        rnd.RndShuffle(vertex_order);
        trovato_miglioramento = false;
        // Cerco la migliore modifica che migliora clusters_result spostando fondendo i cluster iniziali
        for (std::uint_fast64_t k = 0; k < vertex_order.size(); ++k) {
            auto nodo = vertex_order.at(k);
            auto clusters_connected_to_nodo = edge_vertex_clusters.at(nodo);
            auto nodo_to_clusters = edge_vertex_clusters.at(nodo);
            auto nodo_where = where_i_am[nodo];
            std::list<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>::iterator> clusters_best;
            for (auto id_cluster : clusters_connected_to_nodo) {
                auto it_current = clusters_local->find(id_cluster.first);
                if (nodo_where != it_current->first) {
                    auto misura_at = misura.LMAt(nodo, nodo_where, it_current->first, *clusters_local);
                    auto misura_running = misura_at.first;

                    if (misura_running - misura_best >= TOLLERANZA) {
                        misura_best = misura_running;
                        trovato_miglioramento = true;
                        clusters_best.clear();
                        clusters_best.push_back(it_current);
                    } else if (0 < misura_running - misura_best && misura_running - misura_best < TOLLERANZA) {
                        clusters_best.push_back(it_current);
                    }
                }
            }
            if (clusters_best.size() > 0) {
                // trovato il cluster eseguo la modifica

                auto position = rnd.RndNextInt(0, clusters_best.size() - 1);
                auto cluster_best = clusters_best.begin();
                std::advance(cluster_best, position);
                EdgeVertexClusterUpdate(nodo, nodo_where, (*cluster_best)->first, edge_between_start_cluster, *clusters_local, edge_vertex_clusters, edge_clusters_vertex);
                auto t = clusters_local->find(nodo_where);
                t->second->erase(nodo);
                (*cluster_best)->second->insert(nodo);
                where_i_am.at(nodo) = (*cluster_best)->first;
                changed_global = true;

                if (t->second->size() == 0) {
                    edge_clusters_vertex.erase(nodo_where);
                    clusters_local->erase(t);
                }
                misura.LMUpdate(nodo, nodo_where, (*cluster_best)->first, *clusters_local);

            }
        }
    }
    std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters_final = nullptr;
    if (changed_global) {
        clusters_final = std::make_shared<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>();
        {
            auto it_edge_vertex_clusters = edge_vertex_clusters.begin();
            auto it_edge_clusters_vertex = edge_clusters_vertex.begin();
            for (auto it = clusters_local->begin(); it != clusters_local->end(); ++it) {

                if (it->second->size() == 0) {
                    std::string err_str = "Wrong size ";
                    throw_line(err_str);
                }
                while (it_edge_vertex_clusters->first != it->first) {
                    it_edge_vertex_clusters = edge_vertex_clusters.erase(it_edge_vertex_clusters);
                }
                while (it_edge_clusters_vertex->first != it->first) {
                    it_edge_clusters_vertex = edge_vertex_clusters.erase(it_edge_clusters_vertex);
                }
                // non vuoto

                auto valori = std::make_shared<std::set<std::uint_fast64_t>>();
                for (auto v : *(it->second)) {
                    valori->insert(big_node->at(v)->begin(), big_node->at(v)->end());
                }
                it_edge_vertex_clusters->second.clear();
                for (auto it_edge_clusters_vertex_in = it_edge_clusters_vertex->second.begin(); it_edge_clusters_vertex_in != it_edge_clusters_vertex->second.end(); ++it_edge_clusters_vertex_in) {
                    auto v_in = where_i_am[it_edge_clusters_vertex_in->first];
                    if (v_in != it_edge_clusters_vertex->first) {
                        it_edge_vertex_clusters->second.insert({v_in, 1});
                    }
                }
                if (valori->size() > 0) {
                    (*clusters_final)[it->first] = valori;
                }
                ++it_edge_vertex_clusters;
                ++it_edge_clusters_vertex;
            }
            while (it_edge_vertex_clusters != edge_vertex_clusters.end()) {
                it_edge_vertex_clusters = edge_vertex_clusters.erase(it_edge_vertex_clusters);
            }
            edge_clusters_vertex=edge_vertex_clusters;

        }
        misura.LMRestart();
    }

    return std::make_tuple(changed_global, clusters_final);
}

// ***********************************************
// ***********************************************
// ***********************************************

void louvainMethod(UGraph& grafo,
                   CommunityMeasure& misura,
                   std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters_start,
                   std::shared_ptr<std::vector<double>> pi_value,
                   std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>& clusters_final,
                   std::pair<double, std::vector<double>>& fo_best,
                   std::fstream& debug_file) {

    static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");

    std::uint_fast64_t seed_rnd = RandomUni::GENERATORE_SEED_RANDOM.RndNextInt(0, std::numeric_limits<std::uint_fast64_t>::max());
    RandomUni rnd(seed_rnd);
    std::vector<std::uint_fast64_t> where_i_am(grafo.Size(), 0);
    std::vector<std::vector<bool>> edge_between_start_cluster(grafo.Size(), std::vector<bool>(grafo.Size(), false));

    std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>> edge_vertex_clusters;
    std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>> edge_clusters_vertex;

    clusters_final = clusters_start;

    misura.LMInit(clusters_final);
    BuildEdgeBetweenBigNode(grafo, *clusters_final, edge_between_start_cluster);
    EdgeVertexClusterStart(edge_between_start_cluster, *clusters_final, edge_vertex_clusters, edge_clusters_vertex);

    while (true) {
        if (clusters_final->size() == 1) {
            break;
        }
        auto result_move = moveNodes(where_i_am
                                     ,edge_between_start_cluster
                                     ,edge_vertex_clusters
                                     ,edge_clusters_vertex
                                     ,rnd
                                     ,grafo
                                     ,misura
                                     ,clusters_final
                                     ,debug_file
                                     );
        if (!std::get<0>(result_move)) {
            break;
        }
        clusters_final = std::get<1>(result_move);
    }

    std::vector<std::uint_fast64_t> membership(grafo.Size());
    BuildMembership(*clusters_final, membership);
    auto final_values = misura.globalValue(membership, nullptr);
    fo_best.first = 0;
    fo_best.second.resize(final_values.size());
    for (std::uint_fast64_t v = 0; v < final_values.size(); ++v) {
        auto value = final_values.at(v).first + final_values.at(v).second;
        fo_best.first += value;
        fo_best.second.at(v) = value;
    }

    return;
}

