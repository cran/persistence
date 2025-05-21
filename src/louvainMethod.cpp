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

using CV = std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, std::uint_fast64_t>>;

const double TOLLERANZA = 1e-6;

//std::uint_fast64_t conta_diversi = 0;

// ***********************************************
// ***********************************************
// ***********************************************

std::shared_ptr<std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>> RenameEdges(std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>& old_edges_list, std::vector<std::uint_fast64_t>& vertex_rename) {
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
                                                                                                       std::map<std::uint_fast64_t, std::uint_fast64_t>& vertex_rename) {
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
                         std::map<std::uint_fast64_t, std::set<std::uint_fast64_t>>& edge_between_big_node) {
    edge_between_big_node.clear();
    for (auto it_c1 = clusters.begin(); it_c1 != clusters.end(); ++it_c1) {
        auto vc_v1_iter = edge_between_big_node.insert({it_c1->first, std::set<std::uint_fast64_t>()});
        for (auto it_c2 = std::next(it_c1); it_c2 != clusters.end(); ++it_c2) {
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
                vc_v1_iter.first->second.insert(it_c2->first);
            }
        }
    }
}



// ***********************************************
// ***********************************************
// ***********************************************

void EdgeVertexClusterStart(const std::map<std::uint_fast64_t, std::set<std::uint_fast64_t>>& edge_between_big_node,
                            const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters,
                            CV& edge_from_vertex_to_clusters,
                            CV& edge_to_clusters_from_vertex) {
    edge_from_vertex_to_clusters.clear();
    edge_to_clusters_from_vertex.clear();
    for (auto it_c1 = clusters.begin(); it_c1 != clusters.end(); ++it_c1) {
        auto vc_v1_iter = edge_from_vertex_to_clusters.insert({it_c1->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
        auto cv_c1_iter = edge_to_clusters_from_vertex.insert({it_c1->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
        for (auto it_c2 = std::next(it_c1); it_c2 != clusters.end(); ++it_c2) {
            auto vc_v2_iter = edge_from_vertex_to_clusters.insert({it_c2->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
            auto cv_c2_iter = edge_to_clusters_from_vertex.insert({it_c2->first, std::map<std::uint_fast64_t, std::uint_fast64_t>()});
            bool arco_trovato = edge_between_big_node.at(it_c1->first).find(it_c2->first) != edge_between_big_node.at(it_c1->first).end();
            if (arco_trovato) {
                vc_v1_iter.first->second.insert({it_c2->first, 1});
                cv_c1_iter.first->second.insert({it_c2->first, 1});
                vc_v2_iter.first->second.insert({it_c1->first, 1});
                cv_c2_iter.first->second.insert({it_c1->first, 1});
            }
        }
    }
}

// ***********************************************
// ***********************************************
// ***********************************************

void EdgeVertexClusterUpdate(std::uint_fast64_t v,
                             std::uint_fast64_t cluster_from,
                             std::uint_fast64_t cluster_to,
                             std::map<std::uint_fast64_t, std::set<std::uint_fast64_t>>& edge_between_clusters,
                             std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& clusters,
                             CV& edge_from_vertex_to_clusters,
                             CV& edge_to_clusters_from_vertex) {

    {
        auto it_cv_cluster_from = edge_to_clusters_from_vertex.find(cluster_from);
        //auto it_clusters_cluster_from = clusters.find(cluster_from);
        for (auto it_cv_cluster_from_vert = it_cv_cluster_from->second.begin(); it_cv_cluster_from_vert != it_cv_cluster_from->second.end(); ) {
            auto v1 = it_cv_cluster_from_vert->first;

            if (
                (v < v1 && edge_between_clusters.at(v).find(v1) != edge_between_clusters.at(v).end())
                ||
                (v1 < v && edge_between_clusters.at(v1).find(v) != edge_between_clusters.at(v1).end())
                ) {
            //if (grafo.Weight(v, v1).first) {     // controllo
                auto it_vc_v1 = edge_from_vertex_to_clusters.find(v1);
                auto it_vc_v1_cluster_from = it_vc_v1->second.find(cluster_from);
                --it_vc_v1_cluster_from->second;
                if (it_vc_v1_cluster_from->second == 0) {
                    it_vc_v1->second.erase(it_vc_v1_cluster_from);
                }
                --it_cv_cluster_from_vert->second;
                if (it_cv_cluster_from_vert->second == 0) {
                    it_cv_cluster_from->second.erase(it_cv_cluster_from_vert++);
                }
                /*if (it_clusters_cluster_from->second.find(v1) != it_clusters_cluster_from->second.end()) {
                    auto it_cv_cluster_from_v = it_cv_cluster_from->second.find(v);
                    --it_cv_cluster_from_v->second;
                    if (it_cv_cluster_from_v->second == 0) {
                        it_cv_cluster_from->second.erase(it_cv_cluster_from_v);
                    }
                }*/
            } else {
                ++it_cv_cluster_from_vert;
            }
        }
    }
    {
        auto it_cv_cluster_to = edge_to_clusters_from_vertex.find(cluster_to);
        //for (auto v1_iter : grafo.EdgesOf(v)) {
        for (auto v1_iter : edge_from_vertex_to_clusters) {
            auto v1 = v1_iter.first;
            if (!
                (
                 (v < v1 && edge_between_clusters.at(v).find(v1) != edge_between_clusters.at(v).end())
                 ||
                 (v1 < v && edge_between_clusters.at(v1).find(v) != edge_between_clusters.at(v1).end())
                 )
                ) continue;
            //if (!grafo.Weight(v, v1).first) continue;     // controllo

            auto it_cv_cluster_to_v1 = it_cv_cluster_to->second.insert({v1, 1});
            if (it_cv_cluster_to_v1.second == false) {
                ++it_cv_cluster_to_v1.first->second;
            }
            auto it_vc_v1 = edge_from_vertex_to_clusters.find(v1);
            auto it_vc_v1_cluster_to = it_vc_v1->second.insert({cluster_to, 1});
            if (it_vc_v1_cluster_to.second == false) {
                ++it_vc_v1_cluster_to.first->second;
            }
        }
        ;
    }
}

// ***********************************************
// ***********************************************
// ***********************************************

void louvainMethod(UGraph& grafo,
                   CommunityMeasure& misura,
                   std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters_start,
                   std::shared_ptr<std::vector<double>> pi_value,
                   std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>& clusters_final,
                   std::pair<double, std::vector<double>>& fo_best) {

    static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required");

    auto moveNodes = [] (
                         UGraph& grafo,
                         CommunityMeasure& misura,
                         const std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& big_node)
    {
        CV edge_from_vertex_to_clusters;
        CV edge_to_clusters_from_vertex;
        std::map<std::uint_fast64_t, std::set<std::uint_fast64_t>> edge_between_big_node;
        std::map<std::uint_fast64_t, std::uint_fast64_t> where_i_am;

        auto clusters_local = std::make_shared<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>();
        for(auto& it : big_node) {
            (*clusters_local)[it.first] = std::make_shared<std::set<std::uint_fast64_t>>();
            (*clusters_local)[it.first]->insert(it.first);
            where_i_am[it.first] = it.first;
        }

        //auto start_time_move_cpu = std::clock();
        //auto start_time_move_wall = std::chrono::high_resolution_clock::now();

        misura.LMStart(grafo, big_node);

        //auto end_time_move_cpu = std::clock();
        //auto end_time_move_wall = std::chrono::high_resolution_clock::now();
        //std::uint64_t milliseconds_move_cpu = std::max((std::uint64_t) (1000.0 * (((double) end_time_move_cpu) - ((double) start_time_move_cpu)) / CLOCKS_PER_SEC), (std::uint64_t) 0);
        //std::uint64_t milliseconds_move_wall = std::chrono::duration<double, std::milli>(end_time_move_wall - start_time_move_wall).count();
        //std::cout << "\t\tTempo moveNodes " << big_node.size() << ": " << (milliseconds_move_cpu / 1000.0) << "s (cpu), " << (milliseconds_move_wall / 1000.0) << "s (wall), " << "\n";

        BuildEdgeBetweenBigNode(grafo, big_node, edge_between_big_node);
        EdgeVertexClusterStart(edge_between_big_node, big_node, edge_from_vertex_to_clusters, edge_to_clusters_from_vertex);


        bool changed_global = false;
        while (true) {
            bool trovato_miglioramento = false;
            auto misura_best = misura.LMValue();

            // Cerco la migliore modifica che migliora clusters_result spostando un nodo da un cluster ad un altro
            for (auto& it_vertex : edge_from_vertex_to_clusters) {
                auto nodo = it_vertex.first;
                auto nodo_where = where_i_am[it_vertex.first];
                auto cluster_best = clusters_local->end();
                for (auto id_cluster : it_vertex.second) {
                    auto it_current = clusters_local->find(id_cluster.first);
                    if (nodo_where != it_current->first) {
                        auto misura_running = misura.LMAt(nodo, nodo_where, it_current->first, *clusters_local);
                      if (misura_running - misura_best > TOLLERANZA) {
                            misura_best = misura_running;
                            trovato_miglioramento = true;
                            cluster_best = it_current;
                        }
                    }
                }
                if (cluster_best != clusters_local->end()) {
                    // trovato il cluster eseguo la modifica
                    EdgeVertexClusterUpdate(nodo, nodo_where, cluster_best->first, edge_between_big_node, *clusters_local, edge_from_vertex_to_clusters, edge_to_clusters_from_vertex);
                    misura.LMUpdate(nodo, nodo_where, cluster_best->first, *clusters_local);
                    auto t = clusters_local->find(nodo_where);
                    t->second->erase(nodo);
                    clusters_local->at(nodo_where)->erase(nodo);
                    cluster_best->second->insert(nodo);
                    where_i_am.at(nodo) = cluster_best->first;
                    changed_global = true;

                    if (t->second->size() == 0) {
                        edge_to_clusters_from_vertex.erase(t->first);
                        clusters_local->erase(t);
                    }
                    ;
                }
            }
            if (!trovato_miglioramento) {
                // nessuno spostamento di nodi migliora clusters_result
                break;
            }
        }
        std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> clusters_final = nullptr;
        if (changed_global) {
            clusters_final = std::make_shared<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>();
            for (auto it = clusters_local->begin(); it != clusters_local->end(); ++it) {
                if (it->second->size() > 0) {
                    auto minimo = *std::min_element(it->second->begin(), it->second->end());
                    auto valori = std::make_shared<std::set<std::uint_fast64_t>>();
                    for (auto v : *(it->second)) {
                        valori->insert(big_node.at(v)->begin(), big_node.at(v)->end());
                    }
                    (*clusters_final)[minimo] = valori;
                }
            }
        }



        return std::make_tuple(changed_global, clusters_final);
    };

    clusters_final = clusters_start;

    while (true) {

        auto result_move = moveNodes(grafo, misura, *clusters_final);

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

