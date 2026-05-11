/**
 * @file milano_algorithm.h
 */

#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <concepts>
#include <thread>
#include <stop_token>
#include <atomic>
#include <mutex>
#include <functional>
#include "grafi.h"
#include "community_measure.h"
#include "random.h"

namespace Milano {
    
    struct MoveNodeResult {
        Partition partition;
        Partition inner_partition;
        bool improvement;
        std::vector<std::uint64_t> node2com;
        std::uint64_t while_iterations;
    };
    
    struct MilanoResult {
        Partition communities;
        std::vector<double> cluster_scores;
        std::vector<std::uint64_t> move_node_iterations;
    };
    
    struct RawEdge {
        std::uint64_t u;
        std::uint64_t v;
        double wt;
        
        bool operator<(const RawEdge& other) const {
            if (u != other.u) return u < other.u;
            return v < other.v;
        }
    };
    
    template <IsCommunityMeasure MeasureType>
    inline MoveNodeResult moveNode(
                                   const GraphContext& ctx,
                                   MeasureType& measure,
                                   Partition partition,
                                   Partition inner_partition, // <-- Ora viene passata dal chiamante
                                   const Partition& graph_nodes_mapping,
                                   double threshold,
                                   std::uint64_t max_level,
                                   Milano::Random& rng)
    {
        std::uint64_t num_nodes = ctx.graph->Size();
        
        // Inizializza node2com basandosi sulla inner_partition fornita
        std::vector<std::uint64_t> node2com(num_nodes, 0);
        for(size_t c = 0; c < inner_partition.size(); ++c) {
            for(auto u : inner_partition[c]) {
                if (u < num_nodes) node2com[u] = c;
            }
        }
        
        std::vector<std::uint64_t> rand_nodes(num_nodes);
        std::iota(rand_nodes.begin(), rand_nodes.end(), 0);
        rng.RndShuffle(rand_nodes);
        
        measure.start_level(ctx, inner_partition);
        
        bool improvement = false;
        int nb_moves = 1;
        std::uint64_t while_iterations = 0;
        
        std::vector<double> weights_array(num_nodes, 0.0);
        std::vector<std::uint64_t> touched_communities;
        touched_communities.reserve(1024);
        
        while (nb_moves > 0) {
            if (max_level > 0 && while_iterations >= max_level) break;
            while_iterations++;
            nb_moves = 0;
            
            for (std::uint64_t u : rand_nodes) {
                std::uint64_t old_com = node2com[u];
                std::uint64_t best_com = old_com;
                
                double best_partition_value = measure.current_measure();
                double best_gain_for_update = 0.0;
                double best_wt_new = 0.0;
                double degree = ctx.node_degrees[u];
                
                for (const auto& [nbr, wt] : ctx.graph->EdgesOf(u)) {
                    std::uint64_t nbr_com = node2com[nbr];
                    if (weights_array[nbr_com] == 0.0) {
                        touched_communities.push_back(nbr_com);
                    }
                    weights_array[nbr_com] += wt;
                }
                
                double wt_old = weights_array[old_com];
                
                for (std::uint64_t nbr_com : touched_communities) {
                    if (nbr_com == old_com) continue;
                    
                    double wt_new = weights_array[nbr_com];
                    double candidate_value = measure.evaluate_move(u, old_com, nbr_com, degree, wt_old, wt_new, ctx.self_loops[u]);
                    
                    if (candidate_value - best_partition_value > threshold) {
                        best_partition_value = candidate_value;
                        best_com = nbr_com;
                        best_gain_for_update = candidate_value - measure.current_measure();
                        best_wt_new = wt_new;
                    }
                }
                
                for (std::uint64_t com : touched_communities) {
                    weights_array[com] = 0.0;
                }
                touched_communities.clear();
                
                if (best_com != old_com) {
                    measure.update(u, old_com, best_com, degree, wt_old, best_wt_new, ctx.self_loops[u], best_gain_for_update);
                    
                    const auto& com_nodes = graph_nodes_mapping[u];
                    for (auto orig_u : com_nodes) {
                        std::erase(partition[old_com], orig_u);
                        partition[best_com].push_back(orig_u);
                    }
                    
                    std::erase(inner_partition[old_com], u);
                    inner_partition[best_com].push_back(u);
                    
                    improvement = true;
                    nb_moves += 1;
                    node2com[u] = best_com;
                }
            }
        }
        
        Partition filtered_partition;
        Partition filtered_inner;
        std::vector<std::uint64_t> new_node2com(num_nodes, 0);
        std::uint64_t new_idx = 0;
        
        for (size_t i = 0; i < inner_partition.size(); ++i) {
            if (!inner_partition[i].empty()) {
                filtered_partition.push_back(std::move(partition[i]));
                filtered_inner.push_back(std::move(inner_partition[i]));
                for (auto u : filtered_inner.back()) {
                    new_node2com[u] = new_idx;
                }
                new_idx++;
            }
        }
        
        return {std::move(filtered_partition), std::move(filtered_inner), improvement, std::move(new_node2com), while_iterations};
    }
    
    template <IsCommunityMeasure MeasureType>
    inline MilanoResult milano_communities(
                                           const Graph& G,
                                           MeasureType& measure,
                                           const Partition& initial_partition = {}, // <-- Nuovo parametro
                                           double threshold = 0.0000001,
                                           std::uint64_t max_level = 0,
                                           bool collect_stats = false,
                                           std::uint64_t seed = 0)
    {
        if (G.Size() == 0) return {{}, {}, {}};
        
        std::uint64_t num_nodes = G.Size();
        Partition partition;
        Partition inner_partition;
        Partition graph_nodes_mapping(num_nodes);
        
        // Il mapping iniziale (livello 0) è sempre 1 a 1 verso il grafo di base
        for(std::uint64_t u = 0; u < num_nodes; ++u) {
            graph_nodes_mapping[u].push_back(u);
        }
        
        // Gestione della partizione iniziale
        if (initial_partition.empty()) {
            partition.resize(num_nodes);
            inner_partition.resize(num_nodes);
            for(std::uint64_t u = 0; u < num_nodes; ++u) {
                partition[u].push_back(u);
                inner_partition[u].push_back(u);
            }
        } else {
            partition = initial_partition;
            // Pad di sicurezza per garantire che ci sia spazio nel caso un nodo voglia creare una comunità isolata
            partition.resize(num_nodes);
            inner_partition = initial_partition;
            inner_partition.resize(num_nodes);
        }
        
        if (G.TotalWeight() == 0.0) {
            std::vector<double> scores;
            if (collect_stats) scores.assign(G.Size(), 0.0);
            return {partition, scores, {0}};
        }
        
        measure.init_graph(G);

        Milano::Random rng(seed);
        
        std::vector<double> initial_degrees(num_nodes);
        for(std::uint64_t u = 0; u < num_nodes; ++u) {
            initial_degrees[u] = G.Strength(u);
        }
        std::vector<double> initial_self_loops(num_nodes, 0.0);
        
        GraphContext ctx(&G, std::move(initial_degrees), std::move(initial_self_loops));
        std::unique_ptr<Graph> generated_G = nullptr;
        std::vector<std::uint64_t> iterations_history;
        
        auto level_result = moveNode(ctx, measure, partition, inner_partition, graph_nodes_mapping, threshold, max_level, rng);
        partition = level_result.partition;
        
        if (collect_stats) iterations_history.push_back(level_result.while_iterations);
        bool improvement = level_result.improvement;
        
        while (improvement) {
            std::vector<double> new_node_degrees(level_result.inner_partition.size(), 0.0);
            std::vector<double> new_self_loops(level_result.inner_partition.size(), 0.0);
            
            std::vector<RawEdge> raw_edges;
            raw_edges.reserve(ctx.graph->Size() * 2);
            
            for (std::uint64_t u = 0; u < ctx.graph->Size(); ++u) {
                std::uint64_t com_u = level_result.node2com[u];
                new_node_degrees[com_u] += ctx.node_degrees[u];
                new_self_loops[com_u] += ctx.self_loops[u];
                
                for (const auto& [v, wt] : ctx.graph->EdgesOf(u)) {
                    std::uint64_t com_v = level_result.node2com[v];
                    if (com_u != com_v) {
                        if (com_u < com_v) raw_edges.push_back({com_u, com_v, wt});
                    } else {
                        if (u < v) new_self_loops[com_u] += wt;
                    }
                }
            }
            
            std::sort(raw_edges.begin(), raw_edges.end());
            
            EdgeList new_edges;
            if (!raw_edges.empty()) {
                std::uint64_t curr_u = raw_edges[0].u;
                std::uint64_t curr_v = raw_edges[0].v;
                double curr_wt = raw_edges[0].wt;
                
                for (size_t i = 1; i < raw_edges.size(); ++i) {
                    if (raw_edges[i].u == curr_u && raw_edges[i].v == curr_v) {
                        curr_wt += raw_edges[i].wt;
                    } else {
                        new_edges.push_back({curr_u, curr_v, curr_wt});
                        curr_u = raw_edges[i].u;
                        curr_v = raw_edges[i].v;
                        curr_wt = raw_edges[i].wt;
                    }
                }
                new_edges.push_back({curr_u, curr_v, curr_wt});
            }
            
            generated_G = ctx.graph->Make(level_result.inner_partition.size(), new_edges);
            ctx.graph = generated_G.get();
            ctx.node_degrees = std::move(new_node_degrees);
            ctx.self_loops = std::move(new_self_loops);
            graph_nodes_mapping = partition;
            
            // Per i livelli successivi, i nodi partono sempre isolati nel "Super-Grafo"
            Partition next_inner_partition(level_result.inner_partition.size());
            for(std::uint64_t u = 0; u < level_result.inner_partition.size(); ++u) {
                next_inner_partition[u].push_back(u);
            }
            
            level_result = moveNode(ctx, measure, partition, next_inner_partition, graph_nodes_mapping, threshold, max_level, rng);
            partition = level_result.partition;
            
            if (collect_stats) iterations_history.push_back(level_result.while_iterations);
            improvement = level_result.improvement;
        }
        
        std::vector<double> final_cluster_scores;
        if (collect_stats) {
            final_cluster_scores.reserve(partition.size());
            for (const auto& com : partition) {
                final_cluster_scores.push_back(measure.localScore(G, com));
            }
        }
        
        return {std::move(partition), std::move(final_cluster_scores), std::move(iterations_history)};
    }
    
    /**
     * @brief Esegue restart multipli dell'algoritmo Milano e restituisce la partizione migliore.
     * @details Esegue l'algoritmo per ogni seed fornito nel vettore e salva il risultato migliore in memoria.
     */
    template <IsCommunityMeasure MeasureType>
    inline MilanoResult milano_restart(
                                       const Graph& G,
                                       const MeasureType& measure_base,
                                       const Partition& initial_partition = {},
                                       double threshold = 0.0000001,
                                       std::uint64_t max_level = 0,
                                       bool collect_stats = false,
                                       const std::vector<std::uint64_t>& seeds = {0},
                                       std::uint64_t n_threads = 0)
    {
        if (seeds.empty()) return {{}, {}, {}};
        
        // Imposta il numero di thread in base all'hardware se n_threads è 0
        if (n_threads == 0) {
            n_threads = static_cast<std::uint64_t>(std::thread::hardware_concurrency());
            if (n_threads == 0) n_threads = 2; // Fallback di sicurezza
        }
        
        // Assicuriamoci di non lanciare più thread del numero totale di starts
        std::uint64_t total_runs = seeds.size();
        if (n_threads > total_runs) {
            n_threads = total_runs;
        }
        
        MilanoResult global_best_result;
        double global_best_score = -std::numeric_limits<double>::infinity();
        std::mutex best_mutex;
        
        // Contatore atomico per distribuire dinamicamente i seed ai thread
        std::atomic<std::size_t> task_index{0};
        
        // Funzione worker eseguita da ogni thread
        auto worker = [&]() {
            MilanoResult local_best_result;
            double local_best_score = -std::numeric_limits<double>::infinity();
            
            while (true) {
                // Preleva il prossimo seed da elaborare in modo thread-safe
                std::size_t i = task_index.fetch_add(1, std::memory_order_relaxed);
                if (i >= total_runs) break; // Non ci sono più task
                
                std::uint64_t current_seed = seeds[i];
                
                // Creiamo una copia della misura per questo thread
                MeasureType local_measure = measure_base;
                
                // Esecuzione dell'algoritmo
                MilanoResult current_result = milano_communities(
                                                                 G,
                                                                 local_measure,
                                                                 initial_partition,
                                                                 threshold,
                                                                 max_level,
                                                                 collect_stats,
                                                                 current_seed
                                                                 );
                
                // Valutazione dello score globale
                double current_score = 0.0;
                if (!current_result.cluster_scores.empty()) {
                    for (double score : current_result.cluster_scores) current_score += score;
                } else {
                    current_score = 0.0;
                }
                
                // Aggiornamento on-the-fly del best result LOCALE del thread
                if (current_score > local_best_score) {
                    local_best_score = current_score;
                    local_best_result = std::move(current_result);
                }
            }
            
            // Una volta finiti i task, aggiorniamo il best globale (richiede il lock)
            if (local_best_score > -std::numeric_limits<double>::infinity()) {
                std::lock_guard<std::mutex> lock(best_mutex);
                if (local_best_score > global_best_score) {
                    global_best_score = local_best_score;
                    global_best_result = std::move(local_best_result);
                }
            }
        };
        // Lancio dei thread
        std::vector<std::jthread> threads;
        threads.reserve(n_threads);
        for (std::uint64_t t = 0; t < n_threads; ++t) {
            threads.emplace_back(worker);
        }
        
        // std::jthread esegue il join automaticamente alla distruzione.
        // Chiamando .clear() blocchiamo l'esecuzione principale finché
        // tutti i thread non hanno completato il loro lavoro.
        threads.clear();
        
        return global_best_result;
    }
    
} // namespace Milano

