/**
 * @file community_measure.h
 *
 * Misure incluse:
 * 1. ModularityMeasure
 * 2. PersistenceMeasure
 * 3. PersistenceModularityMeasure
 * 4. PersistenceModularityDensityMeasure (Ripristinato alpha_density originale)
 */

#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <concepts>
#include <stdexcept>
#include "grafi.h"

namespace Milano {
    
    /**
     * @typedef Partition
     */
    using Partition = std::vector<std::vector<std::uint64_t>>;
    
    // ========================================================================
    // GRAPH CONTEXT
    // ========================================================================
    
    struct GraphContext {
        const Graph* graph;
        std::vector<double> node_degrees;
        std::vector<double> self_loops;
        
        GraphContext(const Graph* g, std::vector<double> nd, std::vector<double> sl)
        : graph(g), node_degrees(std::move(nd)), self_loops(std::move(sl)) {}
    };
    
    // ========================================================================
    // C++20 CONCEPT FOR MEASURES
    // ========================================================================
    
    template <typename T>
    concept IsCommunityMeasure = requires(T measure,
                                          const Graph& G,
                                          const GraphContext& ctx,
                                          const Partition& partition,
                                          const std::vector<std::uint64_t>& com,
                                          std::uint64_t u,
                                          std::uint64_t old_com,
                                          std::uint64_t new_com,
                                          double degree_u,
                                          double wt_old,
                                          double wt_new,
                                          double self_loop_u,
                                          double gain) {
        { measure.init_graph(G) } -> std::same_as<void>;
        { measure.start_level(ctx, partition) } -> std::same_as<void>;
        { measure.localScore(G, com) } -> std::convertible_to<double>;
        { measure.globalScore(G, partition) } -> std::convertible_to<double>;
        { measure.current_measure() } -> std::convertible_to<double>;
        { measure.evaluate_move(u, old_com, new_com, degree_u, wt_old, wt_new, self_loop_u) } -> std::convertible_to<double>;
        { measure.update(u, old_com, new_com, degree_u, wt_old, wt_new, self_loop_u, gain) } -> std::same_as<void>;
    };
    
    // ========================================================================
    // MEASURE 1: CLASSIC MILANO MODULARITY
    // ========================================================================
    
    class ModularityMeasure {
    private:
        double resolution_;
        double m_, inv_m_, inv_2m2_res_;
        std::vector<double> Stot_;
        double current_q_;
        
    public:
        explicit ModularityMeasure(double resolution = 1.0)
        : resolution_(resolution), m_(0.0), inv_m_(0.0), inv_2m2_res_(0.0), current_q_(0.0) {}
        
        void init_graph(const Graph& original_G) {
            m_ = original_G.TotalWeight();
            if (m_ > 0.0) {
                inv_m_ = 1.0 / m_;
                inv_2m2_res_ = resolution_ / (2.0 * m_ * m_);
            }
        }
        
        void start_level(const GraphContext& ctx, const Partition& partition) {
            std::uint64_t num_coms = partition.size();
            Stot_.assign(num_coms, 0.0);
            
            for (std::uint64_t i = 0; i < num_coms; ++i) {
                for (auto u : partition[i]) Stot_[i] += ctx.node_degrees[u];
            }
            
            current_q_ = 0.0;
            if (m_ > 0.0) {
                for (std::uint64_t i = 0; i < num_coms; ++i) {
                    if (partition[i].empty()) continue;
                    
                    double e_c = 0.0;
                    for (auto u : partition[i]) {
                        e_c += ctx.self_loops[u];
                        for (const auto& [v, wt] : ctx.graph->EdgesOf(u)) {
                            // C++20: ranges::find su vector
                            if (std::ranges::find(partition[i], v) != partition[i].end()) e_c += wt / 2.0;
                        }
                    }
                    current_q_ += (e_c * inv_m_) - resolution_ * std::pow(Stot_[i] / (2.0 * m_), 2);
                }
            }
        }
        
        [[nodiscard]] double localScore(const Graph& original_G, const std::vector<std::uint64_t>& com) const {
            if (m_ == 0.0 || com.empty()) return 0.0;
            double e_c = 0.0, a_c = 0.0;
            for (std::uint64_t u : com) {
                a_c += original_G.Strength(u);
                for (const auto& [v, wt] : original_G.EdgesOf(u)) {
                    if (std::ranges::find(com, v) != com.end()) e_c += wt / 2.0;
                }
            }
            return (e_c * inv_m_) - resolution_ * std::pow(a_c / (2.0 * m_), 2);
        }
        
        [[nodiscard]] double globalScore(const Graph& original_G, const Partition& partition) const {
            if (m_ == 0.0) return 0.0;
            double q = 0.0;
            for (const auto& com : partition) q += localScore(original_G, com);
            return q;
        }
        
        [[nodiscard]] double current_measure() const { return current_q_; }
        
        [[nodiscard]] inline double evaluate_move([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                                                  double degree_u, double wt_old, double wt_new, [[maybe_unused]] double self_loop_u) const {
            if (m_ == 0.0) return current_q_;
            
            double stot_old_res = Stot_[old_com] - degree_u;
            double remove_cost = -wt_old * inv_m_ + (stot_old_res * degree_u) * inv_2m2_res_;
            
            double stot_new_res = Stot_[new_com];
            double add_gain = wt_new * inv_m_ - (stot_new_res * degree_u) * inv_2m2_res_;
            
            return current_q_ + remove_cost + add_gain;
        }
        
        inline void update([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                           double degree_u, [[maybe_unused]] double wt_old, [[maybe_unused]] double wt_new, [[maybe_unused]] double self_loop_u, double gain) {
            Stot_[old_com] -= degree_u;
            Stot_[new_com] += degree_u;
            current_q_ += gain;
        }
    };
    
    // ========================================================================
    // MEASURE 2: CLASSIC PERSISTENCE
    // ========================================================================
    
    class PersistenceMeasure {
    private:
        std::vector<double> in_wt_;
        std::vector<double> tot_str_;
        double current_q_;
        
    public:
        PersistenceMeasure() : current_q_(0.0) {}
        
        void init_graph([[maybe_unused]] const Graph& original_G) {}
        
        void start_level(const GraphContext& ctx, const Partition& partition) {
            std::uint64_t num_coms = partition.size();
            in_wt_.assign(num_coms, 0.0);
            tot_str_.assign(num_coms, 0.0);
            current_q_ = 0.0;
            
            for (std::uint64_t i = 0; i < num_coms; ++i) {
                if (partition[i].empty()) continue;
                double e_c = 0.0, s_c = 0.0;
                for (auto u : partition[i]) {
                    s_c += ctx.node_degrees[u];
                    e_c += ctx.self_loops[u];
                    for (const auto& [v, wt] : ctx.graph->EdgesOf(u)) {
                        if (std::ranges::find(partition[i], v) != partition[i].end()) e_c += wt / 2.0;
                    }
                }
                in_wt_[i] = e_c;
                tot_str_[i] = s_c;
                if (s_c > 0) current_q_ += (2.0 * e_c) / s_c;
            }
        }
        
        [[nodiscard]] double localScore(const Graph& original_G, const std::vector<std::uint64_t>& C) const {
            if (C.empty()) return 0.0;
            double inner_weights = 0.0, inner_strength = 0.0;
            for (auto v1_it = C.begin(); v1_it != C.end(); ++v1_it) {
                auto v1 = *v1_it;
                inner_strength += original_G.Strength(v1);
                for (auto v2_it = std::next(v1_it); v2_it != C.end(); ++v2_it) {
                    inner_weights += original_G.Weight(v1, *v2_it).second;
                }
            }
            return (inner_strength > 0) ? (2.0 * inner_weights) / inner_strength : 0.0;
        }
        
        [[nodiscard]] double globalScore(const Graph& original_G, const Partition& partition) const {
            double q = 0.0;
            for (const auto& com : partition) q += localScore(original_G, com);
            return q;
        }
        
        [[nodiscard]] double current_measure() const { return current_q_; }
        
        [[nodiscard]] inline double evaluate_move([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                                                  double degree_u, double wt_old, double wt_new, double self_loop_u) const {
            
            double old_A_score = (tot_str_[old_com] > 0) ? (2.0 * in_wt_[old_com]) / tot_str_[old_com] : 0.0;
            double old_B_score = (tot_str_[new_com] > 0) ? (2.0 * in_wt_[new_com]) / tot_str_[new_com] : 0.0;
            
            double new_W_A = in_wt_[old_com] - wt_old - self_loop_u;
            double new_S_A = tot_str_[old_com] - degree_u;
            double new_A_score = (new_S_A > 0) ? (2.0 * new_W_A) / new_S_A : 0.0;
            
            double new_W_B = in_wt_[new_com] + wt_new + self_loop_u;
            double new_S_B = tot_str_[new_com] + degree_u;
            double new_B_score = (new_S_B > 0) ? (2.0 * new_W_B) / new_S_B : 0.0;
            
            return current_q_ - old_A_score - old_B_score + new_A_score + new_B_score;
        }
        
        inline void update([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                           double degree_u, double wt_old, double wt_new, double self_loop_u, double gain) {
            in_wt_[old_com] -= (wt_old + self_loop_u);
            tot_str_[old_com] -= degree_u;
            in_wt_[new_com] += (wt_new + self_loop_u);
            tot_str_[new_com] += degree_u;
            current_q_ += gain;
        }
    };
    
    // ========================================================================
    // MEASURE 3: PERSISTENCE WITH MODULARITY PENALTY
    // ========================================================================
    
    class PersistenceModularityMeasure {
    private:
        std::vector<double> in_wt_;
        std::vector<double> tot_str_;
        double total_graph_strength_;
        double inv_total_graph_strength_;
        double current_q_;
        
        [[nodiscard]] inline double calc_cluster_score(double w_in, double s_tot) const {
            if (s_tot <= 0) return 0.0;
            double p_score = (2.0 * w_in) / s_tot;
            double m_penalty = s_tot * inv_total_graph_strength_;
            return p_score - m_penalty;
        }
        
    public:
        PersistenceModularityMeasure()
        : total_graph_strength_(0.0), inv_total_graph_strength_(0.0), current_q_(0.0) {}
        
        void init_graph(const Graph& original_G) {
            total_graph_strength_ = 2.0 * original_G.TotalWeight();
            if (total_graph_strength_ > 0) {
                inv_total_graph_strength_ = 1.0 / total_graph_strength_;
            }
        }
        
        void start_level(const GraphContext& ctx, const Partition& partition) {
            std::uint64_t num_coms = partition.size();
            in_wt_.assign(num_coms, 0.0);
            tot_str_.assign(num_coms, 0.0);
            current_q_ = 0.0;
            
            for (std::uint64_t i = 0; i < num_coms; ++i) {
                if (partition[i].empty()) continue;
                double e_c = 0.0, s_c = 0.0;
                for (auto u : partition[i]) {
                    s_c += ctx.node_degrees[u];
                    e_c += ctx.self_loops[u];
                    for (const auto& [v, wt] : ctx.graph->EdgesOf(u)) {
                        if (std::ranges::find(partition[i], v) != partition[i].end()) e_c += wt / 2.0;
                    }
                }
                in_wt_[i] = e_c;
                tot_str_[i] = s_c;
                current_q_ += calc_cluster_score(e_c, s_c);
            }
        }
        
        [[nodiscard]] double localScore(const Graph& original_G, const std::vector<std::uint64_t>& C) const {
            if (C.empty()) return 0.0;
            double inner_weights = 0.0, inner_strength = 0.0;
            for (auto v1_it = C.begin(); v1_it != C.end(); ++v1_it) {
                auto v1 = *v1_it;
                inner_strength += original_G.Strength(v1);
                for (auto v2_it = std::next(v1_it); v2_it != C.end(); ++v2_it) {
                    inner_weights += original_G.Weight(v1, *v2_it).second;
                }
            }
            return calc_cluster_score(inner_weights, inner_strength);
        }
        
        [[nodiscard]] double globalScore(const Graph& original_G, const Partition& partition) const {
            double q = 0.0;
            for (const auto& com : partition) q += localScore(original_G, com);
            return q;
        }
        
        [[nodiscard]] double current_measure() const { return current_q_; }
        
        [[nodiscard]] inline double evaluate_move([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                                                  double degree_u, double wt_old, double wt_new, double self_loop_u) const {
            
            double old_A_score = calc_cluster_score(in_wt_[old_com], tot_str_[old_com]);
            double old_B_score = calc_cluster_score(in_wt_[new_com], tot_str_[new_com]);
            
            double new_W_A = in_wt_[old_com] - wt_old - self_loop_u;
            double new_S_A = tot_str_[old_com] - degree_u;
            double new_A_score = calc_cluster_score(new_W_A, new_S_A);
            
            double new_W_B = in_wt_[new_com] + wt_new + self_loop_u;
            double new_S_B = tot_str_[new_com] + degree_u;
            double new_B_score = calc_cluster_score(new_W_B, new_S_B);
            
            return current_q_ - old_A_score - old_B_score + new_A_score + new_B_score;
        }
        
        inline void update([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                           double degree_u, double wt_old, double wt_new, double self_loop_u, double gain) {
            in_wt_[old_com] -= (wt_old + self_loop_u);
            tot_str_[old_com] -= degree_u;
            in_wt_[new_com] += (wt_new + self_loop_u);
            tot_str_[new_com] += degree_u;
            current_q_ += gain;
        }
    };
    
    // ========================================================================
    // MEASURE 4: PERSISTENCE + MODULARITY PENALTY + DENSITY FACTOR
    // ========================================================================
    
    class PersistenceModularityDensityMeasure {
    private:
        double alpha_density_;                ///< Exponent dictating the density bias.
        std::vector<double> in_wt_;
        std::vector<double> tot_str_;
        double total_graph_strength_;
        double inv_total_graph_strength_;
        double current_q_;
        
        [[nodiscard]] inline double calc_cluster_score(double w_in, double s_tot) const {
            if (s_tot <= 0 || total_graph_strength_ <= 0) return 0.0;
            
            double p_score = (2.0 * w_in) / s_tot;
            double m_penalty = s_tot * inv_total_graph_strength_;
            double base_result = p_score - m_penalty;
            double density = std::pow(s_tot * inv_total_graph_strength_, alpha_density_);
            
            return base_result * density;
        }
        
    public:
        // Ripristinato l'argomento alpha_density originale
        explicit PersistenceModularityDensityMeasure(double alpha_density = 0.5)
        : alpha_density_(alpha_density),
        total_graph_strength_(0.0),
        inv_total_graph_strength_(0.0),
        current_q_(0.0) {
            if (alpha_density_ < 0.0 || alpha_density_ > 1.0) {
                throw std::runtime_error("alpha_density must be strictly between 0.0 and 1.0");
            }
        }
        
        void init_graph(const Graph& original_G) {
            total_graph_strength_ = 2.0 * original_G.TotalWeight();
            if (total_graph_strength_ > 0) {
                inv_total_graph_strength_ = 1.0 / total_graph_strength_;
            }
        }
        
        void start_level(const GraphContext& ctx, const Partition& partition) {
            std::uint64_t num_coms = partition.size();
            in_wt_.assign(num_coms, 0.0);
            tot_str_.assign(num_coms, 0.0);
            current_q_ = 0.0;
            
            for (std::uint64_t i = 0; i < num_coms; ++i) {
                if (partition[i].empty()) continue;
                double e_c = 0.0, s_c = 0.0;
                for (auto u : partition[i]) {
                    s_c += ctx.node_degrees[u];
                    e_c += ctx.self_loops[u];
                    for (const auto& [v, wt] : ctx.graph->EdgesOf(u)) {
                        if (std::ranges::find(partition[i], v) != partition[i].end()) e_c += wt / 2.0;
                    }
                }
                in_wt_[i] = e_c;
                tot_str_[i] = s_c;
                current_q_ += calc_cluster_score(e_c, s_c);
            }
        }
        
        [[nodiscard]] double localScore(const Graph& original_G, const std::vector<std::uint64_t>& C) const {
            if (C.empty()) return 0.0;
            double inner_weights = 0.0, inner_strength = 0.0;
            for (auto v1_it = C.begin(); v1_it != C.end(); ++v1_it) {
                auto v1 = *v1_it;
                inner_strength += original_G.Strength(v1);
                for (auto v2_it = std::next(v1_it); v2_it != C.end(); ++v2_it) {
                    inner_weights += original_G.Weight(v1, *v2_it).second;
                }
            }
            return calc_cluster_score(inner_weights, inner_strength);
        }
        
        [[nodiscard]] double globalScore(const Graph& original_G, const Partition& partition) const {
            double q = 0.0;
            for (const auto& com : partition) q += localScore(original_G, com);
            return q;
        }
        
        [[nodiscard]] double current_measure() const { return current_q_; }
        
        [[nodiscard]] inline double evaluate_move([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                                                  double degree_u, double wt_old, double wt_new, double self_loop_u) const {
            
            double old_A_score = calc_cluster_score(in_wt_[old_com], tot_str_[old_com]);
            double old_B_score = calc_cluster_score(in_wt_[new_com], tot_str_[new_com]);
            
            double new_W_A = in_wt_[old_com] - wt_old - self_loop_u;
            double new_S_A = tot_str_[old_com] - degree_u;
            double new_A_score = calc_cluster_score(new_W_A, new_S_A);
            
            double new_W_B = in_wt_[new_com] + wt_new + self_loop_u;
            double new_S_B = tot_str_[new_com] + degree_u;
            double new_B_score = calc_cluster_score(new_W_B, new_S_B);
            
            return current_q_ - old_A_score - old_B_score + new_A_score + new_B_score;
        }
        
        inline void update([[maybe_unused]] std::uint64_t u, std::uint64_t old_com, std::uint64_t new_com,
                           double degree_u, double wt_old, double wt_new, double self_loop_u, double gain) {
            in_wt_[old_com] -= (wt_old + self_loop_u);
            tot_str_[old_com] -= degree_u;
            in_wt_[new_com] += (wt_new + self_loop_u);
            tot_str_[new_com] += degree_u;
            current_q_ += gain;
        }
    };
    
} // namespace Milano

