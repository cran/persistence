/**
 * @file grafi.h
 */

#pragma once

#include <cstdint>
#include <vector>
#include <tuple>
#include <utility>
#include <memory>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace Milano {
    
    using EdgeList = std::vector<std::tuple<std::uint64_t, std::uint64_t, double>>;
    using AdjacencyMap = std::vector<std::pair<std::uint64_t, double>>;
    using AdjacencyStructure = std::vector<AdjacencyMap>;
    
    class Graph {
    public:
        virtual ~Graph() = default;
        
        virtual std::uint64_t Size() const = 0;
        virtual std::uint64_t EdgeCount() const = 0;
        virtual bool IsDirected() const = 0;
        
        virtual bool HasEdge(std::uint64_t v1, std::uint64_t v2) const = 0;
        virtual std::pair<bool, double> Weight(std::uint64_t v1, std::uint64_t v2) const = 0;
        virtual double GetWeight(std::uint64_t v1, std::uint64_t v2) const = 0;
        
        virtual std::uint64_t Degree() const = 0;
        virtual std::uint64_t Degree(std::uint64_t v) const = 0;
        virtual double Strength() const = 0;
        virtual double Strength(std::uint64_t v) const = 0;
        virtual double TotalWeight() const = 0;
        virtual double OutStrength(std::uint64_t v) const = 0;
        virtual double InStrength(std::uint64_t v) const = 0;
        virtual double EdgeScale() const = 0;
        
        virtual const AdjacencyMap& EdgesOf(std::uint64_t v) const = 0;
        virtual AdjacencyMap& EdgesOf(std::uint64_t v) = 0;
        
        virtual bool Traversal() const = 0;
        virtual void ForEachEdge(const std::function<void(std::uint64_t, std::uint64_t, double)>& f) const = 0;
        virtual void ForEachInEdge(std::uint64_t v, const std::function<void(std::uint64_t, double)>& f) const = 0;
        
        virtual std::unique_ptr<Graph> Make(std::uint64_t num_vertices, EdgeList& edges) const = 0;
    };
    
    class UGraph final : public Graph {
    private:
        AdjacencyStructure __neighbors;
        std::vector<std::uint64_t> __degree;
        std::vector<double> __strength;
        std::uint64_t __total_degree;
        double __total_strength;
        std::uint64_t __edge_count;
        
    public:
        UGraph(std::uint64_t num_vertices, const EdgeList& edges)
        : __total_degree(0), __total_strength(0.0), __edge_count(0) {
            __neighbors.resize(num_vertices);
            __degree.assign(num_vertices, 0);
            __strength.assign(num_vertices, 0.0);
            
            for (const auto& [v1, v2, w] : edges) {
                if (v1 == v2) continue;
                if (v1 < num_vertices) { __degree[v1]++; __degree[v2]++; }
            }
            
            for (std::uint64_t i = 0; i < num_vertices; ++i) {
                __neighbors[i].reserve(__degree[i]);
                __degree[i] = 0;
            }
            
            for (const auto& [v1, v2, weight] : edges) {
                double p = std::isnan(weight) ? 1.0 : weight;
                __neighbors[v1].emplace_back(v2, p);
                __neighbors[v2].emplace_back(v1, p);
                
                __degree[v1]++;
                __degree[v2]++;
                __strength[v1] += p;
                __strength[v2] += p;
                
                __total_degree += 2;
                __total_strength += (2.0 * p);
                __edge_count++;
            }
        }
        
        std::uint64_t Size() const override { return __neighbors.size(); }
        std::uint64_t EdgeCount() const override { return __edge_count; }
        bool IsDirected() const override { return false; }
        
        bool HasEdge(std::uint64_t v1, std::uint64_t v2) const override {
            for (const auto& pair : __neighbors[v1]) {
                if (pair.first == v2) return true;
            }
            return false;
        }
        
        std::pair<bool, double> Weight(std::uint64_t v1, std::uint64_t v2) const override {
            for (const auto& pair : __neighbors[v1]) {
                if (pair.first == v2) return {true, pair.second};
            }
            return {false, 0.0};
        }
        
        double GetWeight(std::uint64_t v1, std::uint64_t v2) const override {
            auto result = Weight(v1, v2);
            if (result.first) return result.second;
            throw std::runtime_error("Edge not found.");
        }
        
        std::uint64_t Degree() const override { return __total_degree; }
        std::uint64_t Degree(std::uint64_t v) const override { return __degree[v]; }
        double Strength() const override { return __total_strength; }
        double Strength(std::uint64_t v) const override { return __strength[v]; }
        double TotalWeight() const override { return __total_strength / 2.0; }
        double OutStrength(std::uint64_t v) const override { return __strength[v]; }
        double InStrength(std::uint64_t v) const override { return __strength[v]; }
        double EdgeScale() const override { return 2.0; }
        
        const AdjacencyMap& EdgesOf(std::uint64_t v) const override { return __neighbors[v]; }
        AdjacencyMap& EdgesOf(std::uint64_t v) override { return __neighbors[v]; }
        
        bool Traversal() const override { return true; } // Implementazione semplificata (connessione)
        
        void ForEachEdge(const std::function<void(std::uint64_t, std::uint64_t, double)>& f) const override {
            for (std::uint64_t u = 0; u < Size(); ++u) {
                for (const auto& [v, w] : __neighbors[u]) {
                    f(u, v, w);
                }
            }
        }
        
        void ForEachInEdge(std::uint64_t v, const std::function<void(std::uint64_t, double)>& f) const override {
            for (const auto& [u, w] : __neighbors[v]) {
                f(u, w);
            }
        }
        
        std::unique_ptr<Graph> Make(std::uint64_t num_vertices, EdgeList& edges) const override {
            return std::make_unique<UGraph>(num_vertices, edges);
        }
    };
    
} // namespace Milano

