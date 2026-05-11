/**
 * @file random.h
 * @brief Header-only, HPC-optimised random number generation utilities (C++20).
 *
 * @details
 * Integrazione C++20: Utilizzo di std::ranges per un codice più pulito e veloce.
 * Wrapper zero-overhead attorno a `std::mt19937_64`.
 */

#pragma once

#include <random>
#include <vector>
#include <cstdint>
#include <chrono>
#include <unordered_set>
#include <numeric>
#include <algorithm>
#include <ranges>
#include <cstddef>

#include "eccezioni.h"

namespace Milano {
    
    class Random final {
    private:
        std::uint64_t seed_;
        mutable std::mt19937_64 generatore_;
        
        static std::uint64_t genera_seed_robusto() noexcept {
            std::random_device rd;
            auto t = std::chrono::high_resolution_clock::now().time_since_epoch().count();
            return static_cast<std::uint64_t>(rd()) ^ static_cast<std::uint64_t>(t);
        }
        
    public:
        static Random GENERATORE_SEED_RANDOM;
        
        explicit Random(std::uint64_t s) : seed_(s), generatore_(s) {}
        
        [[nodiscard]] inline std::uint64_t Seed() const noexcept { return seed_; }
        
        inline void Restart() { generatore_.seed(seed_); }
        
        [[nodiscard]] inline double RndNext() const {
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            return uniform(generatore_);
        }
        
        [[nodiscard]] inline std::uint64_t RndNextInt(std::uint64_t min, std::uint64_t max) const {
            std::uniform_int_distribution<std::uint64_t> dis(min, max);
            return dis(generatore_);
        }
        
        [[nodiscard]] inline std::vector<std::uint64_t> RndSample(
                          std::uint64_t min, std::uint64_t max, std::uint64_t k = 1) const
        {
            if (min > max) throw Milano::eccezioni("RndSample: min > max");
            
            std::uint64_t range = max - min + 1;
            if (k > range) k = range;
            if (k == 0) return {};
            
            std::vector<std::uint64_t> result;
            result.reserve(k);
            
            if (k > range / 4 && range < 5000000) {
                std::vector<std::uint64_t> population(range);
                std::iota(population.begin(), population.end(), min);
                std::ranges::shuffle(population, generatore_);
                result.assign(population.begin(), population.begin() + static_cast<std::ptrdiff_t>(k));
            } else {
                std::unordered_set<std::uint64_t> picked;
                picked.reserve(k);
                std::uniform_int_distribution<std::uint64_t> dis(min, max);
                while (picked.size() < k) picked.insert(dis(generatore_));
                result.assign(picked.begin(), picked.end());
            }
            return result;
        }
        
        inline void RndShuffle(std::vector<std::uint64_t>& v) const {
            std::ranges::shuffle(v, generatore_);
        }
        
        inline std::mt19937_64& getEngine() const { return generatore_; }
    };
    
    inline Random Random::GENERATORE_SEED_RANDOM{Random::genera_seed_robusto()};
    
} // namespace Milano

