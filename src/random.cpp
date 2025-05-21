#include <iostream>
#include <chrono>
#include <memory>
#include <random>
#include <limits>
#include <list>
#include <vector>
#include <set>

#include "random.h"

std::uint_fast64_t RandomUni::START_SEED = (std::uint_fast64_t) (std::chrono::high_resolution_clock::now().time_since_epoch().count());
RandomUni RandomUni::GENERATORE_SEED_RANDOM(RandomUni::START_SEED);

// ***********************************************
// ***********************************************
// ***********************************************

RandomUni::RandomUni(std::uint_fast64_t s) {
    seed = s;
    generatore = std::make_shared<std::mt19937_64>(s);
    /*if (this->seed == nullptr) {
        this->seed = std::make_shared<unsigned int>(Random::GENERATORE_SEED_RANDOM->RndNextInt(0, std::numeric_limits<unsigned int>::max()));
    }
    */
}

// ***********************************************
// ***********************************************
// ***********************************************

void RandomUni::Restart() {
    generatore->seed(seed);
}

// ***********************************************
// ***********************************************
// ***********************************************

double RandomUni::RndNext() const {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double risultato = uniform(*generatore);
    return risultato;
}

// ***********************************************
// ***********************************************
// ***********************************************

std::uint_fast64_t RandomUni::RndNextInt(std::uint_fast64_t min, std::uint_fast64_t max) const {
    std::uniform_int_distribution<std::uint_fast64_t> dis(min, max);

    std::uint_fast64_t risultato = dis(*generatore);
    return risultato;
}

// ***********************************************
// ***********************************************
// ***********************************************

std::shared_ptr<std::list<std::uint_fast64_t>> RandomUni::RndSample(std::uint_fast64_t min, std::uint_fast64_t max, std::uint_fast64_t k) const {
    auto result = std::make_shared<std::list<std::uint_fast64_t>>();
    std::vector<std::uint_fast64_t> polulation_vett(max - min + 1);
    std::iota(polulation_vett.begin(), polulation_vett.end(), 0);
    std::set<std::uint_fast64_t> polulation_set(polulation_vett.begin(), polulation_vett.end());

    while (k > 0 && polulation_set.size() > 0) {
        std::discrete_distribution<> dis(polulation_set.begin(), polulation_set.end());
        std::uint_fast64_t r = dis(*this->generatore);
        polulation_set.erase(r);
        result->push_back(r);
        --k;
    }

    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

void RandomUni::RndShuffle(std::vector<std::uint_fast64_t>& v) const {
    std::shuffle(v.begin(), v.end(), *this->generatore);
}
// ***********************************************
// ***********************************************
// ***********************************************
