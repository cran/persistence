#ifndef random_hpp
#define random_hpp

#include <memory>
#include <random>
#include <list>
#include <vector>


// ***********************************************
// ***********************************************
// ***********************************************
class Random {
protected:
    std::uint_fast64_t seed;
public:
    virtual ~Random() {};
    std::uint_fast64_t Seed() {return seed;};
    virtual void Restart() = 0;
    virtual double RndNext() const = 0;
    virtual std::uint_fast64_t RndNextInt(std::uint_fast64_t, std::uint_fast64_t) const = 0;
    virtual std::shared_ptr<std::list<std::uint_fast64_t>> RndSample(std::uint_fast64_t, std::uint_fast64_t, std::uint_fast64_t) const = 0;
    virtual void RndShuffle(std::vector<std::uint_fast64_t>&) const = 0;

};

// ***********************************************
// ***********************************************
// ***********************************************

class RandomUni: public Random {
private:
    std::shared_ptr<std::mt19937_64> generatore;
public:
    RandomUni(std::uint_fast64_t s);
    ~RandomUni() {}
    void Restart();
    double RndNext() const;
    std::uint_fast64_t RndNextInt(std::uint_fast64_t min, std::uint_fast64_t max) const;
    std::shared_ptr<std::list<std::uint_fast64_t>> RndSample(std::uint_fast64_t min, std::uint_fast64_t max, std::uint_fast64_t k=1) const;
    void RndShuffle(std::vector<std::uint_fast64_t>&) const;
public:
    static std::uint_fast64_t START_SEED;
    static RandomUni GENERATORE_SEED_RANDOM;
};

#endif /* random_hpp */
