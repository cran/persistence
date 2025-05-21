#ifndef modularity_h
#define modularity_h

#include <memory>
#include <set>
#include <vector>
#include <list>

class MatriceDouble;
class MatriceBool;
class Random;
class CommunityMeasure;
class CommunityMeasureV2;
class UGraph;

enum class GCType {normal, no_rnd, only_rnd};

std::shared_ptr<std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>>
RenameEdges(std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>& old_edges_list, std::vector<std::uint_fast64_t>& vertex_rename);

std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>
RenamePartition(std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>& old_partition,
                std::map<std::uint_fast64_t, std::uint_fast64_t>& vertex_rename);


void louvainMethod(UGraph& matrice_adiacenza,
                   CommunityMeasure& misura,
                   std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> start_clusters,
                   std::shared_ptr<std::vector<double>> pi_value,
                   std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>& final_clusters,
                   std::pair<double, std::vector<double>>& fo_best);


std::string communities_to_string(std::map<std::uint_fast64_t, std::set<std::uint_fast64_t>>& c, std::uint_fast64_t n);

template <typename T>
std::string to_string(std::vector<T>& s) {
    std::string r = "";
    for (auto v : s) {
        r += std::to_string(v) + " ";
    }

    return r;
}



#endif /* modularity_h */
