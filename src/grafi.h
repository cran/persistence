#ifndef Grafo_h
#define Grafo_h

#include <cstdint>
#include <list>
#include <set>
#include <utility>
#include <map>
#include <vector>
#include <stack>
#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t
#include <cmath>

#include "eccezioni.h"

//************************************
//************************************
//************************************

class UGraph {
private:
    using ESetType = std::map<std::uint_fast64_t, std::map<std::uint_fast64_t, double>>;
    ESetType __edge_set;
    std::vector<std::uint_fast64_t> __degree;
    std::uint_fast64_t __total_degree;
    std::vector<double> __strength;
public:
    UGraph(std::uint_fast64_t numero_vertici, std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>& archi) {
        for (std::uint_fast64_t v = 0; v < numero_vertici; ++v) {
            __edge_set.insert(std::make_pair(v, std::map<std::uint_fast64_t, double>()));
        }

        bool nan_found = false;
        bool nan_all = true;
        for (auto l : archi) {
            auto v1 = std::get<0>(l);
            auto v2 = std::get<1>(l);
            auto p = std::get<2>(l);
            if (v1 == v2) {
                throw_line("Constructor: wrong edge!");
            }
            if (!std::isnan(p) && p <= 0) {
                throw_line("Constructor: wrong weight!");
            }
            if (std::isnan(p)) {
                p = 1.0;
                nan_found = true;
            } else {
                nan_all = false;
            }
            auto& ir1 = __edge_set.at(v1);
            auto& ir2 = __edge_set.at(v2);
            ir1.insert(std::make_pair(v2, p));
            ir2.insert(std::make_pair(v1, p));
        }

        if (nan_found && !nan_all) {
            throw_line("Constructor: some weight are NaN but not all!");
        }
        BuildDegree();
        BuildStrength();
    }

    double Size() const {
        return __edge_set.size();
    }

    std::pair<bool, double> Weight(std::uint_fast64_t v1, std::uint_fast64_t v2) const {
        if (v1 == v2) {
            return std::make_pair(false, 0);
        }
        const auto it = __edge_set.at(v1).find(v2);
        if (it == __edge_set.at(v1).end()) {
            return std::make_pair(false, 0);
        }
        return std::make_pair(true, it->second);
    }

    bool Traversal() const {
        std::set<std::uint_fast64_t> check_set;
        for (std::uint_fast64_t i = 0; i < __edge_set.size(); ++i) {
            check_set.insert(check_set.end(), i);
        }

        return Traversal(check_set);
    }

    bool Traversal(const std::set<std::uint_fast64_t>& check_set) const {
        std::uint_fast64_t start = *(check_set.begin());

        std::stack<std::uint_fast64_t> depth_first_search_stack;

        std::set<std::uint_fast64_t> discovered;
        depth_first_search_stack.push(start);
        while (!depth_first_search_stack.empty()) {
            std::uint_fast64_t v = depth_first_search_stack.top();
            depth_first_search_stack.pop();
            if (discovered.find(v) == discovered.end()) {
                discovered.insert(v);
                for (auto w : check_set) {
                    if ((v != w) && Weight(v, w).first) {
                        depth_first_search_stack.push(w);
                    }
                }
            }
        }

        if (check_set.size() != discovered.size()) {
            return false;
        }
        for (auto v : check_set) {
            if (discovered.find(v) == discovered.end()) {
                return false;
            }
        }
        return true;
    }

    std::uint_fast64_t Degree() const {
        return __total_degree;
    }

    std::uint_fast64_t Degree(std::uint_fast64_t v) const {
        return __degree.at(v);
    }

    double Strength(std::uint_fast64_t v) const {
        return __strength.at(v);
    }

    std::map<std::uint_fast64_t, double>& EdgesOf(std::uint_fast64_t v) {
        return __edge_set.at(v);
    }


    struct EdgeIterator {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using first             = ESetType::iterator;
        using second            = std::map<std::uint_fast64_t, double>::iterator;
        using value_type        = std::pair<std::uint_fast64_t, std::uint_fast64_t>;
        using pointer           = value_type*;  // or also value_type*
        using reference         = value_type&;  // or also value_type&

        EdgeIterator(UGraph& g, first f, second s) :__first(f), __second(s),  __graph(g) {
            setVal();
        }

        reference operator*() { return edge; }
        pointer operator->() { return &edge; }

        // Prefix increment
        EdgeIterator& operator++() {
            ++__second;
            while (__second == __first->second.end()) {
                ++__first;
                if (__first != __graph.__edge_set.end()) {
                    __second = __first->second.begin();
                } else {
                    break;
                }
            }
            setVal();
            return *this;
        }

        EdgeIterator operator++(int) { EdgeIterator tmp = *this; ++(*this); return tmp; }

        friend bool operator== (const EdgeIterator& a, const EdgeIterator& b) {
            return ((a.edge.first == std::numeric_limits<std::uint_fast64_t>::max() && b.edge.first == std::numeric_limits<std::uint_fast64_t>::max())
            || (a.edge.first == b.edge.first && a.edge.second == b.edge.second));

        };
        friend bool operator!= (const EdgeIterator& a, const EdgeIterator& b) {
            return !(a == b);
        };

        value_type edge;
    private:
        first __first;
        second __second;
        UGraph& __graph;

        void setVal() {
            if (__first != __graph.__edge_set.end()) {
                edge.first = __first->first;
                edge.second = __second->first;
            } else {
                edge.first = std::numeric_limits<std::uint_fast64_t>::max();
                edge.second = std::numeric_limits<std::uint_fast64_t>::max();
            }
        }
    };

    EdgeIterator eBegin() {
        return EdgeIterator(*this, __edge_set.begin(), __edge_set.begin()->second.begin());
    }
    EdgeIterator eEnd()   {
        return EdgeIterator(*this, __edge_set.end(), __edge_set.begin()->second.end());
    }


private:
    void BuildDegree() {
        __degree.assign(__edge_set.size(), 0);
        __total_degree = 0;
        for (auto& v1_iter : __edge_set) {
            auto v1 = v1_iter.first;
            __degree.at(v1) += v1_iter.second.size();
            __total_degree += __degree.at(v1);
        }
        return;
    }

    void BuildStrength() {
        __strength.assign(Size(), 0.0);
        for (auto& v1_iter : __edge_set) {
            auto v1 = v1_iter.first;
            for (auto& v2_iter : v1_iter.second) {
              __strength.at(v1) += v2_iter.second;
                //auto v2 = v2_iter.first;
                //__strength.at(v2) += v2_iter.second;
            }
        }
    }
};


#endif /* Grafo_h */
