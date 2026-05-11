#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cstdint>
#include <utility>
#include <cmath>
#include <tuple>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include <R_ext/Parse.h>

#include "eccezioni.h"
#include "rwrapper.h"
#include "random.h"
#include "grafi.h"
#include "community_measure.h"
#include "milano_algorithm.h"


[[noreturn]] static void forwardExceptionToR(const char* message) {
    Rf_error("%s", message);
}

// ***********************************************
// Classe helper per gestire i dati del grafo
// ***********************************************

class GraphData {
public:
    std::vector<std::string> position_to_vertex;
    std::unordered_map<std::string, std::uint64_t> vertex_to_position;
    Milano::EdgeList edge_list;
    
    void buildVertexMaps(SEXP vertex_r) {
        const R_xlen_t vertex_size = Rf_xlength(vertex_r);
        position_to_vertex.resize(vertex_size);
        vertex_to_position.reserve(vertex_size);
        
        for (R_xlen_t k = 0; k < vertex_size; ++k) {
            const char* name = R_CHAR(STRING_ELT(vertex_r, k));
            position_to_vertex[k] = name;
            vertex_to_position.emplace(name, static_cast<std::uint64_t>(k));
        }
    }
    
    
    void loadEdgeList(SEXP edge_list_r, SEXP weights_r) {
        SEXP dim_r = Rf_getAttrib(edge_list_r, R_DimSymbol);
        const R_xlen_t nrow = INTEGER(dim_r)[0];
        const R_xlen_t ncol = INTEGER(dim_r)[1];
        
        if (ncol != 2) {
            throw Milano::eccezioni("Error edge list: expected 2 columns, got " + std::to_string(ncol));
        }
        
        const double* weights_ptr = REAL(weights_r);
        edge_list.clear();
        edge_list.reserve(nrow);
        
        for (R_xlen_t k = 0; k < nrow; ++k) {
            const char* v1_str = R_CHAR(STRING_ELT(edge_list_r, k));
            const char* v2_str = R_CHAR(STRING_ELT(edge_list_r, k + nrow));
            
            auto it1 = vertex_to_position.find(v1_str);
            auto it2 = vertex_to_position.find(v2_str);
            
            if (it1 == vertex_to_position.end() || it2 == vertex_to_position.end()) {
                throw Milano::eccezioni("Vertex not found in vertex list");
            }
            
            edge_list.emplace_back(it1->second, it2->second, weights_ptr[k]);
        }
    }
    
    Milano::UGraph buildGraph() const {
        return Milano::UGraph(position_to_vertex.size(), edge_list);
    }
    
    std::size_t vertexCount() const {
        return position_to_vertex.size();
    }
};

// ***********************************************
// Helper per estrarre H0 da SEXP
// Ritorna: -1.0 se NULL (persistence probability)
//           0.0 se null-adjusted persistence
//          >0.0 se null-adjusted persistence density
// ***********************************************

static double extractH0(SEXP h0_r) {
    if (Rf_isNull(h0_r) || Rf_xlength(h0_r) == 0) {
        return -1.0;
    }
    return REAL(h0_r)[0];
}

// ***********************************************
// Convertitori Membership R -> C++20 Partition
// ***********************************************

static std::vector<std::uint64_t> loadMembershipVector(SEXP membership_r, std::size_t expected_size) {
    if (Rf_isNull(membership_r) || Rf_xlength(membership_r) == 0) return {};
    
    std::vector<std::uint64_t> membership(expected_size);
    const int* membership_ptr = INTEGER(membership_r);
    
    for (std::size_t p = 0; p < expected_size; ++p) {
        membership[p] = static_cast<std::uint64_t>(membership_ptr[p]);
    }
    return membership;
}



static Milano::Partition buildInitialPartition(const std::vector<std::uint64_t>& membership) { // Prefissato
    if (membership.empty()) return {};
    
    std::uint64_t max_cluster_id = 0;
    for (auto c : membership) {
        if (c > 0 && (c - 1) > max_cluster_id) {
            max_cluster_id = c - 1; // R è 1-based index
        }
    }
    
    Milano::Partition partition(max_cluster_id + 1);
    for (std::uint64_t k = 0; k < membership.size(); ++k) {
        if (membership[k] > 0) {
            partition[membership[k] - 1].push_back(k);
        }
    }
    return partition;
}

// ***********************************************
// Concept-Based Dispatcher
// ***********************************************

template <typename Callable>
decltype(auto) dispatch_measure(double h0, Callable&& func) {
    if (h0 < 0.0) {
        Milano::PersistenceMeasure m;
        return func(m);
    } else if (h0 == 0.0) {
        Milano::PersistenceModularityMeasure m;
        return func(m);
    } else {
        Milano::PersistenceModularityDensityMeasure m(h0);
        return func(m);
    }
}



// ***********************************************
// MAIN R-EXPORTS (extern "C")
// ***********************************************

extern "C" {
    SEXP global_persistence(SEXP vertex_r, SEXP edge_list_r, SEXP weights_r,
                           SEXP membership_r, SEXP h0_r) {
        RProtectGuard guard;
        try {
            // ==========================================================
            // FASE 1: Costruzione dei dati C++ con i dati R
            // ==========================================================
            GraphData data;
            data.buildVertexMaps(vertex_r);
            data.loadEdgeList(edge_list_r, weights_r);
            
            Milano::UGraph grafo = data.buildGraph(); // Prefissato
            auto membership_vec = loadMembershipVector(membership_r, data.vertexCount());
            Milano::Partition partition = buildInitialPartition(membership_vec); // Prefissato
            const double h0 = extractH0(h0_r);
           
            // ==========================================================
            // FASE 2: Calcolo C++ Puro
            // ==========================================================
            auto [total_value, cluster_values] = dispatch_measure(h0, [&](auto& measure) {
                measure.init_graph(grafo);
                
                std::vector<double> scores(partition.size());
                double tv = 0.0;
                for (std::size_t i = 0; i < partition.size(); ++i) {
                    double s = measure.localScore(grafo, partition[i]);
                    scores[i] = s;
                    if (std::isfinite(s)) tv += s;
                }
                return std::make_pair(tv, std::move(scores));
            });
            
            // ==========================================================
            // FASE 3: Costruzione oggetto R risultato
            // ==========================================================
            SEXP result_cm_r = RHelper::createRealVector(cluster_values.size(), guard);
            double* result_cm_ptr = REAL(result_cm_r);
            for (std::size_t i = 0; i < cluster_values.size(); ++i) {
                result_cm_ptr[i] = cluster_values[i];
            }
            
            SEXP result_value_r = RHelper::createScalarReal(total_value, guard);
            SEXP result_r = RHelper::createNamedList(2, guard);
            RHelper::setListElement(result_r, 0, result_value_r, "score");
            RHelper::setListElement(result_r, 1, result_cm_r, "clusters_value");
            
            return result_r;
             
        } catch(std::exception &ex) {
            forwardExceptionToR(ex.what());
        } catch(...) {
            forwardExceptionToR("C++ exception (unknown reason)");
        }
    }
    
    SEXP local_persistence(SEXP vertex_r, SEXP edge_list_r, SEXP weights_r,
                          SEXP cluster_r, SEXP h0_r) {
        RProtectGuard guard;
        
        try {
            // ==========================================================
            // FASE 1: Costruzione dei dati C++ con i dati R
            // ==========================================================
            GraphData data;
            data.buildVertexMaps(vertex_r);
            data.loadEdgeList(edge_list_r, weights_r);
            
            Milano::UGraph grafo = data.buildGraph();
            
            const R_xlen_t cluster_size = Rf_xlength(cluster_r);
            const int* cluster_ptr = INTEGER(cluster_r);
            
            std::vector<std::uint64_t> target_cluster;
            target_cluster.reserve(cluster_size);
            
            for (R_xlen_t p = 0; p < cluster_size; ++p) {
                if (cluster_ptr[p] == 1) {
                    target_cluster.push_back(static_cast<std::uint64_t>(p));
                }
            }
            
            const double h0 = extractH0(h0_r);
            
            // ==========================================================
            // FASE 2: Calcolo C++ Puro
            // ==========================================================
            double s = dispatch_measure(h0, [&](auto& measure) {
                measure.init_graph(grafo);
                return measure.localScore(grafo, target_cluster);
            });
            
            // ==========================================================
            // FASE 3: Costruzione oggetto R risultato
            // ==========================================================
            return RHelper::createScalarReal(s, guard);
            
        } catch (const std::exception& ex) {
            forwardExceptionToR(ex.what());
        } catch (...) {
            forwardExceptionToR("C++ exception (unknown reason)");
        }
    }
    
    SEXP cluster_milano(SEXP vertex_r, SEXP edge_list_r, SEXP weights_r,
                       SEXP membership_r, SEXP h0_r, SEXP seed_r, SEXP tol_r, SEXP max_level_r) {
    
        RProtectGuard guard;
        
        try {
            // ==========================================================
            // FASE 1: Costruzione dei dati C++ con i dati R
            // ==========================================================
            GraphData data;
            data.buildVertexMaps(vertex_r);
            data.loadEdgeList(edge_list_r, weights_r);
            
            const auto vertex_size = data.vertexCount();
            Milano::UGraph grafo = data.buildGraph();
            
            auto membership_vec = loadMembershipVector(membership_r, vertex_size);
            Milano::Partition start_partition = buildInitialPartition(membership_vec);
            
            std::uint32_t seed_rnd;
            if (Rf_xlength(seed_r) > 0) {
                seed_rnd = static_cast<std::uint32_t>(INTEGER(seed_r)[0]);
            } else {
                seed_rnd = static_cast<std::uint32_t>(Milano::Random::GENERATORE_SEED_RANDOM.RndNextInt(0, 0xFFFFFFFF));
            }
            
            double threshold;
            if (Rf_isNull(tol_r) || Rf_xlength(tol_r) == 0) {
                // Calcolo dinamico: 1 ordine di grandezza INFERIORE al reciproco
                // della Total Strength. (Total Strength = 2 * Total Weight)
                double total_strength = grafo.TotalWeight() * 2.0;
                if (total_strength > 0.0) {
                    threshold = 1.0 / (10.0 * total_strength);
                    // Limita il threshold per risolvere entrambi i problemi:
                    // 1. Max 1e-6 : Evita che sui grafi piccoli il threshold sia troppo alto (bloccando mosse valide)
                    // 2. Min 1e-10: Evita che sui grafi grandi scenda sotto la soglia del rumore numerico (ferma il ping-pong)
                    threshold = std::clamp(threshold, 1e-10, 1e-6);
                } else {
                    threshold = 1e-7; // Fallback per grafi vuoti
                }
                
            } else {
                // Valore esplicito scelto dall'utente in R
                threshold = REAL(tol_r)[0];
            }
            
            std::uint32_t max_level = 0; // Se NULL o vuoto, il default è 0 (illimitato)
            if (!Rf_isNull(max_level_r) && Rf_xlength(max_level_r) > 0) {
                max_level = static_cast<std::uint32_t>(INTEGER(max_level_r)[0]);
            }
            
            const double h0 = extractH0(h0_r);
            
            // ==========================================================
            // FASE 2: Calcolo C++ Puro
            // ==========================================================
            auto [milano_res, best_value] = dispatch_measure(h0, [&](auto& measure) {
                Milano::MilanoResult res = Milano::milano_communities(
                                                grafo, measure, start_partition, threshold, max_level, false, seed_rnd
                                            );
                double bv = measure.globalScore(grafo, res.communities);
                return std::make_pair(std::move(res), bv);
            });
            
            // ==========================================================
            // FASE 3: Costruzione oggetto R risultato
            // ==========================================================
            SEXP membership_result_r = RHelper::createIntVector(vertex_size, guard);
            int* membership_result_ptr = INTEGER(membership_result_r);
            
            int cluster_label = 1;
            for (const auto& comm : milano_res.communities) {
                if (comm.empty()) continue;
                for (auto v : comm) {
                    membership_result_ptr[v] = cluster_label;
                }
                ++cluster_label;
            }
            
            SEXP best_value_r = RHelper::createScalarReal(best_value, guard);
            SEXP seed_result_r = RHelper::createScalarString(std::to_string(seed_rnd), guard);
            
            SEXP result_r = RHelper::createNamedList(3, guard);
            RHelper::setListElement(result_r, 0, membership_result_r, "membership");
            RHelper::setListElement(result_r, 1, best_value_r, "score");
            RHelper::setListElement(result_r, 2, seed_result_r, "seed");
            
            return result_r;
            
        } catch (const std::exception& ex) {
            forwardExceptionToR(ex.what());
        } catch (...) {
            forwardExceptionToR("C++ exception (unknown reason)");
        }
    }
    
} // extern "C"


// ***********************************************
// Registrazione metodi per R
// ***********************************************

extern "C" {
    
    static const R_CallMethodDef CallEntries[] = {
        {"cluster_milano",     (DL_FUNC) &cluster_milano,     8},
        {"global_persistence", (DL_FUNC) &global_persistence, 5},
        {"local_persistence",  (DL_FUNC) &local_persistence,  5},
        {NULL, NULL, 0}
    };
    
    void R_init_persistence(DllInfo *dll) {
        R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
        R_useDynamicSymbols(dll, FALSE);
        R_forceSymbols(dll, TRUE);
    }
    
} // extern "C"














