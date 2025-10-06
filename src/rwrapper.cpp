#include <typeinfo>
#include <limits>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <memory>
#include <numeric>
#include <ctime>
#include <chrono>
#include <variant>
#include <thread>


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
#include "matriciCommunity.h"
#include "grafi.h"
#include "communityMeasure.h"
#include "louvainMethod.h"


// ***********************************************
// ***********************************************
// ***********************************************

SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);
    auto len = Rf_length(names);
    for (int i = 0; i < len; i++)
        if(strcmp(R_CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
        }
    return elmt;
}

// ***********************************************
// ***********************************************
// ***********************************************

inline void forward_exception_to_r(const std::string message) {
    int conta_protect = 0;

    auto stop_sym  = Proteggi(Rf_install("stop"), conta_protect) ;
    auto message_r = Proteggi(Rf_allocVector(STRSXP, 1), conta_protect);
    SET_STRING_ELT(message_r, 0, Rf_mkChar(message.c_str()));
    auto expr = Proteggi(Rf_lang2(stop_sym, message_r), conta_protect);
    Rf_eval(expr, R_GlobalEnv);
    if (conta_protect > 0) UNPROTECT(conta_protect);
}

template <class T>
void chanFinalizer(SEXP ptr) {
  if(!R_ExternalPtrAddr(ptr)) return;
  std::cout << "chanFinalizer" << std::endl;
  T *x = static_cast<T*>(R_ExternalPtrAddr(ptr));
  delete x;
  R_ClearExternalPtr(ptr);
}


//************************************
//************************************
//************************************

void LoadEdgeListR(SEXP edge_list_r,
                   SEXP weights_r,
                   std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>>& edge_list,
                   std::map<std::string, std::uint_fast64_t>& vertici_to_position) {
    auto edge_list_dim_r = Rf_getAttrib(edge_list_r, R_DimSymbol);
    long edge_list_nrow = INTEGER(edge_list_dim_r)[0];
    long edge_list_ncol = INTEGER(edge_list_dim_r)[1];

    if (edge_list_ncol != 2) {
        std::string err_str = "Error edge list: " + std::to_string(edge_list_ncol);
        throw_line(err_str);
    }
    std::list<std::tuple<std::string, std::string, double>> dati;

    for (long k = 0; k < edge_list_nrow; ++k) {
        std::string v1 = R_CHAR(STRING_ELT(edge_list_r, k + edge_list_nrow * 0));
        std::string v2 = R_CHAR(STRING_ELT(edge_list_r, k + edge_list_nrow * 1));
        double peso = REAL(weights_r)[k];
        edge_list.push_back({vertici_to_position[v1], vertici_to_position[v2], peso});
    }
}

extern "C" {
    SEXP globalPersistence(SEXP vertex_r, SEXP edge_list_r, SEXP weights_r, SEXP membership_r, SEXP h0_r) {
        int conta_protect = 0;
        SEXP result_r = R_NilValue;
        try {
            long vertex_size = Rf_length(vertex_r);
            std::vector<std::string> position_to_vertici(vertex_size);
            std::map<std::string, std::uint_fast64_t> vertici_to_position;
            for (long k = 0; k < vertex_size; ++k) {
                std::string e = R_CHAR(STRING_ELT(vertex_r, k));
                position_to_vertici.at((std::uint_fast64_t) k) = e;
                vertici_to_position[e] = (std::uint_fast64_t) k;
            }
            std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>> edge_list;
            LoadEdgeListR(edge_list_r, weights_r, edge_list, vertici_to_position);

            auto grafo = UGraph(position_to_vertici.size(), edge_list);

            std::vector<std::uint_fast64_t> membership(vertex_size);
            std::uint_fast64_t membership_size = Rf_length(membership_r);
            std::set<std::uint_fast64_t> used;

            for (std::uint_fast64_t p = 0; p < membership_size; ++p) {
                std::uint_fast64_t value = (std::uint_fast64_t) INTEGER(membership_r)[p];
                membership.at(p) = value;
            }

            bool h0 = LOGICAL(h0_r)[0];
            std::shared_ptr<CommunityMeasure> measure = nullptr;
            std::fstream result_file;
            if (h0) {
                measure = std::make_shared<PersistenceModularityMeasure>(grafo, result_file);
            } else {
                measure = std::make_shared<PersistenceMeasure>(grafo);
            }
            auto cm = measure->globalValue(membership, nullptr);
            auto result_cm_r = Proteggi(Rf_allocVector(REALSXP, (int) cm.size()), conta_protect);

            double value = 0.0;
            for (std::uint_fast64_t colonna = 0; colonna < cm.size(); ++colonna) {
                double s = cm.at(colonna).first + cm.at(colonna).second;
                REAL(result_cm_r)[colonna] = s;
                value += (std::isnan(s) ? 0.0: s);
            }

            auto result_value_r = Proteggi(Rf_allocVector(REALSXP, 1), conta_protect);
            REAL(result_value_r)[0] = value;

            result_r = Proteggi(Rf_allocVector(VECSXP, 2), conta_protect);
            auto result_r_names = Proteggi(Rf_allocVector(VECSXP, 2), conta_protect);
            SET_VECTOR_ELT(result_r, 0, result_value_r);
            SET_VECTOR_ELT(result_r_names, 0, Rf_mkChar("value"));
            SET_VECTOR_ELT(result_r, 1, result_cm_r);
            SET_VECTOR_ELT(result_r_names, 1, Rf_mkChar("clusters_value"));
            Rf_setAttrib(result_r, R_NamesSymbol, result_r_names);

        } catch(std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch(...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        if (conta_protect > 0) UNPROTECT(conta_protect);
        return result_r;
    }
}

extern "C" {
    SEXP localPersistence(SEXP vertex_r, SEXP edge_list_r, SEXP weights_r, SEXP cluster_r, SEXP h0_r) {
        int conta_protect = 0;
        SEXP result_r = R_NilValue;
        try {
            long vertex_size = Rf_length(vertex_r);
            std::vector<std::string> position_to_vertici(vertex_size);
            std::map<std::string, std::uint_fast64_t> vertici_to_position;
            for (long k = 0; k < vertex_size; ++k) {
              std::string e = R_CHAR(STRING_ELT(vertex_r, k));
              position_to_vertici.at((std::uint_fast64_t) k) = e;
              vertici_to_position[e] = (std::uint_fast64_t) k;
            }
            std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>> edge_list;
            LoadEdgeListR(edge_list_r, weights_r, edge_list, vertici_to_position);

            auto grafo = UGraph(position_to_vertici.size(), edge_list);

            std::set<std::uint_fast64_t> cluster;
            std::uint_fast64_t cluster_r_size = Rf_length(cluster_r);
            std::set<std::uint_fast64_t> used;

            for (std::uint_fast64_t p = 0; p < cluster_r_size; ++p) {
                auto value = INTEGER(cluster_r)[p];
                if (value == 1) {
                    cluster.insert(p);
                }
            }

            bool h0 = LOGICAL(h0_r)[0];
            std::shared_ptr<CommunityMeasure> measure = nullptr;
            std::fstream result_file;
            if (h0) {
                measure = std::make_shared<PersistenceModularityMeasure>(grafo, result_file);
            } else {
                measure = std::make_shared<PersistenceMeasure>(grafo);
            }

            auto cm = measure->localValue(cluster, nullptr);
            result_r = Proteggi(Rf_allocVector(REALSXP, (int) 1), conta_protect);
            REAL(result_r)[0] = cm.first - cm.second;
        } catch(std::exception &ex) {
            forward_exception_to_r(ex.what());
        } catch(...) {
            forward_exception_to_r("c++ exception (unknown reason)");
        }
        if (conta_protect > 0) UNPROTECT(conta_protect);
        return result_r;
    }
}

extern "C" {
  SEXP clusterMilano(SEXP vertex_r, SEXP edge_list_r, SEXP weights_r, SEXP membership_r, SEXP seed_r) {
      int conta_protect = 0;
      SEXP result_r = R_NilValue;
      try {
          long vertex_size = Rf_length(vertex_r);
          std::vector<std::string> position_to_vertici(vertex_size);
          std::map<std::string, std::uint_fast64_t> vertici_to_position;
          for (long k = 0; k < vertex_size; ++k) {
            std::string e = R_CHAR(STRING_ELT(vertex_r, k));
            position_to_vertici.at((std::uint_fast64_t) k) = e;
            vertici_to_position[e] = (std::uint_fast64_t) k;
          }
          std::list<std::tuple<std::uint_fast64_t, std::uint_fast64_t, double>> edge_list;
          LoadEdgeListR(edge_list_r, weights_r, edge_list, vertici_to_position);

          std::vector<std::uint_fast64_t> membership(vertex_size);

          for (long p = 0; p < vertex_size; ++p) {
              std::uint_fast64_t value = (std::uint_fast64_t) INTEGER(membership_r)[p];
              membership.at(p) = value;
          }

          long seed_r_length = Rf_length(seed_r);
          std::uint_fast64_t seed_rnd = 0;

          if (seed_r_length > 0) {
              seed_rnd = (std::uint_fast64_t) INTEGER(seed_r)[0];
          } else {
              seed_rnd = RandomUni::GENERATORE_SEED_RANDOM.RndNextInt(0, std::numeric_limits<std::uint_fast64_t>::max());
          }

          auto rnd = std::make_shared<RandomUni>(seed_rnd);


          auto grafo = UGraph(position_to_vertici.size(), edge_list);

          auto start_partition = std::make_shared<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>>();

          for (std::uint_fast64_t k = 0; k < membership.size(); ++k) {
              auto cluster_id = membership.at(k);
              auto cluster_set = start_partition->insert({cluster_id, std::make_shared<std::set<std::uint_fast64_t>>()});
              cluster_set.first->second->insert(k);
          }

          {
              std::map<std::uint_fast64_t, std::uint_fast64_t> cluster_map;
              for (std::uint_fast64_t k = 0; k < membership.size(); ++k) {
                  auto membership_value = membership.at(k);
                  auto find_iter = cluster_map.find(membership_value);
                  std::uint_fast64_t cluster_id = cluster_map.size();
                  if (find_iter == cluster_map.end()) {
                      cluster_map.insert({membership_value, cluster_map.size()});
                  } else {
                      cluster_id = find_iter->second;
                  }
                  auto cluster_set = start_partition->insert({cluster_id, std::make_shared<std::set<std::uint_fast64_t>>()});
                  cluster_set.first->second->insert(k);
              }
          }

          std::fstream result_file;
          auto measure = std::make_shared<PersistenceModularityMeasure>(grafo, result_file);

          std::shared_ptr<std::vector<double>> pi_value = nullptr;

          std::shared_ptr<std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>> partizione_best;
          std::pair<double, std::vector<double>> fo_best = std::make_pair(0.0, std::vector<double>());
          std::fstream debug_file;

          louvainMethod(grafo,
                        *measure,
                        start_partition,
                        pi_value,
                        partizione_best,
                        fo_best,
                        debug_file);

          auto membership_r = Proteggi(Rf_allocVector(INTSXP, position_to_vertici.size()), conta_protect);
          int p = 1;

          for (auto& c: *partizione_best) {
              for (auto& v: *c.second){
                  INTEGER(membership_r)[v] = p;
              }
              ++p;
          }

          auto best_value_r = Proteggi(Rf_allocVector(REALSXP, 1), conta_protect);
          REAL(best_value_r)[0] = fo_best.first;

          auto seed_r = Proteggi(Rf_allocVector(STRSXP, 1), conta_protect);
          SET_STRING_ELT(seed_r, 0, Rf_mkChar(std::to_string(rnd->Seed()).c_str()));

          result_r = Proteggi(Rf_allocVector(VECSXP, 3), conta_protect);
          auto result_r_names = Proteggi(Rf_allocVector(VECSXP, 3), conta_protect);
          SET_VECTOR_ELT(result_r, 0, membership_r);
          SET_VECTOR_ELT(result_r_names, 0, Rf_mkChar("membership"));
          SET_VECTOR_ELT(result_r, 1, best_value_r);
          SET_VECTOR_ELT(result_r_names, 1, Rf_mkChar("measure"));
          SET_VECTOR_ELT(result_r, 2, seed_r);
          SET_VECTOR_ELT(result_r_names, 2, Rf_mkChar("seed"));
          Rf_setAttrib(result_r, R_NamesSymbol, result_r_names);

      } catch(std::exception &ex) {
          forward_exception_to_r(ex.what());
      } catch(...) {
          forward_exception_to_r("c++ exception (unknown reason)");
      }
      if (conta_protect > 0) UNPROTECT(conta_protect);
      return result_r;
    }
}



extern "C" {
  static const R_CallMethodDef CallEntries[] = {
      {"_cluster_milano", (DL_FUNC) &clusterMilano, 5},
      {"_global_persistence", (DL_FUNC) &globalPersistence, 5},
      {"_local_persistence", (DL_FUNC) &localPersistence, 5},
      {NULL, NULL, 0}
  };

  void R_init_persistence(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
  }
}













