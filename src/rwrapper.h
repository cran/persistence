#ifndef poset_wrapper_hpp
#define poset_wrapper_hpp

#include <string>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include <R_ext/Parse.h>

// ***********************************************
// ***********************************************
// ***********************************************


inline SEXP Proteggi(SEXP s, int& conta) {
    ++conta;
    return PROTECT(s);
}


// ***********************************************
// ***********************************************
// ***********************************************

SEXP getFromEnv(SEXP env, const std::string& name);
SEXP getListElement(SEXP list, const char *str);

// ***********************************************
// ***********************************************
// ***********************************************

extern "C" {
    SEXP clusterMilano(SEXP, SEXP, SEXP);
    SEXP globalPersistence(SEXP, SEXP, SEXP, SEXP);
    SEXP localPersistence(SEXP, SEXP, SEXP, SEXP);
}

// ***********************************************
// ***********************************************
// ***********************************************



#endif /* poset_wrapper_hpp */
