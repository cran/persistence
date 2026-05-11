#pragma once

#include <string>
#include <cstdint>

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>
#include <R_ext/Parse.h>

// ***********************************************
// RAII wrapper per PROTECT/UNPROTECT automatico
// ***********************************************

class RProtectGuard {
private:
    int count_ = 0;
    
public:
    RProtectGuard() = default;
    
    // Non copiabile, non movibile
    RProtectGuard(const RProtectGuard&) = delete;
    RProtectGuard& operator=(const RProtectGuard&) = delete;
    RProtectGuard(RProtectGuard&&) = delete;
    RProtectGuard& operator=(RProtectGuard&&) = delete;
    
    ~RProtectGuard() {
        if (count_ > 0) {
            UNPROTECT(count_);
        }
    }
    
    SEXP protect(SEXP s) {
        ++count_;
        return PROTECT(s);
    }
    
    void release() {
        if (count_ > 0) {
            UNPROTECT(count_);
            count_ = 0;
        }
    }
    
    int count() const { return count_; }
};

// ***********************************************
// Helper per creare oggetti R comuni
// ***********************************************

namespace RHelper {
    inline SEXP createRealVector(R_xlen_t size, RProtectGuard& guard) {
        return guard.protect(Rf_allocVector(REALSXP, size));
    }
    
    inline SEXP createIntVector(R_xlen_t size, RProtectGuard& guard) {
        return guard.protect(Rf_allocVector(INTSXP, size));
    }
    
    inline SEXP createStrVector(R_xlen_t size, RProtectGuard& guard) {
        return guard.protect(Rf_allocVector(STRSXP, size));
    }
    
    inline SEXP createNamedList(int size, RProtectGuard& guard) {
        SEXP list = guard.protect(Rf_allocVector(VECSXP, size));
        SEXP names = guard.protect(Rf_allocVector(STRSXP, size));
        Rf_setAttrib(list, R_NamesSymbol, names);
        return list;
    }
    
    inline void setListElement(SEXP list, int index, SEXP value, const char* name) {
        SET_VECTOR_ELT(list, index, value);
        SEXP names = Rf_getAttrib(list, R_NamesSymbol);
        SET_STRING_ELT(names, index, Rf_mkChar(name));
    }
    
    inline SEXP createScalarReal(double value, RProtectGuard& guard) {
        SEXP result = guard.protect(Rf_allocVector(REALSXP, 1));
        REAL(result)[0] = value;
        return result;
    }
    
    inline SEXP createScalarString(const std::string& value, RProtectGuard& guard) {
        SEXP result = guard.protect(Rf_allocVector(STRSXP, 1));
        SET_STRING_ELT(result, 0, Rf_mkChar(value.c_str()));
        return result;
    }
}

// ***********************************************
// Funzioni utility
// ***********************************************

SEXP getFromEnv(SEXP env, const std::string& name);
SEXP getListElement(SEXP list, const char* str);

// ***********************************************
// Interfaccia C per R
// ***********************************************

extern "C" {
    SEXP cluster_milano(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
    SEXP global_persistence(SEXP, SEXP, SEXP, SEXP, SEXP);
    SEXP local_persistence(SEXP, SEXP, SEXP, SEXP, SEXP);
}

