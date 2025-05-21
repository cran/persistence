//
//  matriciCommunity.cpp
//  communityDetection
//
//  Created by Alessandro Avellone on 16/02/24.
//

#include <iostream>
#include <cstdint>
#include <memory>
#include <utility>
#include <numeric>
#include <vector>
#include <map>
#include "eccezioni.h"
#include "random.h"
#include "matriciCommunity.h"

// ***********************************************
// ***********************************************
// ***********************************************

//std::uint_fast64_t VettoreBool::NUMBER_OF_BITS = sizeof(std::uint_fast64_t) * 8;



// ***********************************************
// ***********************************************
// ***********************************************

MatriceBool::MatriceBool(std::uint_fast64_t lato) {
    this->lato = lato;
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        dati.push_back(std::make_unique<std::vector<bool>>(lato, false));
    }
}

MatriceBool::~MatriceBool() {
}

std::string MatriceBool::to_string() const {
    std::string DELIMETER = ";";
    std::string result = "";
    
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        std::string r = "";
        for (std::uint_fast64_t h = 0; h < lato; ++h) {
            auto v = dati.at(k)->at(h);
            if (h < lato - 1) {
                r += std::to_string(v) + ";";
            } else {
                r += std::to_string(v);
            }
        }
        result += r + "\n";
    }
    
    return result;
}

std::string MatriceBool::to_string(std::uint_fast64_t row) const {
    std::string r = "";
    for (std::uint_fast64_t h = 0; h < lato; ++h) {
        auto v = dati.at(row)->at(h);
        if (h < lato - 1) {
            r += std::to_string(v) + ";";
        } else {
            r += std::to_string(v);
        }
    }
    return r;
}

void MatriceBool::set(std::uint_fast64_t row, std::uint_fast64_t col, bool val) {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceBool error: set";
        throw_line(err_str);
    }
    dati.at(row)->at(col) = val;
}

void MatriceBool::setDiagonal(bool val) {
    for (std::uint_fast64_t h = 0; h < this->lato; ++h) {
        dati.at(h)->at(h) = val;
    }
}

bool MatriceBool::at(std::uint_fast64_t row, std::uint_fast64_t col) const {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceBool error: at";
        throw_line(err_str);
    }
    return dati.at(row)->at(col);
}

std::uint_fast64_t MatriceBool::Lato() const {
    return lato;
}

void MatriceBool::assign(bool v) {
    for (std::uint_fast64_t h = 0; h < this->lato; ++h) {
        std::fill(dati.at(h)->begin(), dati.at(h)->end(), v);

    }
}

void MatriceBool::copy(MatriceBool& copia_da) {
    for (std::uint_fast64_t h = 0; h < this->lato; ++h) {
        std::copy(copia_da.dati.at(h)->begin(), copia_da.dati.at(h)->end(), dati.at(h)->begin());
    }
}

void MatriceBool::makeSymmetric(){
    for (std::uint_fast64_t riga = 0; riga < lato; riga++) {
        for (std::uint_fast64_t colonna = riga + 1; colonna < lato; colonna++) {
            set(colonna, riga, at(riga, colonna));
        }
    }
}

std::uint_fast64_t MatriceBool::sum() {
    std::uint_fast64_t result = 0;
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        result += std::accumulate(this->dati.at(k)->begin(), this->dati.at(k)->end(), 0);
    }
    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

MatriceUInt::MatriceUInt(std::uint_fast64_t lato) {
    this->lato = lato;
    dati.resize(lato, std::vector<std::uint_fast64_t>(lato));
}

MatriceUInt::~MatriceUInt() {
}

std::string MatriceUInt::to_string() const {
    std::string DELIMETER = ";";
    std::string r = "";
    
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        for (std::uint_fast64_t h = 0; h < lato; ++h) {
            r += std::to_string(this->at(k, h)) + (h + 1 <  lato ? DELIMETER: "");
        }
        r += "\n";
    }
    
    return r;
}

std::string MatriceUInt::to_string(std::uint_fast64_t row) const {
    const char DELIMETER = ';';
    std::string r = "";
    
    for (std::uint_fast64_t h = 0; h < lato; ++h) {
        r += std::to_string(this->at(row, h)) + (h + 1 <  lato ? DELIMETER : ' ');
    }
    return r;
}

void MatriceUInt::set(std::uint_fast64_t row, std::uint_fast64_t col, std::uint_fast64_t val) {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceUInt error: set";
        throw_line(err_str);
    }
    dati.at(row).at(col) = val;
}

void MatriceUInt::setDiagonal(std::uint_fast64_t val) {
    for (std::uint_fast64_t h = 0; h < this->lato; ++h) {
        dati.at(h).at(h) = val;
    }
}

void MatriceUInt::add(std::uint_fast64_t row, std::uint_fast64_t col, std::uint_fast64_t val) {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceUInt error: set";
        throw_line(err_str);
    }
    dati.at(row).at(col) += val;
}

std::uint_fast64_t MatriceUInt::at(std::uint_fast64_t row, std::uint_fast64_t col) const {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceUInt error: at";
        throw_line(err_str);
    }
    return dati.at(row).at(col);
}

std::uint_fast64_t MatriceUInt::Lato() const {
    return lato;
}

void MatriceUInt::assign(std::uint_fast64_t v) {
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        dati.at(k).assign(lato, v);
    }
}

void MatriceUInt::copy(MatriceUInt& copia_da) {
    this->dati = copia_da.dati;
}

void MatriceUInt::copy(MatriceBool& copia_da) {
    for (std::uint_fast64_t riga = 0; riga < this->lato; ++riga) {
        for (std::uint_fast64_t colonna = 0; colonna < this->lato; ++colonna) {
            set(riga, colonna, copia_da.at(riga, colonna));
        }
    }
}

void MatriceUInt::makeSymmetric(){
    for (std::uint_fast64_t riga = 0; riga < lato; riga++) {
        for (std::uint_fast64_t colonna = riga + 1; colonna < lato; colonna++) {
            set(colonna, riga, at(riga, colonna));
        }
    }
}

std::uint_fast64_t MatriceUInt::sum() {
    std::uint_fast64_t result = 0;
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        result += std::accumulate(this->dati.at(k).begin(), this->dati.at(k).end(), 0);
    }
    
    return result;
}

// ***********************************************
// ***********************************************
// ***********************************************

MatriceDouble::MatriceDouble(std::uint_fast64_t lato) {
    this->lato = lato;
    dati.resize(lato, std::vector<double>(lato));
}

MatriceDouble::~MatriceDouble() {
}

std::string MatriceDouble::to_string() const {
    std::string DELIMETER = ";";
    std::string r = "";
    
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        for (std::uint_fast64_t h = 0; h < lato; ++h) {
            r += std::to_string(this->at(k, h)) + (h + 1 <  lato ? DELIMETER: "");
        }
        r += "\n";
    }
    
    return r;
}

std::string MatriceDouble::to_string(std::uint_fast64_t row) const {
    const char DELIMETER = ';';
    std::string r = "";
    
    for (std::uint_fast64_t h = 0; h < lato; ++h) {
        r += std::to_string(this->at(row, h)) + (h + 1 <  lato ? DELIMETER : ' ');
    }
    return r;
}

void MatriceDouble::set(std::uint_fast64_t row, std::uint_fast64_t col, double val) {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceUInt error: set";
        throw_line(err_str);
    }
    dati.at(row).at(col) = val;
}


double MatriceDouble::at(std::uint_fast64_t row, std::uint_fast64_t col) const {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceUInt error: at";
        throw_line(err_str);
    }
    return dati.at(row).at(col);
}

void MatriceDouble::add(std::uint_fast64_t row, std::uint_fast64_t col, double val) {
    if (row >= lato || col >= lato) {
        std::string err_str = "MatriceUInt error: set";
        throw_line(err_str);
    }
    dati.at(row).at(col) += val;
}

std::uint_fast64_t MatriceDouble::Lato() const {
    return lato;
}

void MatriceDouble::assign(double v) {
    for (std::uint_fast64_t k = 0; k < lato; ++k) {
        dati.at(k).assign(lato, v);
    }
}

void MatriceDouble::getDiagonal(std::vector<double>& diag) {
    for (std::uint_fast64_t i = 0; i < lato; i++){
        diag.at(i) = at(i, i);
    }
}

void MatriceDouble::copy(MatriceDouble& copia_da) {
    this->dati = copia_da.dati;
}

// ***********************************************
// ***********************************************
// ***********************************************

BoolMatrixEncoding::BoolMatrixEncoding(std::uint_fast64_t vsize) {
    vector_size = vsize;
    encoding_block_size = sizeof(std::uint_fast64_t) * 8;
    encoding_vector_size = vector_size / encoding_block_size + 1;

    encoding = std::make_shared<EncType>();
}

BoolMatrixEncoding::~BoolMatrixEncoding() {
}

std::pair<bool, std::uint_fast64_t> BoolMatrixEncoding::add(std::shared_ptr<std::vector<bool>> vett) {
    if (dati.size() > 0 && vett->size() != dati.begin()->second->size()) {
        std::string err_str = "MatriceBool error: at";
    }
    auto e = encode(*vett);
    auto r = check(*e);
    if (r == std::numeric_limits<std::uint_fast64_t>::max()) {
        auto key = (dati.size() == 0 ? 0 : std::prev(dati.end())->first + 1);
        dati[key] = vett;
        dati_code[key] = e;
        addEncoding(key, e);
        
        auto vett_set = std::make_shared<std::set<std::uint_fast64_t>>();
        for (std::uint_fast64_t h = 0; h < vett->size(); ++h) {
            if (vett->at(h)) {
                vett_set->insert(h);
            }
        }
        dati_set[key] = vett_set;
        
        return std::make_pair(true, key);
    }
    return std::make_pair(false, r);
}

std::pair<bool, std::uint_fast64_t> BoolMatrixEncoding::addNoCheck(std::shared_ptr<std::vector<bool>> vett) {
    if (dati.size() > 0 && vett->size() != dati.begin()->second->size()) {
        std::string err_str = "MatriceBool error: at";
    }
    if (dati.size() == 0) {
        encoding_vector_size = vett->size() / encoding_block_size + 1;
    }
    auto e = encode(*vett);
    
    auto key = (dati.size() == 0 ? 0 : std::prev(dati.end())->first + 1);
    dati[key] = vett;
    dati_code[key] = e;
    addEncoding(key, e);
    
    auto vett_set = std::make_shared<std::set<std::uint_fast64_t>>();
    for (std::uint_fast64_t h = 0; h < vett->size(); ++h) {
        if (vett->at(h)) {
            vett_set->insert(h);
        }
    }
    dati_set[key] = vett_set;
    
    return std::make_pair(true, key);
}

std::uint_fast64_t BoolMatrixEncoding::size() const {
    return dati.size();
}


std::pair<std::uint_fast64_t, std::uint_fast64_t> BoolMatrixEncoding::shape() const {
    return std::make_pair(vector_size, dati.size());
}

std::tuple<std::shared_ptr<std::vector<bool>>, std::shared_ptr<std::set<std::uint_fast64_t>>, std::shared_ptr<std::vector<std::uint_fast64_t>>>
BoolMatrixEncoding::at(std::uint_fast64_t position) const {
    return std::make_tuple(dati.at(position), dati_set.at(position), dati_code.at(position));
}


std::shared_ptr<std::vector<std::uint_fast64_t>> BoolMatrixEncoding::encode(std::vector<bool>& v) {
    auto risultato = std::make_shared<std::vector<std::uint_fast64_t>>(encoding_vector_size);
    std::uint_fast64_t k = 0;
    for (std::uint_fast64_t h = 0; h < encoding_vector_size && k < v.size(); ++h) {
        std::uint_fast64_t valore = 0;
        for (std::uint_fast64_t p = 0; p < encoding_block_size && k < v.size(); ++k, ++p) {
            valore = (valore << 1) | v.at(k);
        }
        risultato->at(h) = valore;
    }
    return risultato;
}

std::uint_fast64_t BoolMatrixEncoding::check(std::vector<std::uint_fast64_t>& e) const {
    std::shared_ptr<EncType> a = encoding;
    
    if (a->d.find(e.at(0)) == a->d.end()) {
        return std::numeric_limits<std::uint_fast64_t>::max();
    }
    
    std::set<std::uint_fast64_t> result(a->d.at(e.at(0)).first->begin(), a->d.at(e.at(0)).first->end());
    a = a->d[e.at(0)].second;

    for (std::uint_fast64_t p = 1; p < e.size(); ++p) {
        if (a->d.find(e.at(p)) == a->d.end()) {
            return std::numeric_limits<std::uint_fast64_t>::max();
        }
        auto current_Columns = a->d[e.at(p)].first;
        for (auto it = result.begin(); it != result.end(); ) {
            if (current_Columns->find(*it) == current_Columns->end()) {
                it = result.erase(it);
            }
            else {
                ++it;
            }
        }
        a = a->d[e.at(p)].second;
    }
    if (result.size() == 1) {
        return *result.begin();
    } else {
        return std::numeric_limits<std::uint_fast64_t>::max();
    }
}

void BoolMatrixEncoding::addEncoding(std::uint_fast64_t key, std::shared_ptr<std::vector<std::uint_fast64_t>> e)  {
    std::shared_ptr<EncType> a = encoding;
    for (std::uint_fast64_t p = 0; p < e->size(); ++p) {
        if (a->d.size() == 0) {
            a->d[e->at(p)] = std::make_pair(std::make_shared<std::set<std::uint_fast64_t>>(), std::make_shared<EncType>());
        } else {
            if (a->d.find(e->at(p)) == a->d.end()) {
                a->d[e->at(p)] = std::make_pair(std::make_shared<std::set<std::uint_fast64_t>>(), a->d.begin()->second.second);
            }
        }
        a->d[e->at(p)].first->insert(key);
        a = a->d[e->at(p)].second;
    }
    return;
}

/*
 void BoolMatrixEncoding::addEncoding(std::shared_ptr<std::vector<std::uint_fast64_t>> e)  {
     std::shared_ptr<EncType> a = encoding;
     for (std::uint_fast64_t p = 0; p < e->size(); ++p) {
         if (a->d.size() == 0) {
             a->d[e->at(p)] = std::make_pair(std::set<std::shared_ptr<std::vector<std::uint_fast64_t>>>(), std::make_shared<EncType>());
         }
         if (a->d.find(e->at(p)) == a->d.end()) {
             a->d[e->at(p)] = std::make_pair(std::set<std::shared_ptr<std::vector<std::uint_fast64_t>>>(), a->d.begin()->second.second);
         }
         a->d[e->at(p)].first.insert(e);
         a = a->d[e->at(p)].second;
     }
 }
 
 */


