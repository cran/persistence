#ifndef MatriceUInt_h
#define MatriceUInt_h


#include <cstdint>
#include <memory>
#include <utility>
#include <numeric>
#include <vector>
#include <set>
#include <list>
#include <map>

class MatriceDouble;


// ***********************************************
// ***********************************************
// ***********************************************

template <class T>
std::string vector_to_string(std::vector<T>& v, char DELIMETER = ';') {
    std::string r = "";
    bool first = true;
    for (size_t h = 0; h < v.size(); ++h) {
        if (first) {
            r += std::to_string(v.at(h));
            first = false;
        }
        else {
            r += DELIMETER + std::to_string(v.at(h));
        }
    }
    return r;
}


// ***********************************************
// ***********************************************
// ***********************************************

class MatriceBool {
private:
    std::vector<std::unique_ptr<std::vector<bool>>> dati;
    std::uint_fast64_t lato;
public:
    MatriceBool(std::uint_fast64_t lato);
    ~MatriceBool();
    std::string to_string() const;
    std::string to_string(std::uint_fast64_t row) const;
    void set(std::uint_fast64_t row, std::uint_fast64_t col, bool val);
    void setDiagonal(bool val);
    bool at(std::uint_fast64_t row, std::uint_fast64_t col) const;
    std::uint_fast64_t Lato() const;
    void assign(bool);
    void copy(MatriceBool& copia_da);
    std::uint_fast64_t sum();
    void makeSymmetric();
};


// ***********************************************
// ***********************************************
// ***********************************************

class MatriceUInt {
private:
    std::vector<std::vector<std::uint_fast64_t>> dati;
    std::uint_fast64_t lato;
public:
    MatriceUInt(std::uint_fast64_t lato);
    ~MatriceUInt();
    std::string to_string() const;
    std::string to_string(std::uint_fast64_t row) const;
    void set(std::uint_fast64_t row, std::uint_fast64_t col, std::uint_fast64_t val);
    void add(std::uint_fast64_t row, std::uint_fast64_t col, std::uint_fast64_t val);
    void setDiagonal(std::uint_fast64_t val);
    std::uint_fast64_t at(std::uint_fast64_t row, std::uint_fast64_t col) const;
    std::uint_fast64_t Lato() const;
    void assign(std::uint_fast64_t);
    void copy(MatriceUInt& copia_da);
    void copy(MatriceBool& copia_da);
    std::uint_fast64_t sum();
    void makeSymmetric();
};

// ***********************************************
// ***********************************************
// ***********************************************

class MatriceDouble {
private:
    std::vector<std::vector<double>> dati;
    std::uint_fast64_t lato;
public:
    MatriceDouble(std::uint_fast64_t lato);
    ~MatriceDouble();
    std::string to_string() const;
    std::string to_string(std::uint_fast64_t row) const;
    void set(std::uint_fast64_t row, std::uint_fast64_t col, double val);
    double at(std::uint_fast64_t row, std::uint_fast64_t col) const;
    void add(std::uint_fast64_t row, std::uint_fast64_t col, double val);
    std::uint_fast64_t Lato() const;
    void assign(double v);
    void getDiagonal(std::vector<double>& diag);
    void copy(MatriceDouble& copia_da);
};

// ***********************************************
// ***********************************************
// ***********************************************
class BoolMatrixEncoding {
public:
    struct EncType {
        std::map<std::uint_fast64_t, std::pair<std::shared_ptr<std::set<std::uint_fast64_t>>, std::shared_ptr<EncType>>> d;
    };
private:
    using DATI = std::map<std::uint_fast64_t, std::shared_ptr<std::vector<bool>>>;
    using DATISET = std::map<std::uint_fast64_t, std::shared_ptr<std::set<std::uint_fast64_t>>>;
    using DATICODE = std::map<std::uint_fast64_t, std::shared_ptr<std::vector<std::uint_fast64_t>>>;

    DATI dati;
    DATISET dati_set;
    DATICODE dati_code;
    std::shared_ptr<EncType> encoding;
    std::uint_fast64_t vector_size;
    std::uint_fast64_t encoding_vector_size;
    std::uint_fast64_t encoding_block_size;
public:
    BoolMatrixEncoding(std::uint_fast64_t);
    ~BoolMatrixEncoding();
    std::pair<bool, std::uint_fast64_t> add(std::shared_ptr<std::vector<bool>> vett);
    std::pair<bool, std::uint_fast64_t> addNoCheck(std::shared_ptr<std::vector<bool>> vett);
    std::uint_fast64_t size() const;
    std::pair<std::uint_fast64_t, std::uint_fast64_t> shape() const;
    std::tuple<std::shared_ptr<std::vector<bool>>, std::shared_ptr<std::set<std::uint_fast64_t>>, std::shared_ptr<std::vector<std::uint_fast64_t>>>
    at(std::uint_fast64_t) const;
    
    struct Iterator {
        using value_type = std::tuple<DATI::iterator, DATISET::iterator, DATICODE::iterator>;

        Iterator(DATI::iterator itd, DATISET::iterator itds, DATICODE::iterator itdc)
            : itdati(itd), itdatiset(itds), itdaticode(itdc) {}
        Iterator(const Iterator& m)
            : itdati(m.itdati), itdatiset(m.itdatiset), itdaticode(m.itdaticode) {}
        
        
        value_type operator*() const {
            return std::make_tuple(itdati, itdatiset, itdaticode);
        }

        // Prefix increment
        Iterator& operator++() {
            ++itdati;
            ++itdatiset;
            ++itdaticode;
            return *this;
        }

        // Postfix increment
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

        friend bool operator== (const Iterator& a, const Iterator& b) {
            return (a.itdati == b.itdati &&
                    a.itdatiset == b.itdatiset &&
                    a.itdaticode == b.itdaticode);
        };
        friend bool operator!= (const Iterator& a, const Iterator& b) {
            return !(a.itdati == b.itdati && a.itdatiset == b.itdatiset && a.itdaticode == b.itdaticode);
        };
    private:
        DATI::iterator itdati;
        DATISET::iterator itdatiset;
        DATICODE::iterator itdaticode;
    };
    
    typedef Iterator iterator;
    
    iterator begin() {
        return iterator(dati.begin(), dati_set.begin(), dati_code.begin());
    }
    
    iterator end() {
        return iterator(dati.end(), dati_set.end(), dati_code.end());
    }
    
private:
    std::shared_ptr<std::vector<std::uint_fast64_t>> encode(std::vector<bool>& v);
    std::uint_fast64_t check(std::vector<std::uint_fast64_t>& e) const;
    void addEncoding(std::uint_fast64_t, std::shared_ptr<std::vector<std::uint_fast64_t>> e);
};

#endif /* MatriceUIntaria_h */
