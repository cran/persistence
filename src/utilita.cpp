//
//  utilita.cpp
//  POSet
//
//  Created by Alessandro Avellone on 09/04/2019.
//  Copyright © 2019 Alessandro Avellone. All rights reserved.
//

#include <fstream>
#include <filesystem>
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <mutex>
#include <map>
#include <list>
#include <set>
#include <algorithm>

//#include "types.h"
#include "utilita.h"
//#include "insieme.h"

extern const char DELIMETER;

// ***********************************************
// ***********************************************
// ***********************************************

std::string FindAndReplaceAll(const std::string& data, const std::string toSearch, const std::string replaceStr) {
    // Get the first occurrence
    std::string risultato = data;
    size_t pos = risultato.find(toSearch);

    // Repeat till end is reached
    while( pos != std::string::npos)
    {
        // Replace this occurrence of Sub String
        risultato.replace(pos, toSearch.size(), replaceStr);
        // Get the next occurrence from the current position
        pos = risultato.find(toSearch, pos + replaceStr.size());
    }
    return risultato;
}

// ***********************************************
// ***********************************************
// ***********************************************

std::string matrix_to_string(std::vector<std::vector<double>>& v, char DELIMETER) {
    std::string r = "";
    for (size_t k = 0; k <  v.size(); ++k) {
        bool first = true;
        for (size_t h = 0; h < v.at(k).size(); ++h) {
            if (first) {
                r += std::to_string(v.at(k).at(h));
                first = false;
            }
            else {
                r += DELIMETER + std::to_string(v.at(k).at(h));
            }
        }
        r += "\n";
    }

    return r;
}



//************************************
//************************************
//************************************

std::pair<std::shared_ptr<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>, std::shared_ptr<std::set<std::uint64_t>>> LoadPesiFromFile(std::string nomeFile, char dataSetColSpe) {
    std::ifstream f(nomeFile);
    if (!f.good()) {
        throw std::invalid_argument{"File not found: " + std::string(nomeFile)};
    }

    std::string line;
    std::ifstream fp(nomeFile);

    auto dati = std::make_shared<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>();
    auto nomi_vertici = std::make_shared<std::set<std::uint64_t>>();

    int riga = 0;
    while (getline(fp, line)) {
        Trim(line);
        if (line.length() == 0) continue;

        auto tokens = split(line, dataSetColSpe);

        for(std::uint64_t colonna = 0; colonna < tokens.size(); ++colonna) {
            std::string v = tokens.at(colonna);
            double peso = std::stod(v);
            dati->push_back(std::make_tuple(riga, colonna, peso));
        }
        nomi_vertici->insert(riga);
        ++riga;
    }
    fp.close();

    return std::make_pair(dati, nomi_vertici);
}


//************************************
//************************************
//************************************

std::pair<std::shared_ptr<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>, std::shared_ptr<std::map<std::string, std::uint64_t>>> LoadArchiFromFile(std::string nomeFile, char dataSetColSpe) {
    std::ifstream f(nomeFile);
    if (!f.good()) {
        throw std::invalid_argument{"File not found: " + std::string(nomeFile)};
    }

    std::string line;
    std::ifstream fp(nomeFile);

    auto dati = std::make_shared<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>();
    auto nomi_vertici = std::make_shared<std::map<std::string, std::uint64_t>>();
    std::uint64_t indice = 0;
    while (getline(fp, line)) {
        Trim(line);
        if (line.length() == 0) continue;

        auto tokens = split(line, dataSetColSpe);
        if (tokens.size() == 2 || tokens.size() == 3) {
            auto v1 = tokens.at(0);
            std::uint64_t v1_index = 0;
            if (nomi_vertici->find(v1) == nomi_vertici->end()) {
                v1_index = indice++;
                nomi_vertici->insert(std::make_pair(v1, v1_index));
            } else {
                v1_index = nomi_vertici->at(v1);
            }
            auto v2 = tokens.at(1);
            std::uint64_t v2_index = 0;
            if (nomi_vertici->find(v2) == nomi_vertici->end()) {
                v2_index = indice++;
                nomi_vertici->insert(std::make_pair(v2, v2_index));
            } else {
                v2_index = nomi_vertici->at(v2);
            }
            double peso = 1.0;
            if (tokens.size() == 3) {
                peso = std::stod(tokens.at(2));
            }
            dati->push_back(std::make_tuple(v1_index, v2_index, peso));
        }
    }
    fp.close();

    return std::make_pair(dati, nomi_vertici);
}

//************************************
//************************************
//************************************

std::pair<std::shared_ptr<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>, std::shared_ptr<std::map<std::string, std::uint64_t>>> LoadArchiFromFileV2(std::string nomeFile, char dataSetColSpe) {
    std::ifstream f(nomeFile);
    if (!f.good()) {
        throw std::invalid_argument{"File not found: " + std::string(nomeFile)};
    }
    // Open the stream to 'lock' the file.
    std::ifstream fp(nomeFile, std::ios::in | std::ios::binary);
    const auto file_size = std::filesystem::file_size(nomeFile);
    std::string file_str(file_size, '\0');
    fp.read(file_str.data(), file_size);
    fp.close();

    std::uint64_t p = file_str.find('\n');
    std::string line = "";
    if (p == std::string::npos) {
        line = file_str;
    } else {
        line = file_str.substr(0, p);
    }
    Trim(line);
    auto elements_in_row = split(line, dataSetColSpe).size();
    if (elements_in_row != 2 && elements_in_row != 3) {
        throw std::invalid_argument{"File line size: " + std::to_string(elements_in_row)};
    }

    auto dati = std::make_shared<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>();
    auto nomi_vertici = std::make_shared<std::map<std::string, std::uint64_t>>();
    std::uint64_t indice = 0;
    std::uint64_t start_pos = 0;
    bool continua = true;
    while (continua) {
        std::uint64_t end_pos = file_str.find('\n', start_pos);
        std::string line = "";
        if (end_pos == std::string::npos) {
            line = file_str.substr(start_pos);
            continua = false;
        } else {
            line = file_str.substr(start_pos, end_pos - start_pos);
        }
        Trim(line);
        start_pos = end_pos + 1;
        if (line.length() == 0) continue;

        auto tokens = split(line, dataSetColSpe);
        if (tokens.size() != elements_in_row) {
            continue;
        }

        auto v1 = tokens.at(0);
        std::uint64_t v1_index = 0;
        if (nomi_vertici->find(v1) == nomi_vertici->end()) {
            v1_index = indice++;
            nomi_vertici->insert(std::make_pair(v1, v1_index));
        } else {
            v1_index = nomi_vertici->at(v1);
        }
        auto v2 = tokens.at(1);
        std::uint64_t v2_index = 0;
        if (nomi_vertici->find(v2) == nomi_vertici->end()) {
            v2_index = indice++;
            nomi_vertici->insert(std::make_pair(v2, v2_index));
        } else {
            v2_index = nomi_vertici->at(v2);
        }
        double peso = std::numeric_limits<double>::signaling_NaN();
        if (elements_in_row == 3) {
            peso = std::stod(tokens.at(2));
        }
        dati->push_back(std::make_tuple(v1_index, v2_index, peso));

    }

    return std::make_pair(dati, nomi_vertici);
}

//************************************
//************************************
//************************************

void LTrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

//************************************
//************************************
//************************************

void RTrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

//************************************
//************************************
//************************************

void Trim(std::string &s) {
    LTrim(s);
    RTrim(s);
}

//************************************
//************************************
//************************************

std::vector<std::string>& split(const std::string &text, char sep) {
    static std::vector<std::string> tokens;

    tokens.clear();
    std::uint64_t start = 0, end = 0;
    while ((end = text.find(sep, start)) != std::string::npos) {
        tokens.push_back(text.substr(start, end - start));
        start = end + 1;
    }
    tokens.push_back(text.substr(start));
    return tokens;
}

