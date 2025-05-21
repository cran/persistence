#ifndef utilita_hpp
#define utilita_hpp

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <cstdint>
#include <memory>


std::string FindAndReplaceAll(const std::string& data, const std::string toSearch, const std::string replaceStr);
std::string matrix_to_string(std::vector<std::vector<double>>&, char DELIMETER = ';');

void LTrim(std::string &s);
void RTrim(std::string &s);
void Trim(std::string &s);
std::vector<std::string>& split(const std::string &text, char sep);
std::pair<std::shared_ptr<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>, std::shared_ptr<std::set<std::uint64_t>>> LoadPesiFromFile(std::string nomeFile, char dataSetColSpe);
std::pair<std::shared_ptr<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>, std::shared_ptr<std::map<std::string, std::uint64_t>>> LoadArchiFromFile(std::string nomeFile, char dataSetColSpe);
std::pair<std::shared_ptr<std::list<std::tuple<std::uint64_t, std::uint64_t, double>>>, std::shared_ptr<std::map<std::string, std::uint64_t>>> LoadArchiFromFileV2(std::string nomeFile, char dataSetColSpe);

//************************************
//************************************
//************************************


#endif /* utilita_hpp */
