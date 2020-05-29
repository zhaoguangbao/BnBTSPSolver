#pragma once

#include <string>
#include <list>
#include <fstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <algorithm>

template <typename data_type>
class Parser
{
public:
    enum convertType{
        tp_uint_t,
        tp_double_t
    };
public:
    explicit Parser(const std::string& filename);
    ~Parser(){file.close();}
    const std::vector<std::vector<data_type>>& getCitiesMatrix() const;

    void loadCitiesList(convertType tp);
    void loadCitiesMatrix(convertType tp);
    void loadAtspCitiesMatrix(convertType tp);
    void loadLowerDiagonalRow(convertType tp);
private:
    void skipToCoordinates();
    //Strips all colons from a given string
    void stripColons(std::string &x);
    std::vector<std::vector<data_type>> convertToMatrix(const std::list<std::tuple<int, data_type, data_type>>& citiesList);
    std::ifstream file;
    std::vector<std::vector<data_type>> citiesMatrix;
};

//implementation

template <typename data_type>
Parser<data_type>::Parser(const std::string& filename) : file(std::ifstream(filename))
{
    if (!file.is_open())
    {
        throw std::runtime_error("Can't open file");
    }
}

template <typename data_type>
const std::vector<std::vector<data_type> >& Parser<data_type>::getCitiesMatrix() const
{
    return citiesMatrix;
}

template <typename data_type>
void Parser<data_type>::loadCitiesList(convertType tp)
{
    std::list<std::tuple<int, double, double>> citiesList;
    skipToCoordinates();

    int cityNumber;
    while(file.good())
    {
        bool tp_intflag= tp==tp_uint_t ? true:false;
        if(tp_intflag)
        {
            uint32_t firstCoordinate, secondCoordinate;
            file >> cityNumber >> firstCoordinate >> secondCoordinate;
            citiesList.emplace_back(cityNumber, firstCoordinate, secondCoordinate);
        }else{
            double firstCoordinate, secondCoordinate;
            file >> cityNumber >> firstCoordinate >> secondCoordinate;
            citiesList.emplace_back(cityNumber, firstCoordinate, secondCoordinate);
        }
    }
    citiesMatrix = convertToMatrix(citiesList);
}

template <typename data_type>
void Parser<data_type>::loadCitiesMatrix(convertType tp)
{
    bool tp_intflag= tp==tp_uint_t ? true:false;
    decltype(citiesMatrix.size()) dimension = 0;
    std::string line;
    file >> dimension;
    citiesMatrix.resize(dimension);
    for (decltype(citiesMatrix.size()) i = 0; i < dimension; ++i)
    {
        for (decltype(citiesMatrix.size()) j = 0; j < dimension; ++j)
        {
            file >> line;
            if(tp_intflag)
                citiesMatrix[i].push_back(std::stoi(line));
            else {
                citiesMatrix[i].push_back(std::stod(line));
            }
        }
    }
}

template <typename data_type>
void Parser<data_type>::loadAtspCitiesMatrix(convertType tp)
{
    bool tp_intflag= tp==tp_uint_t ? true:false;
    decltype(citiesMatrix.size()) dimension = 0;
    std::string line = "", option = "";
    bool scan = false;
    do {
        getline(file, line);
        stripColons(line);
        std::stringstream strstr;
        strstr << line;
        strstr >> option;

        if (option == "DIMENSION") {
            strstr >> dimension;
        }

        if (option == "EDGE_WEIGHT_SECTION") {
            scan = true;
        }
    } while (!scan && file.good());

    if (!scan)
        throw std::runtime_error("File not in right format");

    citiesMatrix.resize(dimension);
    for (decltype(citiesMatrix.size()) i = 0; i < dimension; ++i)
    {
        for (decltype(citiesMatrix.size()) j = 0; j < dimension; ++j)
        {
            file >> line;
            if(tp_intflag)
                citiesMatrix[i].push_back(std::stoi(line));
            else {
                citiesMatrix[i].push_back(std::stod(line));
            }
        }
    }
}

template <typename data_type>
void Parser<data_type>::loadLowerDiagonalRow(convertType tp)
{
    bool tp_intflag= tp==tp_uint_t ? true:false;
    decltype(citiesMatrix.size()) dimension = 0;
    std::string line;
    file >> dimension;
    citiesMatrix.resize(dimension);
    for (size_t i = 0; i < dimension; ++i)
    {
        for (size_t j = 0; j < i + 1; ++j)
        {
            file >> line;
            if(tp_intflag)
                citiesMatrix[i].push_back(std::stoi(line));
            else {
                citiesMatrix[i].push_back(std::stod(line));
            }
        }
    }
    for (size_t i = 0; i < dimension; ++i)
    {
        for (size_t j = i + 1; j < dimension; ++j)
        {
            citiesMatrix[i].push_back(citiesMatrix[j][i]);
        }
    }
}

//private
template <typename data_type>
void Parser<data_type>::skipToCoordinates()
{
    std::string line;
    while (file.good())
    {
        file >> line;
        if (line == "NODE_COORD_SECTION")
        {
            break;
        }
    }
}

template <typename data_type>
void Parser<data_type>::stripColons(std::string &x) {
    auto it = std::remove_if(std::begin(x), std::end(x), [](char c) { return (c == ':'); });
    x.erase(it, std::end(x));
}


template <typename data_type>
std::vector<std::vector<data_type>> Parser<data_type>::convertToMatrix(const std::list<std::tuple<int, data_type, data_type>>& citiesList)
{
    std::vector<std::vector<double>> matrix(citiesList.size(), std::vector<double>(citiesList.size()));
    double x, y;
    auto it1 = citiesList.begin();
    for (size_t i = 0; i < citiesList.size(); ++i, ++it1)
    {
        auto it2 = citiesList.begin();
        for (size_t k = 0; k < citiesList.size(); ++k, ++it2)
        {
            x = std::get<1>(*it1) - std::get<1>(*it2);
            y = std::get<2>(*it1) - std::get<2>(*it2);
            matrix[i][k] = static_cast<int>(std::lround(std::sqrt((x * x + y * y))));
        }
    }
    return matrix;
}

