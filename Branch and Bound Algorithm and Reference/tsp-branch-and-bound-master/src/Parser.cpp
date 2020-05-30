#include "include/Parser.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>

Parser::Parser(const std::string& filename) : file(std::ifstream(filename))
{
    if (!file.is_open())
    {
        throw std::runtime_error("Can't open file");
    }
}

std::vector<std::vector<int> > Parser::getCitiesMatrix() const
{
    return citiesMatrix;
}

void Parser::loadCitiesList()
{
    std::list<std::tuple<int, int, int>> citiesList;
    skipToCoordinates();

    int cityNumber, firstCoordinate, secondCoordinate;
    for (;;)
    {
        file >> cityNumber >> firstCoordinate >> secondCoordinate;
        if (!file.good())
        {
            break;
        }
        citiesList.emplace_back(cityNumber, firstCoordinate, secondCoordinate);
    }
    citiesMatrix = convertToMatrix(citiesList);
}

void Parser::loadCitiesMatrix()
{
    decltype(citiesMatrix.size()) dimension = 0;
    std::string line;
    file >> dimension;
    citiesMatrix.resize(dimension);
    for (auto i = 0; i < dimension; ++i)
    {
        for (auto j = 0; j < dimension; ++j)
        {
            file >> line;
            citiesMatrix[i].push_back(std::stoi(line));
        }
    }
}

void Parser::loadLowerDiagonalRow()
{
    decltype(citiesMatrix.size()) dimension = 0;
    std::string line;
    file >> dimension;
    citiesMatrix.resize(dimension);
    for (auto i = 0; i < dimension; ++i)
    {
        for (auto j = 0; j < i + 1; ++j)
        {
            file >> line;
            citiesMatrix[i].push_back(std::stoi(line));
        }
    }
    for (auto i = 0; i < dimension; ++i)
    {
        for (auto j = i + 1; j < dimension; ++j)
        {
            citiesMatrix[i].push_back(citiesMatrix[j][i]);
        }
    }
}

void Parser::skipToCoordinates()
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

std::vector<std::vector<int>> Parser::convertToMatrix(const std::list<std::tuple<int, int, int>>& citiesList)
{
    std::vector<std::vector<int>> matrix(citiesList.size(), std::vector<int>(citiesList.size()));
    int x, y;
    auto it1 = citiesList.begin();
    for (auto i = 0; i < citiesList.size(); ++i, ++it1)
    {
        auto it2 = citiesList.begin();
        for (auto k = 0; k < citiesList.size(); ++k, ++it2)
        {
            x = std::get<1>(*it1) - std::get<1>(*it2);
            y = std::get<2>(*it1) - std::get<2>(*it2);
            matrix[i][k] = static_cast<int>(std::lround(std::sqrt((x * x + y * y))));
        }
    }
    return matrix;
}
