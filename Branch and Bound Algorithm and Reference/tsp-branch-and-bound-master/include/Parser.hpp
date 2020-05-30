#pragma once

#include <string>
#include <list>
#include <fstream>
#include <vector>
#include <tuple>

class Parser
{
public:
    explicit Parser(const std::string& filename);
    std::vector<std::vector<int>> getCitiesMatrix() const;

    void loadCitiesList();
    void loadCitiesMatrix();
    void loadLowerDiagonalRow();
private:
    void skipToCoordinates();
    std::vector<std::vector<int>> convertToMatrix(const std::list<std::tuple<int, int, int>>& citiesList);

    std::ifstream file;
    std::vector<std::vector<int>> citiesMatrix;
};
