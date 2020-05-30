#pragma once

#include <vector>
#include <limits>

struct Instance
{
public:
    Instance() = delete;
    explicit Instance(std::vector<std::vector<int>> costMatrix);
    const std::vector<int>& operator[](int pos) const { return matrix[pos]; }

    std::vector<std::vector<int>> matrix;
    std::vector<int> bestPath;
    decltype(matrix.size()) size{0};
    int minCost = std::numeric_limits<int>::max();
};
