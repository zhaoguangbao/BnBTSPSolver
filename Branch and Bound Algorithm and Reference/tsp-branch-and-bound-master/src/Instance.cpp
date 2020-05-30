#include "include/Instance.hpp"

Instance::Instance(std::vector<std::vector<int>> costMatrix) : matrix(std::move(costMatrix))
{
    size = matrix.size();

    for (auto i = 0; i < matrix.size(); ++i)
    {
        matrix[i][i] = std::numeric_limits<int>::max();
    }
}
