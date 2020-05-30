#include "include/Solver.hpp"

#include <algorithm>
#include <iostream>
#include <chrono>

namespace
{
constexpr auto intMax = std::numeric_limits<int>::max();
} // namespace

Solver::Solver(std::vector<std::vector<int>> costMatrix) : instance(std::move(costMatrix))
{
}

void Solver::run()
{
    std::vector<int> path;
    path.reserve(instance.size);
    auto begin = std::chrono::system_clock::now();

    path.push_back(0);
    Node root{path};
    branchAndBound(root);

    auto end = std::chrono::system_clock::now();
    elapsedTime = end - begin;
}

void Solver::printSolution()
{
    std::cout << std::endl;
    std::cout << "Cost: " << instance.minCost << std::endl;
    std::cout << std::endl;

    instance.bestPath.push_back(0);
    std::cout << "Path: ";

    for (const auto& e : instance.bestPath)
    {
        std::cout << e << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(elapsedTime).count() << "ms ";
    std::cout << std::chrono::duration_cast<std::chrono::seconds>(elapsedTime).count() << "s" << std::endl;
}

std::vector<bool> Solver::getVisited(const std::vector<int>& path) const
{
    std::vector<bool> visited(instance.size, false);
    for (auto i = 0; i < instance.size; ++i)
    {
        if (find(path.begin(), path.end(), i) != path.end())
        {
            visited[i] = true;
        }
    }
    return visited;
}

std::vector<int> Solver::getMinCosts(const std::vector<bool>& visited) const
{
    std::vector<int> values;
    for (auto i = 0; i < instance.size; ++i)
    {
        auto min = intMax;
        if (!visited[i])
        {
            for (auto j = 0; j < instance.size; ++j)
            {
                if (!visited[j])
                {
                    min = std::min(min, instance[i][j]);
                }
            }
            values.push_back(min);
        }
    }
    return values;
}

std::vector<int> Solver::getSortedCosts(const std::vector<int>& path) const
{
    std::vector<bool> visited = getVisited(path);
    std::vector<int> values = getMinCosts(visited);
    std::sort(values.begin(), values.end());

    return values;
}

int Solver::getLowestReturnCost(const std::vector<int>& path) const
{
    auto min = intMax;
    for (auto i = 0; i < instance.size; ++i)
    {
        min = std::min(min, instance[path.front()][i]);
    }
    return min;
}

int Solver::getLowerBound(const std::vector<int>& path, int cost) const
{
    auto lowerBound = cost;
    auto minCosts = getSortedCosts(path);
    for (const auto minCost : minCosts)
    {
        lowerBound += minCost;
        if (lowerBound > instance.minCost)
        {
            return lowerBound;
        }
    }
    return lowerBound + getLowestReturnCost(path);
}

void Solver::updateBestSolution(const Node& node)
{
    if (node.cost < instance.minCost)
    {
        instance.minCost = node.cost;
        instance.bestPath = node.currentPath;
    }
}

void Solver::createNode(const Node& node, std::vector<Node>& nodes, int i, int cost) const
{
    nodes.emplace_back(node.currentPath);
    nodes.back().currentPath.push_back(i);
    nodes.back().cost = cost;
}

std::vector<Node> Solver::createPromisingNodes(const Node& node) const
{
    std::vector<Node> nodes;
    nodes.reserve(instance.size - 1);

    for (auto i = 1; i < instance.size; ++i)
    {
        if (find(node.currentPath.begin(), node.currentPath.end(), i) == node.currentPath.end())
        {
            auto cost = node.cost + instance[node.currentPath.back()][i];
            if (cost < instance.minCost)
            {
                createNode(node, nodes, i, cost);
            }
        }
    }

    return nodes;
}

void Solver::branchAndBound(Node& node)
{
    if (node.currentPath.size() == instance.size)
    {
        node.cost += instance[node.currentPath.back()][node.currentPath.front()];
        updateBestSolution(node);
        return;
    }
    if (getLowerBound(node.currentPath, node.cost) <= instance.minCost)
    {
        std::vector<Node> nodes = createPromisingNodes(node);
        for (auto& next : nodes)
        {
            branchAndBound(next);
        }
    }

}
