#pragma once

#include <vector>
#include <chrono>

#include "Instance.hpp"
#include "Node.hpp"

class Solver
{
public:
    explicit Solver(std::vector<std::vector<int>> costMatrix);
    void run();
    void printSolution();
private:
    std::vector<bool> getVisited(const std::vector<int>& path) const;
    std::vector<int> getMinCosts(const std::vector<bool>& visited) const;
    std::vector<int> getSortedCosts(const std::vector<int>& path) const;
    int getLowestReturnCost(const std::vector<int>& path) const;
    int getLowerBound(const std::vector<int>& path, int cost) const;
    void updateBestSolution(const Node& node);
    void createNode(const Node& node, std::vector<Node>& nodes, int i, int cost) const;
    std::vector<Node> createPromisingNodes(const Node& node) const;
    void branchAndBound(Node& node);

    Instance instance;
    std::chrono::duration<double> elapsedTime;
};
