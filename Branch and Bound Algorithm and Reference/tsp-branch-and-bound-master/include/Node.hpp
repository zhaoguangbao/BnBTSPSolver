#pragma once

#include <vector>

struct Node
{
    Node() = delete;
    explicit Node(std::vector<int> path) : currentPath(std::move(path)) {}

    std::vector<int> currentPath;
    int cost{0};
};
