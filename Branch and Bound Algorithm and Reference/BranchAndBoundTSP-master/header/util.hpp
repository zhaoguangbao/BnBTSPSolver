//
// Created by Alex Dyck and Leon Sievers ; last modified 1/18/18.
//


#ifndef BRANCHANDBOUNDTSP_UTIL_HPP
#define BRANCHANDBOUNDTSP_UTIL_HPP
/**
 * @file util.cpp
 *
 * @brief some helper functions in the TSP Framework used
 */
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>



// Helper functions for accessing a two dimensional array, which is saved
// in one dimension in the context of graphs and adjacent matrices

using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;

/**
 * Compute the EdgeId from two NodeIds
 * @param i first Node
 * @param j second Node
 * @param N Number of nodes in the underlying graph/instance
 * @return the id of an edge (i,j) , iff edges are assigned properly
 */
EdgeId to_EdgeId(NodeId i, NodeId j, size_type N) {
    if (i == j)
        throw std::runtime_error("Loops are not contained in this instance");
    if (i > j)
        std::swap(i, j);
    return (i) * (N) + j;
}

/**
 * Compute the NodeIds corresponding to a given edge
 * @param e EdgeId we want we want to extract NodeIds from
 * @param i placeholder for first NodeId
 * @param j placeholder for second NodeId
 * @param N Number of nodes in the underlying graph/instance
 */
void to_NodeId(EdgeId e, NodeId &i, NodeId &j, size_type N) {
    j = e % (N);
    i = (e - j) / (N);
}

/**
 * Strips all colons from a given string
 * @param x String that needs to be stripped from Colons
 */
void stripColons(std::string &x) {
    auto it = std::remove_if(std::begin(x), std::end(x), [](char c) { return (c == ':'); });
    x.erase(it, std::end(x));
}

#endif //BRANCHANDBOUNDTSP_UTIL_HPP

