//
// Created by Alex Dyck and Leon Sievers ; last modified 1/18/18.
//

/**
 * @file tree.hpp
 *
 * @brief Definition of the OneTree class and general thoughts about a tree structure
 */

#ifndef BRANCHANDBOUNDTSP_TREE_HPP
#define BRANCHANDBOUNDTSP_TREE_HPP

#include<cstdlib>
#include <vector>
#include <stdexcept>
#include "graph.hpp"
#include "util.hpp"

namespace TSP {
using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;

/*
 * General Intentions
 * Given, that we are in the namespace TSP one might think about different trees that could be implemented.
 * Here in particular, we only implement the 1-Tree structure needed by Held-Karp Lower Bound Algorithm.
 * In fact, it's akin to a graph class, but with different constructor and a narrow use for further
 * programming. Formally, one would construct 1-trees such that they share inheritage with a tree and/or
 * a graph. Here, only the Node class is used.
 */
/**
 * @class OneTree
 * 1-tree implementation as a container for Nodes and meta information
 */
class OneTree {
 public:
  /**
   * Only constructor. Initializes the meta information and edge vector
   * @param size size of the tree in terms of Nodes
   */
  OneTree(size_t size) : _size(size), _edges(0) {
      _nodes.push_back(Node()); // root pushback. Here, root is always 0
      for (NodeId node = 1; node < size; node++)
          _nodes.push_back(Node());
      num_nodes = _nodes.size();
      num_edges = 0;
  }

  /**
   * Adding an edge by their NodeIds. Requires the to_Edge functionality from @headerfile util.hpp
   * @param i first node
   * @param j second node
   */
  void add_edge(NodeId i, NodeId j) {
      if (i >= num_nodes || j >= num_nodes)
          throw std::runtime_error("Index out of range while adding edge to 1-tree");

      _edges.push_back(to_EdgeId(i, j, _size));
      _nodes.at(i).add_neighbor(j);
      _nodes.at(j).add_neighbor(i);
      num_edges++;
  }
  //getter functions
  const size_type &get_num_edges() const {
      return num_edges;
  }
  const std::vector<EdgeId> &get_edges() const {
      return _edges;
  }

  const std::vector<Node> &get_nodes() const {
      return _nodes;
  }

  const Node &get_node(NodeId id) const {
      return _nodes.at(id);
  }

 private:
  //names are self-explainatory
  size_type _size;

  std::vector<Node> _nodes;
  std::vector<EdgeId> _edges;

  size_type num_edges;
  size_type num_nodes;
};

}

#endif //BRANCHANDBOUNDTSP_TREE_HPP
