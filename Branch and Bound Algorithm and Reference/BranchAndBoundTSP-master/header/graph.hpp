//
// Created by Alex Dyck and Leon Sievers ; last modified 1/18/18.
//

/**
   @file graph.hpp

   @brief This file was used as a graph class header, but someone did not like it and now its only use is to
   serve his owner with a Node class
**/
#ifndef GRAPH_HPP
#define GRAPH_HPP



#include <cstddef> // std::size_t
#include <iosfwd> // std::ostream fwd declare
#include <limits>
#include <vector>

namespace TSP // TSP
{

using size_type = std::size_t;
using NodeId = size_type;
using TSPlibId = size_type;

/**
   @class Node

   @brief A @c Node stores an array of neighbors (via their ids).

   @note The neighbors are not necessarily ordered, so searching for a specific neighbor takes O(degree)-time.
**/
class Node {
 public:
  typedef std::size_t size_type;

  /** @brief Create an isolated node (you can add neighbors later). **/
  Node() {};   //= default;

  /** @return The number of neighbors of this node. **/
  size_type degree() const;

  /** @return The array of ids of the neighbors of this node. **/
  std::vector<NodeId> const &neighbors() const;
  /**
     @brief Adds @c id to the list of neighbors of this node.
     @warning Does not check whether @c id is already in the list of neighbors (a repeated neighbor is legal, and
     models parallel edges).
     @warning Does not check whether @c id is the identity of the node itself (which would create a loop!).
  **/
  void add_neighbor(NodeId const id);
 private:
  friend class OneTree;
  friend class BranchingTree;

  std::vector<NodeId> _neighbors;
}; // class Node

//BEGIN: Inline section

inline
Node::size_type Node::degree() const {
    return neighbors().size();
}

inline
std::vector<NodeId> const &Node::neighbors() const {
    return _neighbors;
}





//END: Inline section



} // namespace ED

#endif /* GRAPH_HPP */
