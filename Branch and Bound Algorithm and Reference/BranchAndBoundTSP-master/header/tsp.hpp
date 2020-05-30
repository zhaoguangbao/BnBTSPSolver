//
// Created by Alex Dyck and Leon Sievers ; last modified 1/18/18.
//

/**
 * @file tsp.hpp
 *
 * @brief header file for the main TSP framework including the TSP Instance and the BranchingNode templates
 */
#ifndef BRANCHANDBOUNDTSP_TSP_HPP
#define BRANCHANDBOUNDTSP_TSP_HPP


#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <climits>
#include <vector>
#include <sstream>
#include <queue>
#include <numeric>
#include <utility>
#include <cassert>
#include "util.hpp"
#include "tree.hpp"

#define EPS 10e-7


/**
 * @namespace TSP defines elements, methods and classes within the Travelling Salesman Context
 */
namespace TSP {

using size_type = std::size_t;
using NodeId = size_type;
using EdgeId = size_type;

template<class coord_type, class dist_type>
class BranchingNode;

template<class coord_type, class dist_type>
class Instance;

/**
 * @class Instance holds the edge weights and reads the instance from a file in TSPLIB format
 * @tparam coord_type Container in which the Coordinates are given. Assumably double
 * @tparam dist_type Container in which the distances  are given. Assumably double
 */
template<class coord_type, class dist_type>
class Instance {
 public:
  /**
   * Constructor of @class Instance which takes the filename as an argument.
   * File has to be in TSPLIB format
   * @param filename
   */
  Instance(const std::string &filename);

  /**
   *  Distance function in the TSP Instance. Could be also made private, since it's only there
   *  at init point.
   * @param x1
   * @param y1
   * @param x2
   * @param y2
   * @return  rounded \f$ \sqrt{(x_1 - x_2)^2  + (y_1 - y_2)^2}\f$
   */
  dist_type distance(coord_type x1, coord_type y1, coord_type x2, coord_type y2) {
      return std::lround(std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
  }

  /**
   * Computes and optimal tour on this Instance and saves it as EdgeIds in _tour
   */
  void compute_optimal_tour();

  /**
   * Output the optimal tour into a file by TSPLIB rules
   * @param filename
   */
  void print_optimal_tour(const std::string &filename);

  //Getter functions

  size_type size() const {
      return dimension;
  }

  size_type num_edges() const {
      return _weights.size();
  }

  dist_type weight(EdgeId id) const {
      return _weights.at(id);
  }

  const std::vector<dist_type> &weights() const {
      return _weights;
  }
  const dist_type & length() {
      return _length;
  }
 private:
  std::vector<NodeId> _nodes;
  std::vector<dist_type> _weights;
  size_type dimension;
  std::vector<NodeId> _tour;
  dist_type _length;
};

/**
 * @class BranchingNode represents a Node in our branch and bound tree. It contains the
 * lower bound for itself and saves temporary information as tree, lambda, required, forbidden.
 * @comment Yeah, for the different branching nodes and their 'to append' edges, we implemented
 * separate constructors. It could be achieved by one, but this was easier, as we are
 * lacking some structures necessary for this. - Alex
 * @tparam coord_type Container in which the Coordinates are given. Assumably double
 * @tparam dist_type Container in which the distances  are given. Assumably double
 */
template<class coord_type, class dist_type>
class BranchingNode {
 public:
  /**
   * First Constructor: Constructs a BranchingNode without any forbidden or required edges,
   * i.e. the root of our B'n'B tree.
   * @param tsp The TSP Instance
   */
  BranchingNode(const Instance<coord_type, dist_type> &tsp
  ) : size(tsp.size()), required(), required_neighbors(size), forbidden(),
      forbidden_neighbors(size), lambda(size, 0), tree(size) {
      HK = Held_Karp(tsp, this->lambda, this->tree, *this, true);
  }
  /**
   * Second Constructor: Constructs a BranchingNode with \f$ F := F \cup e_1 \f$
   * @param BNode predecessor BranchingNode
   * @param tsp The TSP Instance
   * @param e1 additional edge for forbidden edges
   */
  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1
  ) : size(tsp.size()),
      required(BNode.get_required()),
      required_neighbors(BNode.required_neighbors),
      forbidden(BNode.get_forbidden()),
      forbidden_neighbors(BNode.forbidden_neighbors),
      lambda(BNode.get_lambda()),
      tree(tsp.size()) {

      add_forbidden(e1);
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }
  /**
   * Third Constructor: Constructs a BranchingNode with \f$ R := R \cup {e_1} F := F \cup {e_2} \f$
   * @param BNode predecessor BranchingNode
   * @param tsp The TSP Instance
   * @param e1 additional edge for required edges
   * @param e2 additional edge for forbidden edges
   */
  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1,
                EdgeId e2
  ) : size(tsp.size()),
      required(BNode.get_required()),
      required_neighbors(BNode.required_neighbors),
      forbidden(BNode.get_forbidden()),
      forbidden_neighbors(BNode.forbidden_neighbors),
      lambda(BNode.get_lambda()),
      tree(tsp.size()) {
      add_required(e1);
      add_forbidden(e2);
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }
  /**
   * Fourth Constructor: Constructs a BranchingNode with \f$ R := R \cup {e_1 } \cup {e_2}\f$
   * @param BNode predecessor BranchingNode
   * @param tsp The TSP Instance
   * @param e1 additional edge for required edges
   * @param e2 additional edge for forbidden edges
   * @param both_req bool to distinguish third from forth constructor
   */
  BranchingNode(const BranchingNode<coord_type, dist_type> &BNode,
                const Instance<coord_type, dist_type> &tsp,
                EdgeId e1,
                EdgeId e2,
                bool both_req
  ) : size(tsp.size()),
      required(BNode.get_required()),
      required_neighbors(BNode.required_neighbors),
      forbidden(BNode.get_forbidden()),
      forbidden_neighbors(BNode.forbidden_neighbors),
      lambda(BNode.get_lambda()),
      tree(tsp.size()) {

      if (both_req) {
          add_required(e1);
          add_required(e2);
      }
      HK = Held_Karp(tsp, this->lambda, this->tree, *this);
  }
  /**
   * Overloading operator > and comparing lowerbounds of two BranchingNodes.
   * Used by priority_queue
   * @param rhs
   * @return true if >
   */
  bool operator>(const BranchingNode<coord_type, dist_type> &rhs) const;

  /**
   * For an EdgeId e corresponding to an edge (i,j) it returns the edge
   * corresponding to the edge (j,i)
   * @param e \f$ =(i,j)\f$
   * @param n size on the tsp instance
   * @return \f$(j,i) = e' \f$
   */
  EdgeId reverse_edge(EdgeId e, size_type n) const {
      NodeId i = 0, j = 0;
      to_NodeId(e, i, j, n);
      return j * n + i;
  }

  /**
   *
   * @param id EdgeId id
   * @return true, if the undirected edge {i,j} is required, else false
   */
  bool is_required(EdgeId id) const {
      return (std::find(required.begin(), required.end(), id) != required.end()
          && std::find(required.begin(), required.end(), reverse_edge(id, size)) != required.end());
  }
   /**
   *
   * @param id EdgeId id
   * @return true, if the undirected edge {i,j} is forbidden, else false
   */
  bool is_forbidden(EdgeId id) const {
      return (std::find(forbidden.begin(), forbidden.end(), id) != forbidden.end()
          && std::find(forbidden.begin(), forbidden.end(), reverse_edge(id, size)) != forbidden.end());
  }

  /**
   * forbids all edges incident to idx except e1,e2
   * @param idx
   * @param e1
   * @param e2
   */
  void forbid(NodeId idx, EdgeId e1, EdgeId e2) ;

  /**
   * adds non-forbidden edges of incident edges of idx to required
   * @param idx
   */
  void admit(NodeId idx) ;

  /**
   * pushes both, e and reverse edge of e to required
   * @param e
   * @return false, if it was already in required
   */
  bool push_required(EdgeId e);

  /**
   * external function to add an edge e to required
   * @param e
   */
  void add_required(EdgeId e) ;

  /**
   * pushes both, e and reverse edge of e to forbidden
   * @param e
   * @return false, if it was already forbidden
   */
  bool push_forbidden(EdgeId e);
  /**
   * external function to add an edge e to forbidden
   * @param e
   */
  void add_forbidden(EdgeId e);


  //Getter functions
  const std::vector<EdgeId> &get_required() const {
      return required;
  }
  const std::vector<EdgeId> &get_forbidden() const {
      return forbidden;
  }

  const std::vector<dist_type> &get_lambda() const {
      return lambda;
  }
  std::vector<dist_type> &get_lambda() {
      return lambda;
  }

  const OneTree &get_tree() const {
      return tree;
  }
  OneTree &get_tree() {
      return tree;
  }

  const dist_type get_HK() const {
      return this->HK;
  }

  const std::vector<Node> &get_required_neighbors() const {
      return required_neighbors;
  }
  /**
   * checks if the current tree is 2-regular
   * @return true, if T 2-regular
   */
  bool tworegular() {
      for (const auto &el : tree.get_nodes())
          if (el.degree() != 2)
              return false;
      return true;
  }

 private:
  size_type size;
  std::vector<EdgeId> required;
  std::vector<Node> required_neighbors;

  std::vector<EdgeId> forbidden;
  std::vector<Node> forbidden_neighbors;

  std::vector<double> lambda;
  OneTree tree;

  dist_type HK;
};
}

#include "tsp_impl.hpp"

#endif // BRANCHANDBOUNDTSP_TSP_HPP