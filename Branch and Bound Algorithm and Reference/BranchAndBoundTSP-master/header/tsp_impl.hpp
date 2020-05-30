//
// Created by Alex Dyck and Leon Sievers ; last modified 1/18/18.
//

/**
 * @file tsp_impl.hpp
 *
 * @brief Implementation details for the tsp.hpp. One might ask, why we are using a .hpp. Well,
 * templates have to be beforehand which is achieved by putting them into the header. This also holds for
 * their implementation. There is one other way, but this is the way conforming the current ISO
 */
#ifndef BRANCHANDBOUNDTSP_TSP_IMPL_HPP
#define BRANCHANDBOUNDTSP_TSP_IMPL_HPP

#include <cassert>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <sstream>
#include <queue>
#include <numeric>
#include "tree.hpp"

namespace TSP {

/**
 * computes a minimum-1-tree for a given BranchingNode
 * @tparam coord_type
 * @tparam dist_type
 * @param tree space to save the optimal tree
 * @param lambda if we are in the root node, we'll save our holy lambda here, else it's just the root lambda
 * @param tsp The TSP Instance
 * @param BNode the correspnding BranchingNode
 */
template<class coord_type, class dist_type>
void compute_minimal_1_tree(TSP::OneTree &tree,
                            const std::vector<double> &lambda,
                            const TSP::Instance<coord_type, dist_type> &tsp,
                            const TSP::BranchingNode<coord_type, dist_type> &BNode) {

    //compute modified weights c_\lambda and set the weight of required edges
    // to -inf and for forbidden edges to +inf
    std::vector<dist_type> mod_weights(tsp.num_edges(), 0);

    for (size_t edge = 0; edge < tsp.weights().size(); edge++) {
        TSP::NodeId v = 0, w = 0;
        to_NodeId(edge, v, w, tsp.size());
        mod_weights.at(edge) = tsp.weight(edge) + lambda[v] + lambda[w];
    }

    for (const auto &el : BNode.get_forbidden())
        mod_weights.at(el) = std::numeric_limits<dist_type>::max();
    for (const auto &el : BNode.get_required())
        mod_weights.at(el) = -1; //std::numeric_limits<dist_type>::min();


    // computing a MST on {2,..,n} by PRIM MST Algorithm
    typedef std::pair<dist_type, int> Pair;
    //instantiate our priority_queue properly
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair> > pq;

    int src = 1; // Start at the first node != 0
    TSP::size_type n = tsp.size();
    // First, make all nodes unreachable
    std::vector<double> key(n, ::std::numeric_limits<double>::max() / 2.);

    // parent will give access to the second node in an edge for the MST
    std::vector<int> parent(n, -1);

    // included vertices vector
    std::vector<bool> MST_contained(n, false);
    //start with the source....
    pq.push(std::make_pair(0, src));
    key[src] = 0;

    while (!pq.empty()) {
        TSP::NodeId u = pq.top().second;
        pq.pop();

        MST_contained[u] = true;  // Include vertex in MST

        for (TSP::NodeId i = 1; i < n; i++) {
            if (i != u) {
                dist_type weight = mod_weights.at(to_EdgeId(u, i, n));
                if (MST_contained[i] == false && key[i] > weight) {
                    // Updating key of i
                    key[i] = weight;
                    pq.push(std::make_pair(key[i], static_cast<int>(i)));
                    parent[i] = u;
                }
            }
        }
    }
    //Done with MST computation, add them to our tree
    for (TSP::NodeId k = 2; k < n; k++) {
        tree.add_edge(k, static_cast<NodeId >(parent[k]));
    }

    //seek for smallest two edges incident to 0 ..
    TSP::NodeId smallest = 1;
    for (TSP::NodeId k = 2; k < n; k++) {
        if (mod_weights.at(to_EdgeId(0, k, n)) < mod_weights.at(to_EdgeId(0, smallest, n))) {
            smallest = k;
        }
    }

    TSP::NodeId smallest1 = 1;
    if (smallest == 1) smallest1 = 2;
    for (TSP::NodeId k = 2; k < n; k++) {
        if (k != smallest) {
            if (mod_weights.at(to_EdgeId(0, k, n)) < mod_weights.at(to_EdgeId(0, smallest1, n))) {
                smallest1 = k;
            }
        }
    }
    // ..add them
    tree.add_edge(0, smallest);
    tree.add_edge(0, smallest1);
}


/**
 * computes the Held-Karp lower bound
 * @tparam coord_type
 * @tparam dist_type
 * @param tsp The TSP Instance
 * @param lambda lambda set by root BranchingNode (or where to set for not BranchingNode)
 * @param tree container for tree computation
 * @param bn current BranchingNode
 * @param root true, if we are in the root of our B'n'B tree
 * @return
 */
template<class coord_type, class dist_type>
dist_type Held_Karp(const TSP::Instance<coord_type, dist_type> &tsp,
                    std::vector<double> &lambda,
                    TSP::OneTree &tree,
                    const TSP::BranchingNode<coord_type, dist_type> &bn,
                    bool root = false) {
    // Initialization
    TSP::size_type n = tsp.size();
    std::vector<dist_type> sol_vector;
    std::vector<double> lambda_max(lambda.size(), 0), lambda_tmp(lambda);
    TSP::OneTree tree_max(tree), tree_tmp(tree);
    double t_0 = 0., del_0 = 0., deldel = 0.;
    size_t max_el = 0;
    size_t N = std::ceil(n / 4.) + 5;
    if (root) {
        N = std::ceil(n * n / 50.) + n + 15;
    }
    // First tree computation to obtain t_0, del_0 , deldel
    compute_minimal_1_tree<coord_type, dist_type>(tree, lambda_tmp, tsp, bn);
    if (root) {
        dist_type sum = 0;
        for (const auto &el : tree.get_edges())
            sum += tsp.weight(el);
        t_0 = sum / (2. * n);
    } else {
        t_0 = 0;
        for (TSP::NodeId i = 0; i < n; i++) {
            t_0 += fabs(lambda.at(i));
        }
        t_0 *= 1. / (2. * n);
    }

    deldel = t_0 / (N * N - N);
    del_0 = 3. * t_0 / (2. * N);

    for (size_t i = 0; i < N; i++) {
        //Computing the sum we later on want to maximize over
        dist_type sum = 0, sum2 = 0;
        for (const auto &el : tree.get_edges())
            sum += tsp.weight(el);
        for (size_t node = 0; node < n; node++)
            sum2 += (tree.get_node(node).degree() - 2.) * lambda_tmp[node];
        sol_vector.push_back(sum + sum2); //we save all, not necessary, but nice for understanding

        if (i == 0) { // the first iteration is slightly different..
            tree_max = tree;
            lambda_max = lambda_tmp;

            for (size_t j = 0; j < lambda_tmp.size(); j++) {
                lambda_tmp[j] += t_0 * (tree.get_node(j).degree() - 2.);
            }
            t_0 = t_0 - del_0;
            del_0 = del_0 - deldel;
            tree_tmp = tree;
        }

        if (i > 0) { // ..then this one
            if (sol_vector[max_el] < sol_vector[i]) {
                lambda_max = lambda_tmp;
                max_el = i;
                tree_max = tree;
            }
            for (size_t j = 0; j < lambda_tmp.size(); j++) {
                lambda_tmp[j] +=
                    t_0 * (0.6 * (tree.get_node(j).degree() - 2.) + 0.4 * (tree_tmp.get_node(j).degree() - 2.));
            }
            t_0 = t_0 - del_0;
            del_0 = del_0 - deldel;
            tree_tmp = tree;
        }
        tree = TSP::OneTree(n);
        compute_minimal_1_tree<coord_type, dist_type>(tree, lambda_tmp, tsp, bn);
    }
    if (root) { //Setting the holy lambda
        lambda = lambda_max;
    }
    tree = tree_max;
    // Multiplying by 1. - EPS whereas EPS is a Macro defined to 10e-7 since we do not want to
    // obtain a lower bound larger than the optimum solution. This could occur due to
    // floating point computations .
    return std::ceil((1. - EPS) * (*std::max_element(sol_vector.begin(), sol_vector.end())));
}

// ---------------------------------------------------------------------------------
// ---------------    TSP::Instance section ----------------------------------------
// ---------------------------------------------------------------------------------
template<class coord_type, class dist_type>
Instance<coord_type, dist_type>::Instance(const std::string &filename) : _length(0){
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("File " + filename + " could not be opened");

    std::string line = "", option = "";
    bool scan = false;
    do {
        getline(file, line);
        stripColons(line);
        std::stringstream strstr;
        strstr << line;
        strstr >> option;

        if (option == "DIMENSION") {
            strstr >> this->dimension;
        }

        if (option == "NODE_COORD_SECTION") {
            scan = true;
        }
    } while (!scan && file.good());

    if (!scan)
        throw std::runtime_error("File not in right format");

    std::vector<coord_type> x, y;
    x.reserve(dimension), y.reserve(dimension);
    coord_type coord_x = std::numeric_limits<coord_type>::max(), coord_y = std::numeric_limits<coord_type>::max();
    while (file.good()) {
        getline(file, line);
        std::stringstream strstr;
        strstr << line;
        strstr >> option;
        if (option != "EOF") {
            try {
                strstr >> coord_x >> coord_y;
                _nodes.push_back(std::stoi(option) - 1), x.push_back(coord_x), y.push_back(coord_y);
            }
            catch (int e) {
                std::cerr << "An exception occurred while reading the file. Exception Nr. " << e << '\n';
            }
        } else break;
    }
    file.close();
    for (size_t i = 0; i < dimension; i++)
        for (size_t j = 0; j < dimension; j++)
            this->_weights.push_back(
                distance(x[i], y[i], x[j], y[j])
            );
    _tour = std::vector<NodeId>(dimension);
}

template<class coord_type, class dist_type>
void Instance<coord_type, dist_type>::compute_optimal_tour() {
    typedef BranchingNode<coord_type, dist_type> BNode;

    dist_type upperBound = std::numeric_limits<dist_type>::max();
    std::priority_queue<BNode,
                        std::vector<BNode>, std::greater<BNode> > Q;
    Q.push(BranchingNode<coord_type, dist_type>(*this)); // Adding empty node to Q

    while (!Q.empty()) {
        BNode current_BNode(Q.top());
        Q.pop();
        if (current_BNode.get_HK() >= upperBound)
            continue;
        else {
            if (current_BNode.tworegular()) {
                upperBound = current_BNode.get_HK();
                std::cerr << "Upper Bound " << upperBound << std::endl;
                _tour = current_BNode.get_tree().get_edges();
                continue;
            } else {
                size_type gl_i = 0, choice1 = std::numeric_limits<size_type>::max(),
                    choice2 = std::numeric_limits<size_type>::max();
                for (NodeId node = 1; node < current_BNode.get_tree().get_nodes().size(); node++) {
                    if (current_BNode.get_tree().get_node(node).degree() > 2) {
                        gl_i = node;
                        break;
                    }
                }

                assert(gl_i != 0);
                size_t counter = 0;
                for (const auto &el : current_BNode.get_tree().get_node(gl_i).neighbors()) {
                    if (!current_BNode.is_required(to_EdgeId(gl_i, el, size()))) {
                        assert(!current_BNode.is_forbidden(to_EdgeId(gl_i, el, size())));
                        if (counter == 0)
                            choice1 = el;
                        if (counter == 1)
                            choice2 = el;
                        counter++;
                        if (counter > 1)
                            break;
                    }
                }
                assert(choice1 < std::numeric_limits<NodeId>::max());
                assert(choice2 < std::numeric_limits<NodeId>::max());
                BNode q1(current_BNode, *this, to_EdgeId(gl_i, choice1, this->size())),
                    q2(current_BNode,
                       *this,
                       to_EdgeId(gl_i, choice1, this->size()),
                       to_EdgeId(gl_i, choice2, this->size()));
                Q.push(q1);
                Q.push(q2);

                if (!current_BNode.get_required_neighbors().at(gl_i).degree()) {
                    BNode q3(current_BNode, *this,
                             to_EdgeId(gl_i, choice1, this->size()),
                             to_EdgeId(gl_i, choice2, this->size()),
                             true);
                    Q.push(q3);
                }
            }
        }
    }
    std::cerr << "Optimal Length " << upperBound << std::endl;

    this->_length = upperBound;
}

template<class coord_type, class dist_type>
void Instance<coord_type, dist_type>::print_optimal_tour(const std::string &filename) {
    if (this->_tour.size() != this->size())
        throw std::runtime_error("No tour computed yet!");
    size_type n = this->size();
    std::ofstream file_to_print;
    file_to_print.open(filename, std::ios::out);

    file_to_print << "TYPE : TOUR" << std::endl;
    file_to_print << "DIMENSION : " << this->size() << std::endl;
    file_to_print << "TOUR_SECTION" << std::endl;
    std::vector<NodeId> path(n);
    std::vector<bool> visited(n, false);
    std::vector<std::vector<NodeId> > x(n, std::vector<NodeId>());

    for (const auto &el : _tour) {
        NodeId v = 0, w = 0;
        to_NodeId(el, v, w, n);
        x.at(v).push_back(w);
        x.at(w).push_back(v);
    }

    NodeId current = 0;
    file_to_print << 1 << std::endl;
    visited.at(current) = true;
    for (size_t i = 0; i < n - 1; i++) {
        NodeId neighbour = x.at(current).at(0);
        if (visited.at(neighbour))
            neighbour = x.at(current).at(1);
        file_to_print << neighbour + 1 << std::endl;
        current = neighbour;
        visited.at(current) = true;
    }

    file_to_print << "-1" << std::endl;
    file_to_print << "EOF" << std::endl;
}

// end class Instance section

// --------------------------------------------------------
// -------------------- BranchingNode section -------------
// --------------------------------------------------------

template<class coord_type, class dist_type>
bool BranchingNode<coord_type, dist_type>::operator>(const BranchingNode <coord_type, dist_type> &rhs) const {
    return this->get_HK() > rhs.get_HK();
}

template<class coord_type, class dist_type>
void BranchingNode<coord_type, dist_type>::forbid(NodeId idx, EdgeId e1, EdgeId e2) {
    for (NodeId k = 0; k < size; k++) {
        if (idx != k) {
            EdgeId edge = to_EdgeId(idx, k, size);
            if (edge != e1 && edge != e2) {
                add_forbidden(edge);
            }
        }
    }
}

template<class coord_type, class dist_type>
void BranchingNode<coord_type, dist_type>::admit(NodeId idx) {
    for (NodeId k = 0; k < size; k++) {
        if (idx != k) {
            EdgeId edge = to_EdgeId(idx, k, size);
            if (!is_forbidden(edge))
                add_required(edge);
        }
    }
}

template<class coord_type, class dist_type>
bool BranchingNode<coord_type, dist_type>::push_required(EdgeId e) {
    if (is_required(e))
        return false;
    assert(!is_forbidden(e));
    NodeId i = 0, j = 0;
    to_NodeId(e, i, j, size);
    required.push_back(e);
    required_neighbors.at(i).add_neighbor(j);
    required.push_back(reverse_edge(e, size));
    required_neighbors.at(j).add_neighbor(i);
    return true;
}

template<class coord_type, class dist_type>
void BranchingNode<coord_type, dist_type>::add_required(EdgeId e) {
    if (!push_required(e))
        return;
    NodeId i = 0, j = 0;
    to_NodeId(e, i, j, size);

    if (required_neighbors.at(i).degree() == 2)
        forbid(i, to_EdgeId(i, required_neighbors.at(i).neighbors().at(0), size),
               to_EdgeId(i, required_neighbors.at(i).neighbors().at(1), size));
    if (required_neighbors.at(j).degree() == 2)
        forbid(j, to_EdgeId(j, required_neighbors.at(j).neighbors().at(0), size),
               to_EdgeId(j, required_neighbors.at(j).neighbors().at(1), size));
}

template<class coord_type, class dist_type>
bool BranchingNode<coord_type, dist_type>::push_forbidden(EdgeId e) {
    if (is_forbidden(e))
        return false;
    assert(!is_required(e));
    NodeId i = 0, j = 0;
    to_NodeId(e, i, j, size);
    forbidden.push_back(e);
    forbidden_neighbors[i].add_neighbor(j);
    forbidden.push_back(reverse_edge(e, size));
    forbidden_neighbors[j].add_neighbor(i);
    return true;
}

template<class coord_type, class dist_type>
void BranchingNode<coord_type, dist_type>::add_forbidden(EdgeId e) {
    // hier pushen wir die edges entsprechend
    if (!push_forbidden(e))
        return;
    NodeId i = 0, j = 0;
    to_NodeId(e, i, j, size);

    if (forbidden_neighbors.at(i).degree() == size - 3) {
        admit(i);
    }
    if (forbidden_neighbors.at(j).degree() == size - 3) {
        admit(j);
    }

}

// end class BranchingNone section

} //end namespace TSP

#endif //BRANCHANDBOUNDTSP_TSP_IMPL_HPP
