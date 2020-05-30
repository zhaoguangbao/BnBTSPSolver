//
// Modified by Alex Dyck and Leon Sievers ; last modified 1/18/18.
//

/**
 * @file graph.cpp
 *
 * @brief one small function was not defined inline
 */
#include "../header/graph.hpp"


namespace TSP {
/////////////////////////////////////////////
//! \c Node definitions
/////////////////////////////////////////////

    void Node::add_neighbor(NodeId const id) {
        _neighbors.push_back(id);
    }


} // namespace TSP
