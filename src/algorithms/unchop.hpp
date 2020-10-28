#pragma once

#include <handlegraph/types.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>
#include <vector>

#include "simple_components.hpp"

namespace odgi {
namespace algorithms {

using namespace handlegraph;


/**
 * Unchop by gluing abutting handles with just a single edge between them and
 * compatible path steps together.
 * @param graph
 */
void unchop(handlegraph::MutablePathDeletableHandleGraph& graph);

/**
 * Unchop by gluing abutting handles with just a single edge between them and
 * compatible path steps together.
 * @param graph
 * @param show_info
 */

void unchop(handlegraph::MutablePathDeletableHandleGraph& graph,
            const uint64_t& nthreads,
            const bool& show_info);

//std::vector<std::deque<handle_t>> simple_components(PathHandleGraph* graph, int min_size = 1, false);

handle_t concat_nodes(handlegraph::MutablePathDeletableHandleGraph& graph, const std::vector<handle_t>& nodes);
//handle_t concat_nodes(handlegraph::MutablePathDeletableHandleGraph* graph, const std::vector<handle_t>& nodes);

}
}
