//
//  graph.cpp
//  

#include "graph.hpp"

namespace dankgraph {

graph_t::graph_t(void) { }
graph_t::~graph_t(void) { }

/// Look up the handle for the node with the given ID in the given orientation
handle_t graph_t::get_handle(const id_t& node_id, bool is_reverse) const {
    return handle_helper::pack(graph_id_wt.select(0, node_id), is_reverse);
}
    
/// Get the ID from a handle
id_t graph_t::get_id(const handle_t& handle) const {
    return graph_id_wt.at(handle_helper::unpack_number(handle));
}
    
/// Get the orientation of a handle
bool graph_t::get_is_reverse(const handle_t& handle) const {
    return handle_helper::unpack_bit(handle);
}
    
/// Invert the orientation of a handle (potentially without getting its ID)
handle_t graph_t::flip(const handle_t& handle) const {
    return handle_helper::toggle_bit(handle);
}
    
/// Get the length of a node
size_t graph_t::get_length(const handle_t& handle) const {
    uint64_t offset = handle_helper::unpack_number(handle);
    return seq_wt.select(offset+1, 0) - seq_wt.select(offset, 0);
}
    
/// Get the sequence of a node, presented in the handle's local forward orientation.
std::string graph_t::get_sequence(const handle_t& handle) const {
    std::string seq;
    uint64_t offset = handle_helper::unpack_number(handle);
    for (uint64_t i = seq_wt.select(offset, 0); ; ++i) {
        char c = seq_wt.at(i);
        if (c) { seq += c; } else break;
    }
    return seq;
}
    
/// Loop over all the handles to next/previous (right/left) nodes. Passes
/// them to a callback which returns false to stop iterating and true to
/// continue. Returns true if we finished and false if we stopped early.
bool graph_t::follow_edges(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    bool result = true;
    bool offset = handle_helper::unpack_number(handle);
    bool is_rev = handle_helper::unpack_bit(handle);
    // NB edges are stored in canonical orientation, forward to reverse prefered
    if (!go_left && !is_rev || go_left && is_rev) {
        uint64_t edges_begin = edge_fwd_wt.select(offset, 0)+1;
        for (uint64_t i = edges_begin; ; ++i) {
            id_t rank = edge_fwd_wt.at(i);
            if (rank==0) break; // end of record
            bool inv = edge_fwd_inv_bv.at(i);
            handle_t handle = handle_helper::pack(rank-1, (inv ? !is_rev : is_rev));
            result &= iteratee(handle);
            if (!result) break;
        }
    } else {
        assert(go_left && !is_rev || !go_left && is_rev);
        uint64_t edges_begin = edge_rev_wt.select(offset, 0)+1;
        for (uint64_t i = edges_begin; ; ++i) {
            id_t rank = edge_rev_wt.at(i);
            if (rank==0) break; // end of record
            bool inv = edge_rev_inv_bv.at(i);
            handle_t handle = handle_helper::pack(rank-1, (inv ? !is_rev : is_rev));
            result &= iteratee(handle);
            if (!result) break;
        }
    }
    return result;
}
    
/// Loop over all the nodes in the graph in their local forward
/// orientations, in their internal stored order. Stop if the iteratee
/// returns false. Can be told to run in parallel, in which case stopping
/// after a false return value is on a best-effort basis and iteration
/// order is not defined.
void graph_t::for_each_handle(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {
    if (parallel) {
        volatile bool flag=false;
#pragma omp parallel for
        for (uint64_t i = 0; i < graph_id_wt.size(); ++i) {
            if (flag) continue;
            bool result = iteratee(handle_helper::pack(i,false));
#pragma omp atomic
            flag &= result;
        }
    } else {
        for (uint64_t i = 0; i < graph_id_wt.size(); ++i) {
            if (!iteratee(handle_helper::pack(i,false))) break;
        }
    }
}
    
/// Return the number of nodes in the graph
/// TODO: can't be node_count because XG has a field named node_count.
size_t graph_t::node_size(void) const {
    graph_id_wt.size() - graph_id_wt.rank(graph_id_wt.size(), 0);
}
    
/// Return the smallest ID in the graph, or some smaller number if the
/// smallest ID is unavailable. Return value is unspecified if the graph is empty.
id_t graph_t::min_node_id(void) const {
    return _min_node_id;
}
    
/// Return the largest ID in the graph, or some larger number if the
/// largest ID is unavailable. Return value is unspecified if the graph is empty.
id_t graph_t::max_node_id(void) const {
    return _max_node_id;
}
    
////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////
    
/// Get the number of edges on the right (go_left = false) or left (go_left
/// = true) side of the given handle. The default implementation is O(n) in
/// the number of edges returned, but graph implementations that track this
/// information more efficiently can override this method.
size_t graph_t::get_degree(const handle_t& handle, bool go_left) const {
    uint64_t offset = handle_helper::unpack_number(handle);
    bool is_rev = handle_helper::unpack_bit(handle);
    if (!go_left && !is_rev || go_left && is_rev) {
        return edge_fwd_wt.select(offset+1, 0) - edge_fwd_wt.select(offset, 0);
    } else {
        return edge_rev_wt.select(offset+1, 0) - edge_rev_wt.select(offset, 0);
    }
}
    
////////////////////////////////////////////////////////////////////////////
// Concrete utility methods
////////////////////////////////////////////////////////////////////////////
    
/// Get a Protobuf Visit from a handle.
//Visit to_visit(const handle_t& handle) const;
    
/// Get the locally forward version of a handle
handle_t graph_t::forward(const handle_t& handle) const {
    handle_t handle_fwd = handle;
    if (handle_helper::unpack_bit(handle)) handle_helper::toggle_bit(handle_fwd);
    return handle_fwd;
}

/*
/// A pair of handles can be used as an edge. When so used, the handles have a
/// canonical order and orientation.
edge_t graph_t::edge_handle(const handle_t& left, const handle_t& right) const {
    make_pair(left, right);
}
    
/// Such a pair can be viewed from either inward end handle and produce the
/// outward handle you would arrive at.
handle_t graph_t::traverse_edge_handle(const edge_t& edge, const handle_t& left) const {
}
*/
    
/**
 * This is the interface for a handle graph that stores embedded paths.
 */
    
////////////////////////////////////////////////////////////////////////////
// Path handle interface that needs to be implemented
////////////////////////////////////////////////////////////////////////////
    
/// Determine if a path name exists and is legal to get a path handle for.
bool graph_t::has_path(const std::string& path_name) const {
    std::string query_s = "$" + path_name + "$";
    std::vector<uint64_t> query_v(query_s.begin(), query_s.end());
    return path_name_fmi.locate(query_v).size() > 0;
}
    
/// Look up the path handle for the given path name.
/// The path with that name must exist.
path_handle_t graph_t::get_path_handle(const std::string& path_name) const {
    std::string query_s = "$" + path_name + "$";
    std::vector<uint64_t> query_v(query_s.begin(), query_s.end());
    std::vector<uint64_t> occs = path_name_fmi.locate(query_v);
    assert(occs.size() == 1);
    uint64_t offset = occs.front();
    return as_path_handle(path_name_bv.rank1(offset));
}
    
/// Look up the name of a path from a handle to it
std::string graph_t::get_path_name(const path_handle_t& path_handle) const {
    return paths.at(as_integer(path_handle)).name;
}
    
/// Returns the number of node occurrences in the path
size_t graph_t::get_occurrence_count(const path_handle_t& path_handle) const {
    return paths.at(as_integer(path_handle)).occurrence_count();
}

/// Returns the number of paths stored in the graph
size_t graph_t::get_path_count(void) const {
    return _path_count;
}
    
/// Execute a function on each path in the graph
// TODO: allow stopping early?
void graph_t::for_each_path_handle(const std::function<void(const path_handle_t&)>& iteratee) const {
    for (uint64_t i = 0; i < paths.size(); ++i) {
        if (paths.at(i).occurrence_count()) {
            iteratee(as_path_handle(i));
        }
    }
}

void graph_t::for_each_occurrence_on_handle(const handle_t& handle, const std::function<void(const occurrence_handle_t&)>& iteratee) const {
    id_t id = graph_id_wt.at(handle_helper::unpack_number(handle));
    uint64_t i = path_id_wt.select(0, id)+1;
    uint64_t j = path_rank_wt.select(0, id)+1;
    uint64_t p = path_id_wt.at(i);
    uint64_t r = path_rank_wt.at(i);
    while (true) {
        // check for record termination
        if (p == 0 && r == 0) break;
        assert(p != 0 && r != 0); // both should be non-zero
        occurrence_handle_t occurrence_handle;
        as_integers(occurrence_handle)[0] = p-1;
        as_integers(occurrence_handle)[1] = r-1;
        // hit our callback
        iteratee(occurrence_handle);
        // update our values
        p = path_id_wt.at(++i);
        r = path_rank_wt.at(++j);
    }
}

/// Get a node handle (node ID and orientation) from a handle to an occurrence on a path
handle_t graph_t::get_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    //path_handle_t = as_path_handle(occ_handle[0]);
    auto& path = paths.at(occ_handle[0]);
    // get the step
    step_t step = path.get_occurrence(occ_handle[1]);
    // compute the handle
    handle_t handle = handle_helper::pack(graph_id_wt.select(0, step.id), step.strand);
    // TODO assert that it's fully live ?
    //assert(hidden_wt.rank(hidden_wt.size(), step.id) == 0);
    // if it is, we return it
    return handle;
}

/// Get a handle to the first occurrence in a path.
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_first_occurrence(const path_handle_t& path_handle) const {
    auto& path = paths.at(as_integer(path_handle));
    assert(path.occurrence_count());
    occurrence_handle_t result;
    int64_t* r_ints = as_integers(result);
    r_ints[0] = as_integer(path_handle);
    r_ints[1] = 0;
    return result;
}
    
/// Get a handle to the last occurrence in a path
/// The path MUST be nonempty.
occurrence_handle_t graph_t::get_last_occurrence(const path_handle_t& path_handle) const {
    auto& path = paths.at(as_integer(path_handle));
    assert(path.occurrence_count());
    occurrence_handle_t result;
    int64_t* r_ints = as_integers(result);
    r_ints[0] = as_integer(path_handle);
    r_ints[1] = path.occurrence_count()-1;
    return result;
}
    
/// Returns true if the occurrence is not the last occurence on the path, else false
bool graph_t::has_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    return occ_handle[1] < path.occurrence_count()-1;
}
    
/// Returns true if the occurrence is not the first occurence on the path, else false
bool graph_t::has_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    assert(path.occurrence_count());
    return occ_handle[1] > 0;
}

/// Returns a handle to the next occurrence on the path
occurrence_handle_t graph_t::get_next_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    assert(path.occurrence_count());
    uint64_t rank = occ_handle[1];
    occurrence_handle_t next_occ_handle = occurrence_handle;
    int64_t* next_i = as_integers(next_occ_handle);
    next_i[1] = rank+1;
    return next_occ_handle;
}

/// Returns a handle to the previous occurrence on the path
occurrence_handle_t graph_t::get_previous_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    auto& path = paths.at(occ_handle[0]);
    assert(path.occurrence_count());
    uint64_t rank = occ_handle[1];
    occurrence_handle_t prev_occ_handle = occurrence_handle;
    int64_t* prev_i = as_integers(prev_occ_handle);
    prev_i[1] = rank-1;
    return prev_occ_handle;
}
    
/// Returns a handle to the path that an occurrence is on
path_handle_t graph_t::get_path_handle_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    return as_path_handle(occ_handle[0]);
}
    
/// Returns the 0-based ordinal rank of a occurrence on a path
size_t graph_t::get_ordinal_rank_of_occurrence(const occurrence_handle_t& occurrence_handle) const {
    const int64_t* occ_handle = as_integers(occurrence_handle);
    return occ_handle[1];
}

////////////////////////////////////////////////////////////////////////////
// Additional optional interface with a default implementation
////////////////////////////////////////////////////////////////////////////

/// Returns true if the given path is empty, and false otherwise
bool graph_t::is_empty(const path_handle_t& path_handle) const {
    auto& path = paths.at(as_integer(path_handle));
    return path.occurrence_count() == 0;
}

////////////////////////////////////////////////////////////////////////////
// Concrete utility methods
////////////////////////////////////////////////////////////////////////////

/// Loop over all the occurrences along a path, from first through last
void graph_t::for_each_occurrence_in_path(const path_handle_t& path, const std::function<void(const occurrence_handle_t&)>& iteratee) const {
    if (is_empty(path)) return;
    occurrence_handle_t occ = get_first_occurrence(path);
    iteratee(occ); // run the first time
    while (has_next_occurrence(occ)) {
        occ = get_next_occurrence(occ);
        iteratee(occ);
    }
}

/**
 * This is the interface for a handle graph that supports modification.
 */
/*
 * Note: All operations may invalidate path handles and occurrence handles.
 */
    
/// Create a new node with the given sequence and return the handle.
handle_t graph_t::create_handle(const std::string& sequence) {
    // get node id as max+1
    return create_handle(sequence, _max_node_id+1);
}

/// Create a new node with the given id and sequence, then return the handle.
handle_t graph_t::create_handle(const std::string& sequence, const id_t& id) {
    assert(!graph_id_wt.char_exists(id));
    assert(id > 0);
    id_t new_id = id;
    // set new max
    _max_node_id = max(new_id, _max_node_id);
    // add to graph_id_wt
    graph_id_wt.push_back(new_id);
    // set up initial delimiters if the graph is empty
    if (!seq_wt.size()) {
        seq_wt.push_back(0);
        edge_fwd_wt.push_back(0);
        edge_fwd_inv_bv.push_back(0);
        edge_rev_wt.push_back(0);
        edge_rev_inv_bv.push_back(0);
        path_id_wt.push_back(0);
    }
    // append to seq_wt, delimit by 0
    for (auto c : sequence) seq_wt.push_back(c);
    seq_wt.push_back(0);
    // set up delemiters for edges, for later filling
    edge_fwd_wt.push_back(0);
    edge_fwd_inv_bv.push_back(0);
    edge_rev_wt.push_back(0);
    edge_rev_inv_bv.push_back(0);
    // set up delimiters for path_id_wt, to be filled later
    path_id_wt.push_back(0);
    // increment node count
    ++_node_count;
    // return handle
    return as_handle(new_id);
}
    
/// Remove the node belonging to the given handle and all of its edges.
/// Does not update any stored paths.
/// Invalidates the destroyed handle.
/// May be called during serial for_each_handle iteration **ONLY** on the node being iterated.
/// May **NOT** be called during parallel for_each_handle iteration.
/// May **NOT** be called on the node from which edges are being followed during follow_edges.
void graph_t::destroy_handle(const handle_t& handle) {
    uint64_t offset = handle_helper::unpack_number(handle);
    // remove from graph_id_wt
    id_t id = graph_id_wt.at(offset);
    graph_id_wt.remove(offset);
    // remove occs in edge lists
    while (edge_fwd_wt.rank(edge_fwd_wt.size(), id) > 0) {
        uint64_t i = edge_fwd_wt.select(0, id);
        edge_fwd_wt.remove(i);
        edge_fwd_inv_bv.remove(i);
        --_edge_count;
    }
    while (edge_rev_wt.rank(edge_rev_wt.size(), id) > 0) {
        uint64_t i = edge_rev_wt.select(0, id);
        edge_rev_wt.remove(i);
        edge_rev_inv_bv.remove(i);
    }
    // remove outbound edges and records from edge lists
    uint64_t edge_fwd_offset = edge_fwd_wt.select(offset, 0);
    edge_fwd_wt.remove(edge_fwd_offset);
    edge_fwd_inv_bv.remove(edge_fwd_offset);
    for (uint64_t i = edge_fwd_offset; edge_fwd_wt.at(i) != 0; ) {
        edge_fwd_wt.remove(i);
        edge_fwd_inv_bv.remove(i);
        --_edge_count;
    }
    uint64_t edge_rev_offset = edge_rev_wt.select(offset, 0);
    edge_rev_wt.remove(edge_rev_offset);
    edge_rev_inv_bv.remove(edge_rev_offset);
    for (uint64_t i = edge_rev_offset; edge_rev_wt.at(i) != 0; ) {
        edge_rev_wt.remove(i);
        edge_rev_inv_bv.remove(i);
    }
    // save the node sequence for stashing in the paths
    std::string seq = get_sequence(handle);
    // remove the sequence from seq_wt
    uint64_t seq_wt_offset = seq_wt.select(offset, 0);
    for (uint64_t i = seq_wt_offset; seq_wt.at(i) != 0; ++i) {
        seq_wt.remove(i);
    }
    // move the sequence of the node into each path that traverses it
    // remove reference to the node from the paths
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            auto& path = paths.at(as_integers(occ)[0]);
            path.unlink_occurrence(as_integers(occ)[1], seq);
        });
    // remove reference to the paths
    // ... from the node to path id mapping
    uint64_t path_id_offset = path_id_wt.select(offset, 0);
    path_id_wt.remove(path_id_offset);
    for (uint64_t i = path_id_offset; path_id_wt.at(i) != 0; ++i) {
        path_id_wt.remove(i);
    }
    // ... from the node to path rank mapping
    uint64_t path_rank_offset = path_rank_wt.select(offset, 0);
    path_rank_wt.remove(path_rank_offset);
    for (uint64_t i = path_rank_offset; path_rank_wt.at(i) != 0; ++i) {
        path_rank_wt.remove(i);
    }
    --_node_count;
}
    
/// Create an edge connecting the given handles in the given order and orientations.
/// Ignores existing edges.
void graph_t::create_edge(const handle_t& left, const handle_t& right) {
    if (has_edge(left, right)) return; // do nothing if edge exists
    uint64_t left_rank = handle_helper::unpack_number(left);
    bool left_rev = handle_helper::unpack_bit(left);
    uint64_t right_rank = handle_helper::unpack_number(right);
    bool right_rev = handle_helper::unpack_bit(right);
    // canonicalize
    if (left_rev && right_rev) {
        left_rev = !left_rev;
        right_rev = !right_rev;
        std::swap(left_rank, right_rank);
    }
    bool inv = (left_rev != right_rev);
    if (!left_rev) {
        uint64_t edge_fwd_left_offset = edge_fwd_wt.select(0, left_rank);
        edge_fwd_wt.insert(edge_fwd_left_offset+1, right_rank+1);
        edge_fwd_inv_bv.insert(edge_fwd_left_offset+1, inv);
    } else {
        uint64_t edge_rev_left_offset = edge_rev_wt.select(0, left_rank);
        edge_rev_wt.insert(edge_rev_left_offset+1, right_rank+1);
        edge_rev_inv_bv.insert(edge_rev_left_offset+1, inv);
    }
    if (!right_rev) {
        uint64_t edge_rev_right_offset = edge_rev_wt.select(0, right_rank);
        edge_rev_wt.insert(edge_rev_right_offset+1, left_rank+1);
        edge_rev_inv_bv.insert(edge_rev_right_offset+1, inv);
    } else {
        uint64_t edge_fwd_right_offset = edge_fwd_wt.select(0, right_rank);
        edge_rev_wt.insert(edge_fwd_right_offset+1, left_rank+1);
        edge_rev_inv_bv.insert(edge_fwd_right_offset+1, inv);
    }
    ++_edge_count;
}

bool graph_t::has_edge(const handle_t& left, const handle_t& right) {
    bool exists = false;
    follow_edges(left, false, [&right, &exists](const handle_t& next) {
            if (next == right) exists = true;
        });
    return exists;
}
    
/// Remove the edge connecting the given handles in the given order and orientations.
/// Ignores nonexistent edges.
/// Does not update any stored paths.
void graph_t::destroy_edge(const handle_t& left, const handle_t& right) {
    uint64_t left_rank = handle_helper::unpack_number(left);
    bool left_rev = handle_helper::unpack_bit(left);
    uint64_t right_rank = handle_helper::unpack_number(right);
    bool right_rev = handle_helper::unpack_bit(right);
    // canonicalize
    if (left_rev && right_rev) {
        left_rev = !left_rev;
        right_rev = !right_rev;
        std::swap(left_rank, right_rank);
    }
    bool inv = (left_rev != right_rev);
    if (!left_rev) {
        uint64_t edge_fwd_left_offset = edge_fwd_wt.select(0, left_rank);
        uint64_t edge_fwd_left_offset_erase = 0;
        for (uint64_t i = edge_fwd_left_offset+1; ; ++i) {
            uint64_t c = edge_fwd_wt.at(i);
            if (c != 0) break;
            if (c == right_rank+1) {
                edge_fwd_left_offset_erase = i;
                break;
            }
        }
        if (edge_fwd_left_offset_erase) {
            edge_fwd_wt.remove(edge_fwd_left_offset_erase);
            edge_fwd_inv_bv.remove(edge_fwd_left_offset_erase);
        }
    } else {
        uint64_t edge_rev_left_offset = edge_rev_wt.select(0, left_rank);
        uint64_t edge_rev_left_offset_erase = 0;
        for (uint64_t i = edge_rev_left_offset+1; ; ++i) {
            uint64_t c = edge_rev_wt.at(i);
            if (c != 0) break;
            if (c == right_rank+1) {
                edge_rev_left_offset_erase = i;
                break;
            }
        }
        if (edge_rev_left_offset_erase) {
            edge_rev_wt.remove(edge_rev_left_offset_erase);
            edge_rev_inv_bv.remove(edge_rev_left_offset_erase);
        }
    }
    if (!right_rev) {
        uint64_t edge_rev_right_offset = edge_rev_wt.select(0, right_rank);
        uint64_t edge_rev_right_offset_erase = 0;
        for (uint64_t i = edge_rev_right_offset+1; ; ++i) {
            uint64_t c = edge_rev_wt.at(i);
            if (c != 0) break;
            if (c == right_rank+1) {
                edge_rev_right_offset_erase = i;
                break;
            }
        }
        if (edge_rev_right_offset_erase) {
            edge_rev_wt.remove(edge_rev_right_offset_erase);
            edge_rev_inv_bv.remove(edge_rev_right_offset_erase);
        }
    } else {
        uint64_t edge_fwd_right_offset = edge_fwd_wt.select(0, right_rank);
        uint64_t edge_fwd_right_offset_erase = 0;
        for (uint64_t i = edge_fwd_right_offset+1; ; ++i) {
            uint64_t c = edge_fwd_wt.at(i);
            if (c != 0) break;
            if (c == right_rank+1) {
                edge_fwd_right_offset_erase = i;
                break;
            }
        }
        if (edge_fwd_right_offset_erase) {
            edge_fwd_wt.remove(edge_fwd_right_offset_erase);
            edge_fwd_inv_bv.remove(edge_fwd_right_offset_erase);
        }
    }
    --_edge_count;
}
        
/// Remove all nodes and edges. Does not update any stored paths.
void graph_t::clear(void) {
    dyn::wt_string<dyn::suc_bv> null_wt;
    dyn::suc_bv null_bv;
    dyn::wt_fmi null_fmi;
    _max_node_id = 0;
    _min_node_id = 0;
    graph_id_wt = null_wt;
    edge_fwd_wt = null_wt;
    edge_fwd_inv_bv = null_bv;
    edge_rev_wt = null_wt;
    edge_rev_inv_bv = null_bv;
    seq_wt = null_wt;
    path_id_wt = null_wt;
    path_rank_wt = null_wt;
    path_name_fmi = null_fmi;
    path_name_bv = null_bv;
    paths.clear();
    _node_count = 0;
    _edge_count = 0;
    _path_count = 0;
}
    
/// Swap the nodes corresponding to the given handles, in the ordering used
/// by for_each_handle when looping over the graph. Other handles to the
/// nodes being swapped must not be invalidated. If a swap is made while
/// for_each_handle is running, it affects the order of the handles
/// traversed during the current traversal (so swapping an already seen
/// handle to a later handle's position will make the seen handle be visited
/// again and the later handle not be visited at all).
void graph_t::swap_handles(const handle_t& a, const handle_t& b) {
    assert(false);
}
    
/// Alter the node that the given handle corresponds to so the orientation
/// indicated by the handle becomes the node's local forward orientation.
/// Rewrites all edges pointing to the node and the node's sequence to
/// reflect this. Invalidates all handles to the node (including the one
/// passed). Returns a new, valid handle to the node in its new forward
/// orientation. Note that it is possible for the node's ID to change.
/// Updates all stored paths. May change the ordering of the underlying
/// graph.
handle_t graph_t::apply_orientation(const handle_t& handle) {
    // do nothing if we're already in the right orientation
    if (!handle_helper::unpack_bit(handle)) return handle;
    // store edges
    vector<handle_t> edges_fwd;
    vector<handle_t> edges_rev;
    follow_edges(handle, false, [&](const handle_t& h) {
            edges_fwd.push_back(h);
        });
    follow_edges(handle, true, [&](const handle_t& h) {
            edges_rev.push_back(h);
        });
    // save the sequence's reverse complement, which we will use to add the new handle
    const std::string seq = (handle_helper::unpack_bit(handle)
                             ? reverse_complement(get_sequence(handle))
                             : get_sequence(handle));
    // get the path set
    //vector<path_occ
    vector<occurrence_handle_t> occurrences;
    vector<bool> orientations;
    for_each_occurrence_on_handle(handle, [&](const occurrence_handle_t& occ) {
            // save occurrence information
            occurrences.push_back(occ);
            auto& path = paths.at(as_integers(occ)[0]);
            // record relative orientation
            orientations.push_back(path.get_occurrence(as_integers(occ)[1]).strand);
        });
    // destroy the handle!
    destroy_handle(handle);
    // we have the technology. we can rebuild it.
    handle_t new_handle = create_handle(seq);
    // reconnect it to the graph
    for (auto& h : edges_fwd) {
        create_edge(handle_helper::toggle_bit(new_handle), h);
    }
    for (auto& h : edges_rev) {
        create_edge(h, handle_helper::toggle_bit(new_handle));
    }
    // and relink it to the paths
    for (uint64_t j = 0; j < occurrences.size(); ++j) {
        auto& occ = occurrences.at(j);
        bool orient = orientations.at(j);
        handle_t h = orient ? new_handle : handle_helper::toggle_bit(new_handle);
        auto& path = paths.at(as_integers(occ)[0]);
        path.link_occurrence(as_integers(occ)[1], h, seq);
    }
}
    
/// Split a handle's underlying node at the given offsets in the handle's
/// orientation. Returns all of the handles to the parts. Other handles to
/// the node being split may be invalidated. The split pieces stay in the
/// same local forward orientation as the original node, but the returned
/// handles come in the order and orientation appropriate for the handle
/// passed in.
/// Updates stored paths.
std::vector<handle_t> graph_t::divide_handle(const handle_t& handle, const std::vector<size_t>& offsets) {
    
}
    

/**
 * This is the interface for a handle graph with embedded paths where the paths can be modified.
 * Note that if the *graph* can also be modified, the implementation will also
 * need to inherit from MutableHandleGraph, via the combination
 * MutablePathMutableHandleGraph interface.
 * TODO: This is a very limited interface at the moment. It will probably need to be extended.
 */
    
/**
 * Destroy the given path. Invalidates handles to the path and its node occurrences.
 */
void graph_t::destroy_path(const path_handle_t& path) {
    paths[as_integer(path)].clear();
}

/**
 * Create a path with the given name. The caller must ensure that no path
 * with the given name exists already, or the behavior is undefined.
 * Returns a handle to the created empty path. Handles to other paths must
 * remain valid.
 */
path_handle_t graph_t::create_path_handle(const std::string& name) {
    paths.emplace_back(name);
}
    
/**
 * Append a visit to a node to the given path. Returns a handle to the new
 * final occurrence on the path which is appended. Handles to prior
 * occurrences on the path, and to other paths, must remain valid.
 */
occurrence_handle_t graph_t::append_occurrence(const path_handle_t& path, const handle_t& to_append) {
    //void append_occurrence(const id_t& id, bool strand);
    paths[as_integer(path)].append_occurrence(to_append);
}


}