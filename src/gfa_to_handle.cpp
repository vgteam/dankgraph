#include "gfa_to_handle.hpp"

namespace odgi {

std::map<char, uint64_t> gfa_line_counts(const char* filename) {
    int gfa_fd = -1;
    char* gfa_buf = nullptr;
    size_t gfa_filesize = gfak::mmap_open(filename, gfa_buf, gfa_fd);
    if (gfa_fd == -1) {
        cerr << "Couldn't open GFA file " << filename << "." << endl;
        exit(1);
    }
    string line;
    size_t i = 0;
    bool seen_newline = true;
    std::map<char, uint64_t> counts;
    while (i < gfa_filesize) {
        if (i == 0 || gfa_buf[i-1] == '\n') {
            counts[gfa_buf[i]]++;
        }
        ++i;
    }
    gfak::mmap_close(gfa_buf, gfa_fd, gfa_filesize);
    return counts;
}

void gfa_to_handle(const string& gfa_filename,
                   handlegraph::MutablePathMutableHandleGraph* graph,
                   uint64_t n_threads,
                   bool progress) {

    n_threads = (n_threads == 0 ? 1 : n_threads);
    char* filename = (char*) gfa_filename.c_str();
    //std::cerr << "filename is " << filename << std::endl;
    gfak::GFAKluge gg;
    uint64_t i = 0;
    uint64_t min_id = std::numeric_limits<uint64_t>::max();
    uint64_t max_id = std::numeric_limits<uint64_t>::min();
    std::map<char, uint64_t> line_counts;
    // in parallel scan over the file to count edges and sequences
    {
        std::thread x(
            [&](void) {
                gg.for_each_sequence_line_in_file(
                    filename,
                    [&](gfak::sequence_elem s) {
                        uint64_t id = stol(s.name);
                        min_id = std::min(min_id, id);
                        max_id = std::max(max_id, id);
                    });
            });
        line_counts = gfa_line_counts(filename);
        x.join();
    }
    uint64_t id_increment = min_id - 1;
    uint64_t node_count = line_counts['S'];
    uint64_t edge_count = line_counts['L'];
    uint64_t path_count = line_counts['P'];
    // set the min id as an offset
    graph->set_id_increment(id_increment);
    // build the nodes
    {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
        if (progress) {
            progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                node_count, "[odgi::gfa_to_handle] building nodes:");
        }
        gg.for_each_sequence_line_in_file(
            filename,
            [&](gfak::sequence_elem s) {
                const uint64_t id = stol(s.name);
                graph->create_handle(s.sequence, id - id_increment);
                if (progress) progress_meter->increment(1);
            });
        if (progress) {
            progress_meter->finish();
        }
    }
    
    {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
        if (progress) {
            progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                edge_count, "[odgi::gfa_to_handle] building edges:");
        }
        gg.for_each_edge_line_in_file(
            filename,
            [&](gfak::edge_elem e) {
                if (e.source_name.empty()) return;
                const handlegraph::handle_t a = graph->get_handle(stol(e.source_name), !e.source_orientation_forward);
                const handlegraph::handle_t b = graph->get_handle(stol(e.sink_name), !e.sink_orientation_forward);
                graph->create_edge(a, b);
                if (progress) progress_meter->increment(1);
            });
        if (progress) {
            progress_meter->finish();
        }
    }

    if (path_count > 0) {
        std::unique_ptr<algorithms::progress_meter::ProgressMeter> progress_meter;
        if (progress) {
            progress_meter = std::make_unique<algorithms::progress_meter::ProgressMeter>(
                path_count, "[odgi::gfa_to_handle] building paths:");
        }
        gfa_path_queue_t path_queue;
        std::mutex logging_mutex;
        std::atomic<bool> work_todo;
        uint64_t idx = 0;
        auto worker =
            [&](uint64_t tid) {
                while (work_todo.load()) {
                    path_elem_t * p;
                    if (path_queue.try_pop(p)) {
                        uint64_t i = 0;
                        for (auto& s : p->gfak.segment_names) {
                            graph->append_step(p->path,
                                               graph->get_handle(std::stoi(s),
                                                                 // in gfak, true == +
                                                                 !p->gfak.orientations[i++]));
                        }
                        delete p;
                        if (progress) progress_meter->increment(1);
                    } else {
                        std::this_thread::sleep_for(std::chrono::nanoseconds(1));
                    }
                }
            };

        std::vector<std::thread> workers;
        workers.reserve(n_threads);
        work_todo.store(true);
        for (uint64_t t = 0; t < n_threads; ++t) {
            workers.emplace_back(worker, t);
        }

        gg.for_each_path_line_in_file(
            filename,
            [&](const gfak::path_elem& path) {
                handlegraph::path_handle_t p_h = graph->create_path_handle(path.name);
                path_elem_t* p = new path_elem_t({p_h, path});
                path_queue.push(p);
            });

        while (!path_queue.was_empty()) {
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        }
        work_todo.store(false);
        for (uint64_t t = 0; t < n_threads; ++t) {
            workers[t].join();
        }
        if (progress) {
            progress_meter->finish();
        }
    }
}
}
