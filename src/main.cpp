#include <CLI/CLI.hpp>

#ifdef TET_SHELL_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <tetshell/Logger.hpp>
#include <tetshell/Parameters.hpp>
#include <Eigen/Dense>

#include <igl/Timer.h>
#include <igl/write_triangle_mesh.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>

#include<bitset>

using namespace tetshell;
using namespace Eigen;

class GeoLoggerForward: public GEO::LoggerClient {
    std::shared_ptr<spdlog::logger> logger_;

public:
    template<typename T>
    GeoLoggerForward(T logger) : logger_(logger) {}

private:
    std::string truncate(const std::string &msg) {
        static size_t prefix_len = GEO::CmdLine::ui_feature(" ", false).size();
        return msg.substr(prefix_len, msg.size() - 1 - prefix_len);
    }

protected:
    void div(const std::string &title) override {
        logger_->trace(title.substr(0, title.size() - 1));
    }

    void out(const std::string &str) override {
        logger_->info(truncate(str));
    }

    void warn(const std::string &str) override {
        logger_->warn(truncate(str));
    }

    void err(const std::string &str) override {
        logger_->error(truncate(str));
    }

    void status(const std::string &str) override {
        // Errors and warnings are also dispatched as status by geogram, but without
        // the "feature" header. We thus forward them as trace, to avoid duplicated
        // logger info...
        logger_->trace(str.substr(0, str.size() - 1));
    }
};


//extern "C" void exactinit();
int main(int argc, char **argv) {

#ifndef WIN32
    setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

    GEO::initialize();

    Parameters params;

    // Import standard command line arguments, and custom ones
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("algo");

    CLI::App command_line{"tet-shell"};
    command_line.add_option("-i,--input", params.input_path,
                            "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")->check(
            CLI::ExistingFile);

    unsigned int max_threads = std::numeric_limits<unsigned int>::max();
#ifdef TET_SHELL_USE_TBB
    command_line.add_option("--max-threads", max_threads, "maximum number of threads used");
#endif

    try {
        command_line.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return command_line.exit(e);
    }


#ifdef TET_SHELL_USE_TBB
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads = std::min(max_threads, num_threads);
    params.num_threads = num_threads;
    std::cout << "TBB threads " << num_threads << std::endl;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

    Logger::init(!params.is_quiet, params.log_path);
    params.log_level = std::max(0, std::min(6, params.log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    GEO::Logger *geo_logger = GEO::Logger::instance();
    geo_logger->unregister_all_clients();
    geo_logger->register_client(new GeoLoggerForward(logger().clone("geogram")));
    geo_logger->set_pretty(false);


    if (params.output_path.empty())
        params.output_path = params.input_path;
    if (params.log_path.empty())
        params.log_path = params.output_path;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    logger().info("EOP");
    std::cerr << "test";

    return EXIT_SUCCESS;
}
