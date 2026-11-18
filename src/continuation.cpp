#include "gpe/GPE_solver_contd.hpp"
#include <archive.h>
#include <archive_entry.h>
#include <array>
#include <cmath>
#include <fstream>
#include <istream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <streambuf>

namespace {
constexpr std::size_t READ_BUFFER_SIZE = 64 * 1024;
constexpr double GRID_TOLERANCE = 1e-9;
constexpr int MIN_GRID_POINTS = 3;

void check_archive(int status, archive* input) {
    if (status != ARCHIVE_OK) {
        const char* message = archive_error_string(input);
        throw std::runtime_error(message ? message : "Cannot read compressed checkpoint.");
    }
}

class CheckpointBuffer : public std::streambuf {
    archive* input_;
    std::array<char, READ_BUFFER_SIZE> buffer_;
    int_type underflow() override {
        const auto count = archive_read_data(input_, buffer_.data(), buffer_.size());
        if (count < 0) check_archive(ARCHIVE_FATAL, input_);
        if (count == 0) return traits_type::eof();
        setg(buffer_.data(), buffer_.data(), buffer_.data() + count);
        return traits_type::to_int_type(*gptr());
    }
public:
    explicit CheckpointBuffer(archive* input) : input_(input) {}
};

Eigen::VectorXcd read_field(archive* input, int points, double dx, double length) {
    CheckpointBuffer buffer(input);
    std::istream text(&buffer);
    text.exceptions(std::ios::badbit);
    Eigen::VectorXcd field(points);
    for (int i = 0; i < points; ++i) {
        double x, re, im;
        if (!(text >> x >> re >> im) || !std::isfinite(x) || !std::isfinite(re)
            || !std::isfinite(im) || std::abs(x - (i * dx - length / 2.0)) > GRID_TOLERANCE)
            throw std::runtime_error("Invalid checkpoint data or grid mismatch.");
        field(i) = std::complex<double>(re, im);
    }
    text >> std::ws;
    if (!text.eof()) throw std::runtime_error("Checkpoint has extra data.");
    return field;
}
}

gpe::SimulationResult gpe::GPE_solver_contd(
    const SimulationParameters& parameters,
    const std::filesystem::path& checkpoint_file,
    const std::string& checkpoint_member) {
    if (!std::isfinite(parameters.t_initial) || parameters.t_initial <= 0
        || !parameters.initial_velocity)
        throw std::runtime_error("Compressed continuation requires t_initial > 0 and initial_velocity.");
    if (parameters.initial_field.size())
        throw std::runtime_error("Leave initial_field empty when loading a compressed checkpoint.");
    const double points = parameters.length / parameters.dx;
    if (!std::isfinite(parameters.length) || !std::isfinite(parameters.dx)
        || parameters.length <= 0 || parameters.dx <= 0 || !std::isfinite(points)
        || points < MIN_GRID_POINTS || points > std::numeric_limits<int>::max()
        || std::abs(points - std::round(points)) > GRID_TOLERANCE)
        throw std::runtime_error("Invalid checkpoint grid.");

    std::unique_ptr<archive, decltype(&archive_read_free)> input(archive_read_new(), archive_read_free);
    if (!input) throw std::runtime_error("Cannot allocate checkpoint reader.");
    // Require built-in filters; ARCHIVE_WARN would allow an external helper.
    check_archive(archive_read_support_filter_gzip(input.get()), input.get());
    check_archive(archive_read_support_filter_xz(input.get()), input.get());
    check_archive(archive_read_support_format_tar(input.get()), input.get());
    check_archive(archive_read_support_format_raw(input.get()), input.get());
    check_archive(archive_read_open_filename(input.get(), checkpoint_file.string().c_str(), READ_BUFFER_SIZE), input.get());

    SimulationParameters settings = parameters;
    archive_entry* entry = nullptr;
    bool found = false;
    int status;
    while ((status = archive_read_next_header(input.get(), &entry)) != ARCHIVE_EOF) {
        check_archive(status, input.get());
        const int filter = archive_filter_code(input.get(), 0);
        if (filter != ARCHIVE_FILTER_GZIP && filter != ARCHIVE_FILTER_XZ)
            throw std::runtime_error("Checkpoint must be gzip- or XZ-compressed.");
        const bool raw = archive_format(input.get()) == ARCHIVE_FORMAT_RAW;
        if (raw && !checkpoint_member.empty())
            throw std::runtime_error("A single compressed field has no archive member name.");
        if (!raw && checkpoint_member.empty())
            throw std::runtime_error("Specify checkpoint_member for a compressed tar archive.");
        if (raw || checkpoint_member == archive_entry_pathname(entry)) {
            if (found) throw std::runtime_error("Duplicate checkpoint member in archive.");
            if (archive_entry_filetype(entry) != AE_IFREG || archive_entry_hardlink(entry))
                throw std::runtime_error("Checkpoint member must be a regular file.");
            settings.initial_field = read_field(input.get(), static_cast<int>(std::round(points)),
                                                parameters.dx, parameters.length);
            found = true;
        } else {
            check_archive(archive_read_data_skip(input.get()), input.get());
        }
    }
    if (!found) throw std::runtime_error("Checkpoint member not found.");
    check_archive(archive_read_close(input.get()), input.get());
    input.reset();
    auto result = GPE_solver(settings);
    std::ofstream metadata(settings.output_directory / "parameters.txt", std::ios::app);
    metadata.exceptions(std::ios::badbit | std::ios::failbit);
    metadata << "restart " << checkpoint_file.string() << "\n";
    if (!checkpoint_member.empty()) metadata << "restart_member " << checkpoint_member << "\n";
    metadata.close();
    return result;
}
