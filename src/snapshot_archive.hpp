#pragma once
#include <zlib.h>
#include <algorithm>
#include <array>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>

// Stream a POSIX ustar archive through gzip, one snapshot at a time.
class SnapshotArchive {
    static constexpr std::size_t BLOCK_SIZE = 512;
    static constexpr unsigned int WRITE_CHUNK = 1024 * 1024;
    static constexpr unsigned int FILE_MODE = 0644;
    static constexpr std::size_t NAME_SIZE = 100;
    static constexpr unsigned long long MAX_FILE_SIZE = 077777777777ULL;
    gzFile file_;
    void write(const char* data, std::size_t size) {
        while (size) {
            const auto count = static_cast<unsigned int>(std::min<std::size_t>(size, WRITE_CHUNK));
            if (gzwrite(file_, data, count) != static_cast<int>(count))
                throw std::runtime_error("Cannot write compressed snapshots.");
            data += count;
            size -= count;
        }
    }
public:
    explicit SnapshotArchive(const std::filesystem::path& path)
        : file_(gzopen(path.string().c_str(), "wb")) {
        if (!file_) throw std::runtime_error("Cannot open configuration archive.");
    }
    SnapshotArchive(const SnapshotArchive&) = delete;
    SnapshotArchive& operator=(const SnapshotArchive&) = delete;
    ~SnapshotArchive() { if (file_) gzclose(file_); }
    void add(const std::string& name, const std::string& contents) {
        if (name.size() >= NAME_SIZE || contents.size() > MAX_FILE_SIZE)
            throw std::runtime_error("Snapshot exceeds ustar format limits.");
        std::array<char, BLOCK_SIZE> header{};
        std::snprintf(header.data(), NAME_SIZE, "%s", name.c_str());
        std::snprintf(header.data()+100, 8, "%07o", FILE_MODE);
        std::snprintf(header.data()+108, 8, "%07o", 0);
        std::snprintf(header.data()+116, 8, "%07o", 0);
        std::array<char, 32> size_text{};
        std::snprintf(size_text.data(), size_text.size(), "%011llo",
                      static_cast<unsigned long long>(contents.size()));
        std::copy_n(size_text.data(), 12, header.begin()+124);
        std::snprintf(header.data()+136, 12, "%011o", 0);
        std::fill(header.begin()+148, header.begin()+156, ' ');
        header[156]='0';
        std::copy_n("ustar", 5, header.begin()+257);
        header[263]='0'; header[264]='0';
        unsigned int checksum=0;
        for (unsigned char value : header) checksum+=value;
        std::snprintf(header.data()+148, 7, "%06o", checksum);
        header[155]=' ';
        write(header.data(), header.size());
        write(contents.data(), contents.size());
        const std::array<char, BLOCK_SIZE> padding{};
        const auto remainder=contents.size()%BLOCK_SIZE;
        if (remainder) write(padding.data(), BLOCK_SIZE-remainder);
        // Keep completed snapshots readable after interruption where possible.
        if (gzflush(file_, Z_SYNC_FLUSH)!=Z_OK)
            throw std::runtime_error("Cannot flush configuration archive.");
    }
    void close() {
        const std::array<char, BLOCK_SIZE * 2> end{};
        write(end.data(), end.size());
        const int status=gzclose(file_);
        file_=nullptr;
        if (status!=Z_OK) throw std::runtime_error("Cannot finish configuration archive.");
    }
};
