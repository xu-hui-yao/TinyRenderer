#pragma once

#include <algorithm>
#include <cctype>
#include <filesystem/fwd.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <ShlObj.h>
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>

#if defined(__linux)
#include <linux/limits.h>
#endif

NAMESPACE_BEGIN(filesystem)
class path {
public:
    enum path_type {
        windows_path = 0,
        posix_path   = 1,
#if defined(_WIN32)
        native_path = windows_path
#else
        native_path = posix_path
#endif
    };

    path() : m_type(native_path) {}

    path(const path &path) = default;

    path(path &&path) noexcept
        : m_type(path.m_type), m_path(std::move(path.m_path)), m_absolute(path.m_absolute), m_smb(path.m_smb) {}

    explicit path(const char *string) { set(string); }

    explicit path(const std::string &string) { set(string); }

#if defined(_WIN32)
    explicit path(const std::wstring &wstring) { set(wstring); }
    explicit path(const wchar_t *wstring) { set(wstring); }
#endif

    [[nodiscard]] size_t length() const { return m_path.size(); }

    [[nodiscard]] bool empty() const { return m_path.empty(); }

    [[nodiscard]] bool is_absolute() const { return m_absolute; }

    [[nodiscard]] path make_absolute() const {
#if !defined(_WIN32)
        char temp[PATH_MAX];
        if (realpath(str().c_str(), temp) == nullptr)
            throw std::runtime_error("Internal error in realpath(): " + std::string(strerror(errno)));
        return path(temp);
#else
        std::wstring value = wstr(), out(m_max_path_windows, '\0');
        DWORD length       = GetFullPathNameW(value.c_str(), m_max_path_windows, &out[0], nullptr);
        if (length == 0)
            throw std::runtime_error("Internal error in realpath(): " + std::to_string(GetLastError()));
        return path(out.substr(0, length));
#endif
    }

    [[nodiscard]] bool exists() const {
#if defined(_WIN32)
        return GetFileAttributesW(wstr().c_str()) != INVALID_FILE_ATTRIBUTES;
#else
        struct stat sb;
        return stat(str().c_str(), &sb) == 0;
#endif
    }

    [[nodiscard]] size_t file_size() const {
#if defined(_WIN32)
        struct _stati64 sb {};
        if (_wstati64(wstr().c_str(), &sb) != 0)
            throw std::runtime_error("path::file_size(): cannot stat file \"" + str() + "\"!");
#else
        struct stat sb;
        if (stat(str().c_str(), &sb) != 0)
            throw std::runtime_error("path::file_size(): cannot stat file \"" + str() + "\"!");
#endif
        return static_cast<size_t>(sb.st_size);
    }

    [[nodiscard]] bool is_directory() const {
#if defined(_WIN32)
        DWORD result = GetFileAttributesW(wstr().c_str());
        if (result == INVALID_FILE_ATTRIBUTES)
            return false;
        return (result & FILE_ATTRIBUTE_DIRECTORY) != 0;
#else
        struct stat sb;
        if (stat(str().c_str(), &sb))
            return false;
        return S_ISDIR(sb.st_mode);
#endif
    }

    [[nodiscard]] bool is_file() const {
#if defined(_WIN32)
        DWORD attr = GetFileAttributesW(wstr().c_str());
        return (attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY) == 0);
#else
        struct stat sb;
        if (stat(str().c_str(), &sb))
            return false;
        return S_ISREG(sb.st_mode);
#endif
    }

    [[nodiscard]] std::string extension() const {
        const std::string &name = filename();
        size_t pos              = name.find_last_of('.');
        if (pos == std::string::npos)
            return "";
        return name.substr(pos + 1);
    }

    [[nodiscard]] std::string filename() const {
        if (empty())
            return "";
        const std::string &last = m_path[m_path.size() - 1];
        return last;
    }

    [[nodiscard]] path parent_path() const {
        path result;
        result.m_absolute = m_absolute;
        result.m_smb      = m_smb;

        if (m_path.empty()) {
            if (!m_absolute)
                result.m_path.emplace_back("..");
        } else {
            size_t until = m_path.size() - 1;
            for (size_t i = 0; i < until; ++i)
                result.m_path.push_back(m_path[i]);
        }
        return result;
    }

    path operator/(const path &other) const {
        if (other.m_absolute)
            throw std::runtime_error("path::operator/(): expected a relative path!");
        if (m_type != other.m_type)
            throw std::runtime_error("path::operator/(): expected a path of the same type!");

        path result(*this);

        for (const auto &i : other.m_path)
            result.m_path.push_back(i);

        return result;
    }

    [[nodiscard]] std::string str(path_type type = native_path) const {
        std::ostringstream oss;

        if (m_absolute) {
            if (m_type == posix_path)
                oss << "/";
            else {
                size_t length = 0;
                for (const auto &i : m_path)
                    // No special case for the last segment to count the nullptr character
                    length += i.length() + 1;
                if (m_smb)
                    length += 2;

                // Windows requires a \\?\ prefix to handle paths longer than MAX_PATH
                // (including their nullptr character). NOTE: relative paths >MAX_PATH are
                // not supported at all in Windows.
                if (length > m_max_path_windows_legacy) {
                    if (m_smb)
                        oss << R"(\\?\UNC\)";
                    else
                        oss << R"(\\?\)";
                } else if (m_smb)
                    oss << "\\\\";
            }
        }

        for (size_t i = 0; i < m_path.size(); ++i) {
            oss << m_path[i];
            if (i + 1 < m_path.size()) {
                if (type == posix_path)
                    oss << '/';
                else
                    oss << '\\';
            }
        }

        return oss.str();
    }

    void set(const std::string &str, path_type type = native_path) {
        m_type = type;
        if (type == windows_path) {
            std::string tmp = str;

            // Long windows paths (sometimes) begin with the prefix \\?\. It should only
            // be used when the path is >MAX_PATH characters long, so we remove it
            // for convenience and add it back (if necessary) in str()/wstr().
            static const std::string long_path_prefix = R"(\\?\)";
            if (tmp.length() >= long_path_prefix.length() &&
                std::mismatch(std::begin(long_path_prefix), std::end(long_path_prefix), std::begin(tmp)).first ==
                    std::end(long_path_prefix)) {
                tmp.erase(0, long_path_prefix.length());
            }

            // Special-case handling of absolute SMB paths, which start with the prefix "\\".
            if (tmp.length() >= 2 && tmp[0] == '\\' && tmp[1] == '\\') {
                m_path = {};
                tmp.erase(0, 2);

                // Interestingly, there is a special-special case where relative paths may be specified as beginning
                // with a "\\" when a non-SMB file with a more-than-260-characters-long absolute _local_ path is
                // double-clicked. This seems to only happen with single-segment relative paths, so we can check for
                // this condition by making sure no further path separators are present.
                if (tmp.find_first_of("/\\") != std::string::npos)
                    m_absolute = m_smb = true;
                else
                    m_absolute = m_smb = false;

                // Special-case handling of absolute SMB paths, which start with the prefix "UNC\"
            } else if (tmp.length() >= 4 && tmp[0] == 'U' && tmp[1] == 'N' && tmp[2] == 'C' && tmp[3] == '\\') {
                m_path = {};
                tmp.erase(0, 4);
                m_absolute = true;
                m_smb      = true;
                // Special-case handling of absolute local paths, which start with the drive letter and a colon "X:\"
            } else if (tmp.length() >= 3 && std::isalpha(tmp[0]) && tmp[1] == ':' &&
                       (tmp[2] == '\\' || tmp[2] == '/')) {
                m_path = { tmp.substr(0, 2) };
                tmp.erase(0, 3);
                m_absolute = true;
                m_smb      = false;
                // Relative path
            } else {
                m_path     = {};
                m_absolute = false;
                m_smb      = false;
            }

            std::vector<std::string> tokenized = tokenize(tmp, "/\\");
            m_path.insert(std::end(m_path), std::begin(tokenized), std::end(tokenized));
        } else {
            m_path     = tokenize(str, "/");
            m_absolute = !str.empty() && str[0] == '/';
        }
    }

    path &operator=(const path &path) = default;

    path &operator=(path &&path) noexcept {
        if (this != &path) {
            m_type     = path.m_type;
            m_path     = std::move(path.m_path);
            m_absolute = path.m_absolute;
            m_smb      = path.m_smb;
        }
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &os, const path &path) {
        os << path.str();
        return os;
    }

    [[nodiscard]] bool remove_file() const {
#if !defined(_WIN32)
        return std::remove(str().c_str()) == 0;
#else
        return DeleteFileW(wstr().c_str()) != 0;
#endif
    }

    [[nodiscard]] bool resize_file(size_t target_length) const {
#if !defined(_WIN32)
        return ::truncate(str().c_str(), (off_t) target_length) == 0;
#else
        HANDLE handle = CreateFileW(wstr().c_str(), GENERIC_WRITE, 0, nullptr, 0, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (handle == INVALID_HANDLE_VALUE)
            return false;
        LARGE_INTEGER size;
        size.QuadPart = static_cast<LONGLONG>(target_length);
        if (SetFilePointerEx(handle, size, nullptr, FILE_BEGIN) == 0) {
            CloseHandle(handle);
            return false;
        }
        if (SetEndOfFile(handle) == 0) {
            CloseHandle(handle);
            return false;
        }
        CloseHandle(handle);
        return true;
#endif
    }

    static path getcwd() {
#if !defined(_WIN32)
        char temp[PATH_MAX];
        if (::getcwd(temp, PATH_MAX) == nullptr)
            throw std::runtime_error("Internal error in getcwd(): " + std::string(strerror(errno)));
        return path(temp);
#else
        std::wstring temp(m_max_path_windows, '\0');
        if (!_wgetcwd(&temp[0], m_max_path_windows))
            throw std::runtime_error("Internal error in getcwd(): " + std::to_string(GetLastError()));
        return path(temp.c_str());
#endif
    }

#if defined(_WIN32)
    [[nodiscard]] std::wstring wstr(path_type type = native_path) const {
        std::string temp = str(type);
        int size         = MultiByteToWideChar(CP_UTF8, 0, &temp[0], static_cast<int>(temp.size()), nullptr, 0);
        std::wstring result(size, 0);
        MultiByteToWideChar(CP_UTF8, 0, &temp[0], static_cast<int>(temp.size()), &result[0], size);
        return result;
    }

    void set(const std::wstring &wstring, path_type type = native_path) {
        std::string string;
        if (!wstring.empty()) {
            int size = WideCharToMultiByte(CP_UTF8, 0, &wstring[0], static_cast<int>(wstring.size()), nullptr, 0,
                                           nullptr, nullptr);
            string.resize(size, 0);
            WideCharToMultiByte(CP_UTF8, 0, &wstring[0], static_cast<int>(wstring.size()), &string[0], size, nullptr,
                                nullptr);
        }
        set(string, type);
    }

    path &operator=(const std::wstring &str) {
        set(str);
        return *this;
    }
#endif

    bool operator==(const path &p) const { return p.m_path == m_path; }
    bool operator!=(const path &p) const { return p.m_path != m_path; }

protected:
    static std::vector<std::string> tokenize(const std::string &string, const std::string &delim) {
        std::string::size_type last_pos = 0, pos = string.find_first_of(delim, last_pos);
        std::vector<std::string> tokens;

        while (last_pos != std::string::npos) {
            if (pos != last_pos)
                tokens.push_back(string.substr(last_pos, pos - last_pos));
            last_pos = pos;
            if (last_pos == std::string::npos || last_pos + 1 == string.length())
                break;
            pos = string.find_first_of(delim, ++last_pos);
        }

        return tokens;
    }

protected:
#if defined(_WIN32)
    static constexpr size_t m_max_path_windows = 32767;
#endif
    static constexpr size_t m_max_path_windows_legacy = 260;
    path_type m_type{};
    std::vector<std::string> m_path;
    bool m_absolute{};
    bool m_smb{}; // Unused, except for on Windows
};

inline bool create_directory(const path &p) {
#if defined(_WIN32)
    return CreateDirectoryW(p.wstr().c_str(), nullptr) != 0;
#else
    return mkdir(p.str().c_str(), S_IRWXU) == 0;
#endif
}

inline bool create_directories(const path &p) {
#if defined(_WIN32)
    return SHCreateDirectory(nullptr, p.make_absolute().wstr().c_str()) == ERROR_SUCCESS;
#else
    if (create_directory(p.str().c_str()))
        return true;

    if (p.empty())
        return false;

    if (errno == ENOENT) {
        if (create_directory(p.parent_path()))
            return mkdir(p.str().c_str(), S_IRWXU) == 0;
        else
            return false;
    }
    return false;
#endif
}

NAMESPACE_END(filesystem)
