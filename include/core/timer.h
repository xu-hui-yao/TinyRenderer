#pragma once

#include <chrono>
#include <core/common.h>

M_NAMESPACE_BEGIN
/**
 * \brief Simple timer with millisecond precision
 *
 * This class is convenient for collecting performance data
 */
class Timer {
public:
    /// Create a new timer and reset it
    Timer() { reset(); }

    /// Reset the timer to the current time
    void reset() { start = std::chrono::system_clock::now(); }

    /// Return the number of milliseconds elapsed since the timer was last reset
    [[nodiscard]] double elapsed() const {
        auto now = std::chrono::system_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
        return static_cast<double>(duration.count());
    }

    /// Like \ref elapsed(), but return a human-readable string
    [[nodiscard]] std::string elapsed_string(bool precise = false) const {
        return time_string(elapsed(), precise);
    }

    /// Return the number of milliseconds elapsed since the timer was last reset
    /// and then reset it
    double lap() {
        auto now = std::chrono::system_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(now - start);
        start = now;
        return static_cast<double>(duration.count());
    }

    /// Like \ref lap(), but return a human-readable string
    std::string lap_string(bool precise = false) {
        return time_string(lap(), precise);
    }

private:
    std::chrono::system_clock::time_point start;
};

M_NAMESPACE_END
