//
// Created by 159-mrv on 3/5/19.
//

#ifndef INFINIMESH_STOPWATCH_H
#define INFINIMESH_STOPWATCH_H

#include <allstd.h>

class Stopwatch {
public:
    Stopwatch() {
        m_seconds = 0.0;
        m_start = std::chrono::steady_clock::now();
        m_finish = m_start;
    }

    void start() {
        m_start = std::chrono::steady_clock::now();
        m_finish = m_start;
    }

    double stop() {
        m_finish = std::chrono::steady_clock::now();
        double loop_seconds = std::chrono::duration<double>(m_finish - m_start).count();
        m_seconds += loop_seconds;
        return loop_seconds;
    }

    double seconds() {
        return m_seconds;
    }

    double milliseconds() {
        return 1000.0 * m_seconds;
    }

    long milliseconds_l() const {
        return static_cast<long>(1000.0 * m_seconds);
    }

private:
    std::chrono::steady_clock::time_point m_start;
    std::chrono::steady_clock::time_point m_finish;
    double m_seconds;
};

#endif //INFINIMESH_STOPWATCH_H
