//
// Created by Joseph Zhang on 7/25/25.
//

#ifndef SEMICLASSICAL_RADIATION_TIMER_LOGGER_H
#define SEMICLASSICAL_RADIATION_TIMER_LOGGER_H

#pragma once
#include <chrono>
#include <string>
#include <cstdio>

class TimerLogger {
public:
    using clock_t = std::chrono::steady_clock;

    TimerLogger(int total_steps)
            : start_time(clock_t::now()), total_steps(total_steps), step_count(0) {}

    void nextStep() { ++step_count; }

    // 返回"hh:mm:ss:cs"格式字符串
    static std::string format_time(double t) {
        int h = static_cast<int>(t / 3600);
        int m = static_cast<int>(t / 60) % 60;
        int s = static_cast<int>(t) % 60;
        int cs = static_cast<int>((t - static_cast<int>(t)) * 100); // centiseconds
        char buf[32];
        snprintf(buf, sizeof(buf), "%d:%02d:%02d:%02d", h, m, s, cs);
        return std::string(buf);
    }

    // 获取当前已用时间（秒）
    double elapsed_s() const {
        auto t_now = clock_t::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(t_now - start_time).count() / 1000.0;
    }

    // 获取剩余时间（秒）
    double eta_s() const {
        if (step_count == 0) return 0.0;
        double avg = elapsed_s() / step_count;
        int left = total_steps - step_count;
        return avg * left;
    }

    // 返回带格式的完整log字符串
    std::string log(double tNext_tp0) const {
        return "[LOG] steps=" + std::to_string(step_count) +
               " t/tp0=" + std::to_string(tNext_tp0) +
               ", elapsed " + format_time(elapsed_s()) +
               ", ETA " + format_time(eta_s());
    }

    int step() const { return step_count; }

private:
    clock_t::time_point start_time;
    int total_steps;
    int step_count;
};

#endif //SEMICLASSICAL_RADIATION_TIMER_LOGGER_H
