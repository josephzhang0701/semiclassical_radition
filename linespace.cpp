//
// Created by Joseph Zhang on 7/23/25.
//
#include "linespace.h"

linespace::linespace(double start, double end, int num)
        : start_(start), end_(end), num_(num)
{
    if (num_ < 2) {
        throw std::invalid_argument("Linespace requires at least two points");
    }
}

std::vector<double> linespace::generate() const {
    int segments = num_ - 1;
    double step     = (end_ - start_) / segments;
    std::vector<double> pts;
    pts.reserve(num_);

    // push start, intermediate, then end
    for (int i = 0; i <= segments; ++i) {
        pts.push_back(start_ + i * step);
    }
    return pts;
}