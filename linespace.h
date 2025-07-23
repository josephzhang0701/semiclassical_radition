//
// Created by Joseph Zhang on 7/23/25.
//

#ifndef SEMICLASSICAL_RADIATION_LINESPACE_H
#define SEMICLASSICAL_RADIATION_LINESPACE_H

#pragma once
#include <vector>
#include <stdexcept>

/**
 * @class Linspace
 * @brief Generate a sequence of evenly‐spaced values between start and end.
 */

class linespace {
public:
    /**
     * @param start  The first value in the sequence.
     * @param end    The last value in the sequence.
     * @param num    Number of points (must be >= 2).
     * @throws std::invalid_argument if num < 2.
     */
    linespace(double start, double end, int num);

    /**
     * @brief Compute and return the vector of points.
     */
    std::vector<double> generate() const;

private:
    double start_;
    double end_;
    int num_;

};


#endif //SEMICLASSICAL_RADIATION_LINESPACE_H
