#include "roots.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <limits>

// anonymous namespace for internal helper functions and constants
namespace {
    constexpr double TOL = 1e-6;
    constexpr int MAX_ITERS = 1000000;

    bool valid_ptr(double *p) { return p != nullptr; }

    // Helper: if endpoint is already a root (within tolerance), return it
    bool check_endpoint_root(double x, double fx, double *root) {
        if (std::abs(fx) <= TOL) {
            *root = x;
            return true;
        }
        return false;
    }
}

// midpoint
bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    if (!valid_ptr(root)) {
        throw std::invalid_argument("Null pointer for root");
    }

    double fa = f(a);
    double fb = f(b);

    if (check_endpoint_root(a, fa, root) || check_endpoint_root(b, fb, root)) {
        return true;
    }

    if (fa * fb > 0) {
        return false; // No sign change
    }

    for (int iter = 0; iter < MAX_ITERS; ++iter) {
        double c = 0.5 * (a + b);
        double fc = f(c);

        if (std::abs(fc) <= TOL || (b - a) / 2 < TOL) {
            *root = c;
            return true;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    return false; // Max iterations reached
}

// midpoint 
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    if (!valid_ptr(root)) {
        throw std::invalid_argument("Null pointer for root");
    }

    double fa = f(a);
    double fb = f(b);

    if (check_endpoint_root(a, fa, root) || check_endpoint_root(b, fb, root)) {
        return true;
    }

    if (fa * fb > 0) {
        return false; // No sign change
    }

    for (int iter = 0; iter < MAX_ITERS; ++iter) {
        double c = (a * fb - b * fa) / (fb - fa);
        double fc = f(c);

        if (std::abs(fc) <= TOL) {
            *root = c;
            return true;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    return false; // Max iterations reached
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
    if (!valid_ptr(root)) {
        throw std::invalid_argument("Null pointer for root");
    }

    for (int iter = 0; iter < MAX_ITERS; ++iter) {
        double fc = f(c);
        double gc = g(c);

        if (std::abs(fc) <= TOL) {
            *root = c;
            return true;
        }

        if (gc == 0) {
            return false; // Derivative zero
        }

        double c_next = c - fc / gc;

        if (c_next < a || c_next > b) {
            return false; // Out of bounds
        }

        c = c_next;
    }

    return false; // Max iterations reached
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    if (!valid_ptr(root)) {
        throw std::invalid_argument("Null pointer for root");
    }

    double fc = f(c);
    double fb = f(b);

    if (check_endpoint_root(c, fc, root) || check_endpoint_root(b, fb, root)) {
        return true;
    }

    double a_orig = a;
    double b_orig = b;

    for (int iter = 0; iter < MAX_ITERS; ++iter) {
        if (fc - fb == 0) {
            return false; // Prevent division by zero
        }

        double c_next = c - fc * (c - b) / (fc - fb);
        double fc_next = f(c_next);

        if (std::abs(fc_next) <= TOL) {
            *root = c_next;
            return true;
        }

        // Check bounds against original interval
        if (c_next < a_orig || c_next > b_orig) {
            return false; // Out of bounds
        }

        b = c;
        fb = fc;
        c = c_next;
        fc = fc_next;
    }

    return false; // Max iterations reached
}