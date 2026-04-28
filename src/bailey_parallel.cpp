#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <thread>
#include <algorithm>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

// Helper to reverse bits for the standard iterative FFT
int reverseBits(int x, int logn) {
    int result = 0;
    for (int i = 0; i < logn; i++) {
        if (x & (1 << i)) {
            result |= 1 << (logn - 1 - i);
        }
    }
    return result;
}

void bitReverse(std::vector<Complex> &a) {
    int n = a.size();
    if (n <= 1) return;
    int logn = __builtin_ctz(n);
    for (int i = 0; i < n; i++) {
        int j = reverseBits(i, logn);
        if (i < j) std::swap(a[i], a[j]);
    }
}

// Standard Iterative FFT used as a building block for rows
void fftIterative(std::vector<Complex> &a) {
    int n = a.size();
    if (n <= 1) return;
    bitReverse(a);
    for (int len = 2; len <= n; len <<= 1) {
        double angle = -2.0 * PI / len;
        Complex wlen(std::cos(angle), std::sin(angle));
        for (int i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for (int j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// Parallel Matrix Transpose
void transposeParallel(std::vector<Complex> &a, int n1, int n2, int nthreads) {
    std::vector<Complex> temp(a.size());
    auto worker = [&](int start_row, int end_row) {
        for (int i = start_row; i < end_row; ++i) {
            for (int j = 0; j < n2; ++j) {
                temp[j * n1 + i] = a[i * n2 + j];
            }
        }
    };

    std::vector<std::thread> threads;
    int rows_per_thread = n1 / nthreads;
    for (int t = 0; t < nthreads; ++t) {
        int start = t * rows_per_thread;
        int end = (t == nthreads - 1) ? n1 : (t + 1) * rows_per_thread;
        threads.emplace_back(worker, start, end);
    }
    for (auto &th : threads) th.join();
    a = std::move(temp);
}

void baileyFFTParallel(std::vector<Complex> &a, int n, int nthreads) {
    int logn = __builtin_ctz(n);
    int n1 = 1 << (logn / 2);
    int n2 = n / n1;

    // Step 1: Parallel Transpose
    transposeParallel(a, n1, n2, nthreads);

    // Step 2: Parallel Row FFTs (length n2)
    auto rowWorker1 = [&](int start_row, int end_row) {
        for (int i = start_row; i < end_row; ++i) {
            std::vector<Complex> row(n2);
            for (int j = 0; j < n2; ++j) row[j] = a[i * n2 + j];
            fftIterative(row);
            for (int j = 0; j < n2; ++j) a[i * n2 + j] = row[j];
        }
    };

    std::vector<std::thread> threads;
    int rpt1 = n1 / nthreads;
    for (int t = 0; t < nthreads; ++t) {
        threads.emplace_back(rowWorker1, t * rpt1, (t == nthreads - 1) ? n1 : (t + 1) * rpt1);
    }
    for (auto &th : threads) th.join();
    threads.clear();

    // Step 3: Parallel Twiddle Factor Multiplication
    auto twiddleWorker = [&](int start_row, int end_row) {
        for (int i = start_row; i < end_row; ++i) {
            for (int j = 0; j < n2; ++j) {
                double angle = -2.0 * PI * i * j / n;
                a[i * n2 + j] *= Complex(std::cos(angle), std::sin(angle));
            }
        }
    };

    for (int t = 0; t < nthreads; ++t) {
        threads.emplace_back(twiddleWorker, t * rpt1, (t == nthreads - 1) ? n1 : (t + 1) * rpt1);
    }
    for (auto &th : threads) th.join();
    threads.clear();

    // Step 4: Parallel Transpose
    transposeParallel(a, n1, n2, nthreads);

    // Step 5: Parallel Row FFTs (length n1)
    auto rowWorker2 = [&](int start_row, int end_row) {
        for (int i = start_row; i < end_row; ++i) {
            std::vector<Complex> row(n1);
            for (int j = 0; j < n1; ++j) row[j] = a[i * n1 + j];
            fftIterative(row);
            for (int j = 0; j < n1; ++j) a[i * n1 + j] = row[j];
        }
    };

    int rpt2 = n2 / nthreads;
    for (int t = 0; t < nthreads; ++t) {
        threads.emplace_back(rowWorker2, t * rpt2, (t == nthreads - 1) ? n2 : (t + 1) * rpt2);
    }
    for (auto &th : threads) th.join();
    threads.clear();

    // Step 6: Final Parallel Transpose
    transposeParallel(a, n2, n1, nthreads);
}

int main() {
    int n, nthreads;
    if (!(std::cin >> n >> nthreads)) return 0;

    std::vector<Complex> signal(n);
    for (int i = 0; i < n; i++) {
        double x;
        std::cin >> x;
        signal[i] = Complex(x, 0.0);
    }

    auto start = std::chrono::high_resolution_clock::now();
    baileyFFTParallel(signal, n, nthreads);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time: " << elapsed.count() << " seconds\n";

    return 0;
}