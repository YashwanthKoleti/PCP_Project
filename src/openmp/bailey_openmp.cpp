#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <omp.h> // Required for OpenMP

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

// Helper to reverse bits
int reverseBits(int x, int logn) {
    int result = 0;
    for (int i = 0; i < logn; i++) {
        if (x & (1 << i)) {
            result |= 1 << (logn - 1 - i);
        }
    }
    return result;
}

// Bit reverse now operates on raw pointers for zero-allocation performance
void bitReverse(Complex *a, int n) {
    if (n <= 1) return;
    int logn = __builtin_ctz(n);
    for (int i = 0; i < n; i++) {
        int j = reverseBits(i, logn);
        if (i < j) std::swap(a[i], a[j]);
    }
}

// FFT operates directly on raw memory
void fftIterative(Complex *a, int n) {
    if (n <= 1) return;
    bitReverse(a, n);
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

// Cache-Oblivious Blocked Parallel Transpose using OpenMP
void transposeParallelBlocked(const Complex *src, Complex *dst, int n1, int n2) {
    const int BLOCK_SIZE = 32; 
    
    // OpenMP 'collapse(2)' merges the two outer loops into one flat parallel workload
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n1; i += BLOCK_SIZE) {
        for (int j = 0; j < n2; j += BLOCK_SIZE) {
            int max_i = std::min(i + BLOCK_SIZE, n1);
            int max_j = std::min(j + BLOCK_SIZE, n2);
            for (int ii = i; ii < max_i; ++ii) {
                for (int jj = j; jj < max_j; ++jj) {
                    dst[jj * n1 + ii] = src[ii * n2 + jj];
                }
            }
        }
    }
}

void baileyFFTParallel(std::vector<Complex> &a, int n) {
    int logn = __builtin_ctz(n);
    int n1 = 1 << (logn / 2);
    int n2 = n / n1;

    // Allocate ONE buffer for the entire algorithm
    std::vector<Complex> buffer(n);

    // Step 1: Parallel Blocked Transpose (a -> buffer)
    transposeParallelBlocked(a.data(), buffer.data(), n1, n2);

    // Step 2: Parallel Row FFTs on buffer
    #pragma omp parallel for
    for (int i = 0; i < n1; ++i) {
        fftIterative(&buffer[i * n2], n2);
    }

    // Step 3: Parallel Twiddle Factor Multiplication
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            double angle = -2.0 * PI * i * j / n;
            buffer[i * n2 + j] *= Complex(std::cos(angle), std::sin(angle));
        }
    }

    // Step 4: Parallel Blocked Transpose (buffer -> a)
    transposeParallelBlocked(buffer.data(), a.data(), n1, n2);

    // Step 5: Parallel Row FFTs on a
    #pragma omp parallel for
    for (int i = 0; i < n2; ++i) {
        fftIterative(&a[i * n1], n1);
    }

    // Step 6: Final Parallel Transpose (a -> buffer)
    transposeParallelBlocked(a.data(), buffer.data(), n2, n1);
    
    // Move final result back to a
    a = std::move(buffer); 
}

int main() {
    int n, nthreads;
    std::cin >> n >> nthreads;

    // Set the number of threads for OpenMP globally
    omp_set_num_threads(nthreads);

    std::vector<Complex> signal(n);
    for (int i = 0; i < n; i++) {
        double x;
        std::cin >> x;
        signal[i] = Complex(x, 0.0);
    }

    auto start = std::chrono::high_resolution_clock::now();
    baileyFFTParallel(signal, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds\n";

    return 0;
}