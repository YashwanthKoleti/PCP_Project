#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

// Helper function to rearrange elements in bit-reversed order
void bitReverse(std::vector<Complex>& a) {
    int n = a.size();
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }
}

// Generalized Iterative Radix-2 FFT (handles both Forward and Inverse)
void fftRadix2(std::vector<Complex>& a, bool invert) {
    int n = a.size();
    bitReverse(a);

    for (int len = 2; len <= n; len <<= 1) {
        // Angle changes sign depending on forward vs inverse transform
        double angle = (invert ? 2.0 : -2.0) * PI / len;
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

    // Scale the output if we are doing an inverse FFT
    if (invert) {
        for (Complex& x : a) {
            x /= n;
        }
    }
}

// Bluestein's Algorithm for arbitrary size N
void fftBluestein(std::vector<Complex>& a) {
    int N = a.size();
    if (N <= 1) return;

    // Step 1: Find the next power of 2 that is >= 2N - 1 for zero-padding
    int M = 1;
    while (M < 2 * N - 1) {
        M <<= 1;
    }

    std::vector<Complex> A(M, {0.0, 0.0});
    std::vector<Complex> B(M, {0.0, 0.0});

    // Step 2: Initialize sequences A and B based on the chirp sequence
    for (int n = 0; n < N; ++n) {
        double theta = -PI * (long long)n * n / N;
        Complex w_n(std::cos(theta), std::sin(theta));
        
        // Sequence A: Input multiplied by chirp
        A[n] = a[n] * w_n;

        // Sequence B: Inverse chirp, arranged for circular convolution
        Complex w_n_inv(std::cos(-theta), std::sin(-theta));
        B[n] = w_n_inv;
        if (n > 0) {
            B[M - n] = w_n_inv; // Wrap around for negative indices
        }
    }

    // Step 3: Perform Convolution using Radix-2 FFT
    fftRadix2(A, false); // Forward FFT of A
    fftRadix2(B, false); // Forward FFT of B

    // Point-wise multiplication in the frequency domain
    std::vector<Complex> C(M);
    for (int i = 0; i < M; ++i) {
        C[i] = A[i] * B[i];
    }

    fftRadix2(C, true); // Inverse FFT to get back to the time domain

    // Step 4: Multiply by the chirp sequence again to get the final output
    for (int k = 0; k < N; ++k) {
        double theta = -PI * (long long)k * k / N;
        Complex w_k(std::cos(theta), std::sin(theta));
        a[k] = C[k] * w_k;
    }
}

int main() {
    // Example input: Size 5 (NOT a power of 2)
    std::vector<Complex> signal = {
        {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}
    };

    std::cout << "Original Signal (N=" << signal.size() << "):\n";
    for (const auto& val : signal) {
        std::cout << val.real() << " ";
    }
    std::cout << "\n\n";

    // Perform Bluestein's FFT
    fftBluestein(signal);

    std::cout << "FFT Result:\n";
    for (const auto& val : signal) {
        double real = (std::abs(val.real()) < 1e-10) ? 0.0 : val.real();
        double imag = (std::abs(val.imag()) < 1e-10) ? 0.0 : val.imag();
        
        std::cout << "(" << std::fixed << std::setprecision(4) << real 
                  << ", " << imag << ")\n";
    }

    return 0;
}