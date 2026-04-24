#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

// Helper function to rearrange elements in bit-reversed order
int reverseBits(int x, int logn) {
    int result = 0;
    for (int i = 0; i < logn; i++) {
        if (x & (1 << i)) {
            result |= 1 << (logn - 1 - i);
        }
    }
    return result;
}

void bitReverse(std::vector<Complex>& a) {
    int n = a.size();
    int logn = __builtin_ctz(n); // since n is power of 2

    for (int i = 0; i < n; i++) {
        int j = reverseBits(i, logn);

        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }
}

// Iterative Cooley-Tukey Radix-2 FFT algorithm
void fftIterative(std::vector<Complex>& a) {
    int n = a.size();

    // Step 1: Rearrange array elements based on bit-reversed indices
    bitReverse(a);

    // Step 2: Bottom-up butterfly operations
    // 'len' represents the size of the sub-problem we are currently solving (2, 4, 8, 16...)
    for (int len = 2; len <= n; len <<= 1) {
        // Calculate the base twiddle factor for this sub-problem size
        double angle = -2.0 * PI / len;
        Complex wlen(std::cos(angle), std::sin(angle));

        // Iterate through the array in chunks of size 'len'
        for (int i = 0; i < n; i += len) {
            Complex w(1.0, 0.0); // Start with w^0 = 1
            
            // Perform the butterfly operation on the two halves of the chunk
            for (int j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                
                w *= wlen; // Increment the twiddle factor for the next element
            }
        }
    }
}

int main() {
    // Example input: Signal size MUST be a power of 2
    std::vector<Complex> signal = {
        {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0},
        {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}
    };

    std::cout << "Original Signal (Time Domain):\n";
    for (const auto& val : signal) {
        std::cout << val.real() << " ";
    }
    std::cout << "\n\n";

    // Perform Iterative FFT
    fftIterative(signal);

    std::cout << "FFT Result (Frequency Domain):\n";
    for (const auto& val : signal) {
        // Cleaning up output to handle very small floating point inaccuracies
        double real = (std::abs(val.real()) < 1e-10) ? 0.0 : val.real();
        double imag = (std::abs(val.imag()) < 1e-10) ? 0.0 : val.imag();
        
        std::cout << "(" << std::fixed << std::setprecision(4) << real 
                  << ", " << imag << ")\n";
    }

    return 0;
}