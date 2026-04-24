#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

// Define a type alias for complex numbers for cleaner code
using Complex = std::complex<double>;

// Define PI
const double PI = std::acos(-1.0);

// Recursive Cooley-Tukey Radix-2 FFT algorithm
void fft(std::vector<Complex>& a) {
    int n = a.size();
    
    // Base case: if the array has 1 or 0 elements, it's already sorted
    if (n <= 1) {
        return;
    }

    // Divide the array into even and odd indexed elements
    std::vector<Complex> even(n / 2);
    std::vector<Complex> odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    // Conquer: recursively compute the FFT of both halves
    fft(even);
    fft(odd);

    // Combine: merge the results back into the original array
    for (int k = 0; k < n / 2; ++k) {
        // Calculate the "twiddle factor"
        Complex t = std::polar(1.0, -2.0 * PI * k / n) * odd[k];
        
        // Butterfly operation
        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
}

int main() {
    // Example input: A simple signal (size must be a power of 2)
    std::vector<Complex> signal = {
        {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0},
        {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}
    };

    std::cout << "Original Signal (Time Domain):" << std::endl;
    for (const auto& val : signal) {
        std::cout << val << " ";
    }
    std::cout << "\n\n";

    // Perform FFT
    fft(signal);

    std::cout << "FFT Result (Frequency Domain):" << std::endl;
    for (const auto& val : signal) {
        // Cleaning up output to handle very small floating point inaccuracies
        double real = (std::abs(val.real()) < 1e-10) ? 0.0 : val.real();
        double imag = (std::abs(val.imag()) < 1e-10) ? 0.0 : val.imag();
        
        std::cout << "(" << std::fixed << std::setprecision(4) << real 
                  << ", " << imag << ")\n";
    }

    return 0;
}