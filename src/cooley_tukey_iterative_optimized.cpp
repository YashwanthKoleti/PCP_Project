#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

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

struct ThreadData {
    int i;
    int n;
};

void* worker(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    int i = data->i;
    int n = data->n;

    double angle = -2.0 * PI * i / n;
    roots[i] = Complex(std::cos(angle), std::sin(angle));

    return nullptr;
}

std::vector<Complex> roots;
void precompute(int n)
{
    std::vector<pthread_t> threads(n);
    std::vector<ThreadData> args(n);

    for (int i = 0; i < n; i++) {
        args[i] = {i, n};
        pthread_create(&threads[i], nullptr, worker, &args[i]);
    }

    for (int i = 0; i < n; i++) {
        pthread_join(threads[i], nullptr);
    }
}

void opti(std::vector<Complex>& a,int i, int len, int step)
{
    // Perform the butterfly operation on the two halves of the chunk
            for (int j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * roots[j*step];
                
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
            }
}

// Iterative Cooley-Tukey Radix-2 FFT algorithm
void fftIterative(std::vector<Complex>& a) {
    int n = a.size();


    // Precompute all roots of unity

    // Step 1: Rearrange array elements based on bit-reversed indices
    bitReverse(a);

    // Step 2: Bottom-up butterfly operations
    // 'len' represents the size of the sub-problem we are currently solving (2, 4, 8, 16...)
    for (int len = 2; len <= n; len <<= 1) {
        int step = n / len;

        // Iterate through the array in chunks of size 'len'
        for (int i = 0; i < n; i += len) {
            
            // Perform the butterfly operation on the two halves of the chunk
            for (int j = 0; j < len / 2; ++j) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * roots[j*step];
                
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
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
    int n = signal.size();

    // Perform Iterative FFT
    roots.resize(n);
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