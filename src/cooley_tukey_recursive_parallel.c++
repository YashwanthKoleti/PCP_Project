#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <thread>
#include <iomanip>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

void fftParallel(std::vector<Complex>& a, int available_threads) {
    int n = a.size();
    
    if (n <= 1) {
        return;
    }

    std::vector<Complex> even(n / 2);
    std::vector<Complex> odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    if (available_threads >= 2) {

        int threads_for_even = available_threads / 2;
        int threads_for_odd = available_threads - threads_for_even;

        std::thread even_worker(fftParallel, std::ref(even), threads_for_even);

        fftParallel(odd, threads_for_odd);
        even_worker.join();
    } else {
        
        fftParallel(even, 1);
        fftParallel(odd, 1);
    }

    // Combine
    for (int k = 0; k < n / 2; ++k) {
        Complex t = std::polar(1.0, -2.0 * PI * k / n) * odd[k];
        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
}

int main() {
    // Array size MUST be a power of 2
    // std::vector<Complex> signal = {
    //     {1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0},
    //     {5.0, 0.0}, {6.0, 0.0}, {7.0, 0.0}, {8.0, 0.0}
    // };

    int signal_lenght;
    std::cin>>signal_lenght;
    std::vector<Complex> signal(signal_lenght);
    for(int i = 0;i < signal_lenght;i++)
    {
        double x;
        std::cin >> x;
        signal[i] = Complex(x, 0.0);
    }

    // Number of total threads to use (given as input)
    // Try changing this to 1, 2, 4, etc.
    int N_THREADS;
    std :: cin >> N_THREADS;

    std::cout << "Original Signal (Time Domain):\n";
    for (const auto& val : signal) {
        std::cout << val.real() << " ";
    }
    std::cout << "\n\n";

    std::cout << "Running Parallel FFT with " << N_THREADS << " thread(s)...\n\n";

    // Perform Parallel FFT
    auto start = std::chrono::high_resolution_clock::now();
    fftParallel(signal, N_THREADS);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time: " << elapsed.count() << " seconds\n";

    std::cout << "FFT Result (Frequency Domain):\n";
    for (const auto& val : signal) {
        // Clean up output for floating point inaccuracies
        double real = (std::abs(val.real()) < 1e-10) ? 0.0 : val.real();
        double imag = (std::abs(val.imag()) < 1e-10) ? 0.0 : val.imag();
        
        std::cout << "(" << std::fixed << std::setprecision(4) << real 
                  << ", " << imag << ")\n";
    }

    return 0;
}