#include <omp.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

struct StageRoots {
    int offset;  
    int count;   
};

std::vector<Complex> all_roots;
std::vector<StageRoots> stage_info;

void precomputeAllRoots(int n) {
    all_roots.resize(n);
    stage_info.clear();
    
    int offset = 0;
    for (int len = 2; len <= n; len <<= 1) {
        int half = len / 2;
        StageRoots sr;
        sr.offset = offset;
        sr.count = half;
        for (int j = 0; j < half; j++) {
            double angle = -2.0 * PI * j / len;
            all_roots[offset + j] = Complex(std::cos(angle), std::sin(angle));
        }
        stage_info.push_back(sr);
        offset += half;
    }
}

void bitReverse(std::vector<Complex> &a) {
    int n = a.size();
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        while (j & bit) { j ^= bit; bit >>= 1; }
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

void fftIterative(std::vector<Complex> &a, int NUM_THREADS) {
    int n = a.size();
    
    if ((n & (n - 1)) != 0) {
        std::cerr << "Error: Signal size N must be a power of 2." << std::endl;
        return;
    }

    bitReverse(a);
    precomputeAllRoots(n);

    int num_stages = stage_info.size();
    Complex* a_ptr = a.data(); // Get raw pointer for fast memory access

    // OpenMP automatically creates the thread pool and handles the joins at the end.
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        // These variables are inherently private to each thread because 
        // they are declared inside the parallel region.
        int tid = omp_get_thread_num();
        int T = omp_get_num_threads();
        
        const int total_butterflies = n / 2;  
        int bpt = total_butterflies / T;
        int rem = total_butterflies % T;
        int b_start, b_end;
        
        // Static workload distribution (exactly as you designed it)
        if (tid < rem) {
            b_start = tid * (bpt + 1);
            b_end   = b_start + bpt + 1;
        } else {
            b_start = tid * bpt + rem;
            b_end   = b_start + bpt;
        }

        for (int s = 0; s < num_stages; s++) {
            const int half = stage_info[s].count;       
            const int len  = half << 1;
            const int roff = stage_info[s].offset;      

            // OPTIMIZATION 1: Calculate 2D indices outside the loop
            int g = b_start / half;
            int j = b_start % half;

            for (int b = b_start; b < b_end; b++) {
                int base = g * len;
                int idx = base + j;
                
                Complex u = a_ptr[idx];
                Complex v = a_ptr[idx + half] * all_roots[roff + j];
                a_ptr[idx]        = u + v;
                a_ptr[idx + half] = u - v;

                // Increment local coordinates directly
                j++;
                if (j == half) {
                    j = 0;
                    g++;
                }
            }

            // Native OpenMP barrier replaces the custom Apple pthread_barrier
            #pragma omp barrier
        }
    }
}

int main() {
    int n;
    if (!(std::cin >> n)) return 0;
    
    int nthreads;
    if (!(std::cin >> nthreads)) return 0;

    std::vector<Complex> signal(n);
    std::mt19937 gen(42); 
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    
    for (int i = 0; i < n; i++) {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = std::chrono::high_resolution_clock::now();
    fftIterative(signal, nthreads);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds\n";

    return 0;
}