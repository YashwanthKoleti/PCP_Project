#include <iostream>
#include <chrono>
#include <random>
#include <thread>
#include <fftw3.h>

using namespace std;

int main()
{
    int n;
    cin >> n;
    
    int nthreads;
    cin >> nthreads;

    

    // 1. Initialize FFTW threads before any FFTW calls
    if (fftw_init_threads() == 0) {
        cerr << "Error: FFTW thread initialization failed.\n";
        return 1;
    }

    int num_threads = nthreads > 0 ? nthreads : 2;
    fftw_plan_with_nthreads(num_threads);

    // 2. Allocate aligned memory for maximum SIMD performance
    fftw_complex* signal = fftw_alloc_complex(n);

    // 3. CREATE THE PLAN FIRST
    // Using FFTW_ESTIMATE makes the planning instant. 
    // If you plan to run the same FFT thousands of times in a loop, 
    // change this to FFTW_MEASURE, but ALWAYS do it before filling the data!
    fftw_plan plan = fftw_plan_dft_1d(
        n, signal, signal, FFTW_FORWARD, FFTW_ESTIMATE 
    );

    // 4. FILL THE ARRAY AFTER PLANNING
    // This guarantees the planner hasn't overwritten your data with garbage.
    mt19937 gen(42);
    uniform_real_distribution<double> dist(-100.0, 100.0);

    for (int i = 0; i < n; i++)
    {
        signal[i][0] = dist(gen); // Real
        signal[i][1] = 0.0;       // Imaginary
    }

    cout << "\n\n";

    // 5. EXECUTE AND MEASURE
    auto start = chrono::high_resolution_clock::now();
    
    fftw_execute(plan);
    
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << elapsed.count() << " seconds\n";
    cout << "Threads used: " << num_threads << "\n\n";

    // 6. CLEANUP
    fftw_destroy_plan(plan);
    fftw_cleanup_threads();
    fftw_cleanup();
    fftw_free(signal); 

    return 0;
}