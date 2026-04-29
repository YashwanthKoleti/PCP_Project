#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>
#include <omp.h> // Include OpenMP

using namespace std;

using Complex = complex<double>;
const double PI = acos(-1.0);

// We no longer need to pass 'available_threads' because OpenMP handles the thread pool.
void fftParallel(vector<Complex> &a, int n)
{
    if (n <= 1)
    {
        return;
    }

    vector<Complex> even(n / 2);
    vector<Complex> odd(n / 2);
    for (int i = 0; i < n / 2; ++i)
    {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    // Using your optimal threshold to avoid task creation overhead on small sizes
    if (n > 32768)
    {
        // Tell OpenMP to create a task for this function call. 
        // 'shared' ensures the task modifies the correct vectors in memory.
        #pragma omp task shared(even)
        fftParallel(even, n / 2);

        #pragma omp task shared(odd)
        fftParallel(odd, n / 2);

        // This is crucial: wait for both of the above tasks to finish 
        // before proceeding to the combination step.
        #pragma omp taskwait
    }
    else
    {
        // Serial execution below the threshold
        fftParallel(even, n / 2);
        fftParallel(odd, n / 2);
    }       

    // Combine
    for (int k = 0; k < n / 2; ++k)
    {
        Complex t = polar(1.0, -2.0 * PI * k / n) * odd[k];
        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
}

int main()
{
    int n;
    cin >> n;
    vector<Complex> signal(n);

    int nthreads;
    cin >> nthreads;
    
    // Configure OpenMP to use the requested number of threads
    omp_set_num_threads(nthreads);

    mt19937 gen(42);
    uniform_real_distribution<double> dist(-100.0, 100.0);

    for (int i = 0; i < n; i++)
    {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = chrono::high_resolution_clock::now();
    
    // To use Tasks in OpenMP, we must first spin up the thread pool.
    // The 'parallel' block wakes up the threads.
    #pragma omp parallel
    {
        // The 'single' block ensures only ONE thread starts the root of the recursion.
        // That single thread will generate tasks, waking up the others to help.
        #pragma omp single nowait
        {
            fftParallel(signal, n);
        }
    }
    
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << elapsed.count() << " seconds\n";

    return 0;
}