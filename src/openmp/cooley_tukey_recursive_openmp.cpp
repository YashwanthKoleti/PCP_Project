#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>
#include <omp.h>

using namespace std;

using namespace std;

using Complex = complex<double>;
const double PI = acos(-1.0);

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

    if (n > 2048)
    {

#pragma omp task shared(even)
        fftParallel(even, n / 2);

#pragma omp task shared(odd)
        fftParallel(odd, n / 2);

#pragma omp taskwait
    }
    else
    {

        fftParallel(even, n / 2);
        fftParallel(odd, n / 2);
    }

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

    // for (int i = 0; i < n; i++) {
    //     double x; cin >> x;
    //     signal[i] = Complex(x, 0.0);
    // }

    int nthreads;
    cin >> nthreads;

    omp_set_num_threads(nthreads);

    mt19937 gen(42);
    uniform_real_distribution<double> dist(-100.0, 100.0);

    for (int i = 0; i < n; i++)
    {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = chrono::high_resolution_clock::now();

#pragma omp parallel
    {

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