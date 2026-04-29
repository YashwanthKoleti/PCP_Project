#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <thread>
#include <chrono>
#include <random>
#include <iomanip>
using namespace std;

using Complex = complex<double>;
const double PI = acos(-1.0);

void fftParallel(vector<Complex> &a, int n, int available_threads)
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

    if (available_threads >= 2 && n > 32768)
    {

        int threads_for_even = available_threads / 2;
        int threads_for_odd = available_threads - threads_for_even;

        thread even_worker(fftParallel, ref(even), n / 2, threads_for_even);

        fftParallel(odd, n / 2, threads_for_odd);
        even_worker.join();
    }
    else
    {

        fftParallel(even, n / 2, 1);
        fftParallel(odd, n / 2, 1);
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
    mt19937 gen(42);
    uniform_real_distribution<double> dist(-100.0, 100.0);

    for (int i = 0; i < n; i++)
    {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = chrono::high_resolution_clock::now();
    fftParallel(signal, n, nthreads);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << elapsed.count() << " seconds\n";

    // for (const auto& val : signal) {
    //     // Clean up output for floating point inaccuracies
    //     double real = (abs(val.real()) < 1e-10) ? 0.0 : val.real();
    //     double imag = (abs(val.imag()) < 1e-10) ? 0.0 : val.imag();

    //     cout << "(" << fixed << setprecision(4) << real
    //               << ", " << imag << ")\n";
    // }

    return 0;
}
