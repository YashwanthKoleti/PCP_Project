// Recursive Cooley Tukey

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <random>
#include <chrono>
#include <iomanip>
#include <cmath>

using namespace std;

using Complex = complex<double>;

const double PI = acos(-1.0);

void fft(vector<Complex> &a)
{
    int n = a.size();

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

    fft(even);
    fft(odd);

    for (int k = 0; k < n / 2; ++k)
    {

        Complex t = polar(1.0, -2.0 * PI * k / n) * odd[k];

        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
}

int main() {
    int signal_lenght;
    cin>>signal_lenght;
    vector<Complex> signal(signal_lenght);

    // for (int i = 0; i < n; i++) {
    //     double x; cin >> x;
    //     signal[i] = Complex(x, 0.0);
    // }

    mt19937 gen(42); 
    uniform_real_distribution<double> dist(-100.0, 100.0);
    
    for (int i = 0; i < signal_lenght; i++) {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = chrono::high_resolution_clock::now();
    fft(signal);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << elapsed.count() << " seconds\n";

    // for (const auto& val : signal) {

    //     double real = (abs(val.real()) < 1e-10) ? 0.0 : val.real();
    //     double imag = (abs(val.imag()) < 1e-10) ? 0.0 : val.imag();

    //     cout << "(" << fixed << setprecision(4) << real
    //               << ", " << imag << ")\n";
    // }

    return 0;
}