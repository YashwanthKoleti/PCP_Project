//  Iterative Cooley Tukey

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

int reverseBits(int x, int logn)
{
    int result = 0;
    for (int i = 0; i < logn; i++)
    {
        if (x & (1 << i))
        {
            result |= 1 << (logn - 1 - i);
        }
    }
    return result;
}

void bitReverse(vector<Complex> &a)
{
    int n = a.size();
    int logn = __builtin_ctz(n);

    for (int i = 0; i < n; i++)
    {
        int j = reverseBits(i, logn);

        if (i < j)
        {
            swap(a[i], a[j]);
        }
    }
}

void fftIterative(vector<Complex> &a)
{
    int n = a.size();

    bitReverse(a);

    for (int len = 2; len <= n; len <<= 1)
    {

        double angle = -2.0 * PI / len;
        Complex wlen(cos(angle), sin(angle));

        for (int i = 0; i < n; i += len)
        {
            Complex w(1.0, 0.0);

            for (int j = 0; j < len / 2; ++j)
            {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;

                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;

                w *= wlen;
            }
        }
    }
}

int main()
{

   int signal_lenght;
    cin>>signal_lenght;

    // for (int i = 0; i < n; i++) {
    //     double x; cin >> x;
    //     signal[i] = Complex(x, 0.0);
    // }

    vector<Complex> signal(signal_lenght);
    mt19937 gen(42); 
    uniform_real_distribution<double> dist(-100.0, 100.0);
    
    for (int i = 0; i < signal_lenght; i++) {
        signal[i] = Complex(dist(gen), 0.0);
    }

    cout << "\n\n";

    auto start = chrono::high_resolution_clock::now();
    fftIterative(signal);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << elapsed.count() << " seconds\n";

    // for (const auto &val : signal)
    // {

    //     double real = (abs(val.real()) < 1e-10) ? 0.0 : val.real();
    //     double imag = (abs(val.imag()) < 1e-10) ? 0.0 : val.imag();

    //     cout << "(" << fixed << setprecision(4) << real
    //               << ", " << imag << ")\n";
    // }

    return 0;
}