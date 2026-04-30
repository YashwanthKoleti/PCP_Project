// Sequential Bailey’s Algorithm

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <random>
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
    if (n <= 1)
        return;
    int logn = __builtin_ctz(n);
    for (int i = 0; i < n; i++)
    {
        int j = reverseBits(i, logn);
        if (i < j)
            swap(a[i], a[j]);
    }
}

void fftIterative(vector<Complex> &a)
{
    int n = a.size();
    if (n <= 1)
        return;
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

void transpose(vector<Complex> &a, int n1, int n2)
{
    vector<Complex> temp(a.size());
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            temp[j * n1 + i] = a[i * n2 + j];
        }
    }
    a = temp;
}

void baileyFFT(vector<Complex> &a, int n)
{

    int logn = __builtin_ctz(n);
    int n1 = 1 << (logn / 2);
    int n2 = n / n1;

    transpose(a, n1, n2);

    for (int i = 0; i < n1; ++i)
    {
        vector<Complex> row(n2);
        for (int j = 0; j < n2; ++j)
            row[j] = a[i * n2 + j];
        fftIterative(row);
        for (int j = 0; j < n2; ++j)
            a[i * n2 + j] = row[j];
    }

    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            double angle = -2.0 * PI * i * j / n;
            Complex twiddle(cos(angle), sin(angle));
            a[i * n2 + j] *= twiddle;
        }
    }

    transpose(a, n1, n2);

    for (int i = 0; i < n2; ++i)
    {
        vector<Complex> row(n1);
        for (int j = 0; j < n1; ++j)
            row[j] = a[i * n1 + j];
        fftIterative(row);
        for (int j = 0; j < n1; ++j)
            a[i * n1 + j] = row[j];
    }

    transpose(a, n2, n1);
}

int main()
{
    int signal_lenght;
    cin >> signal_lenght;
    vector<Complex> signal(signal_lenght);

    // for (int i = 0; i < n; i++) {
    //     double x; cin >> x;
    //     signal[i] = Complex(x, 0.0);
    // }

    mt19937 gen(42);
    uniform_real_distribution<double> dist(-100.0, 100.0);

    for (int i = 0; i < signal_lenght; i++)
    {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = chrono::high_resolution_clock::now();
    baileyFFT(signal, signal_lenght);
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