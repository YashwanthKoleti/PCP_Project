#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <omp.h>

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

void bitReverse(Complex *a, int n)
{
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

void fftIterative(Complex *a, int n)
{
    if (n <= 1)
        return;
    bitReverse(a, n);
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

void transposeParallelBlocked(const Complex *src, Complex *dst, int n1, int n2)
{
    const int BLOCK_SIZE = 32;

#pragma omp parallel for collapse(2)
    for (int i = 0; i < n1; i += BLOCK_SIZE)
    {
        for (int j = 0; j < n2; j += BLOCK_SIZE)
        {
            int max_i = min(i + BLOCK_SIZE, n1);
            int max_j = min(j + BLOCK_SIZE, n2);
            for (int ii = i; ii < max_i; ++ii)
            {
                for (int jj = j; jj < max_j; ++jj)
                {
                    dst[jj * n1 + ii] = src[ii * n2 + jj];
                }
            }
        }
    }
}

void baileyFFTParallel(vector<Complex> &a, int n)
{
    int logn = __builtin_ctz(n);
    int n1 = 1 << (logn / 2);
    int n2 = n / n1;

    vector<Complex> buffer(n);

    transposeParallelBlocked(a.data(), buffer.data(), n1, n2);

#pragma omp parallel for
    for (int i = 0; i < n1; ++i)
    {
        fftIterative(&buffer[i * n2], n2);
    }

#pragma omp parallel for collapse(2)
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            double angle = -2.0 * PI * i * j / n;
            buffer[i * n2 + j] *= Complex(cos(angle), sin(angle));
        }
    }

    transposeParallelBlocked(buffer.data(), a.data(), n1, n2);

#pragma omp parallel for
    for (int i = 0; i < n2; ++i)
    {
        fftIterative(&a[i * n1], n1);
    }

    transposeParallelBlocked(a.data(), buffer.data(), n2, n1);

    a = move(buffer);
}

int main()
{
    int n, nthreads;
    cin >> n;


    vector<Complex> signal(n);

    // for (int i = 0; i < n; i++) {
    //     double x; cin >> x;
    //     signal[i] = Complex(x, 0.0);
    // }


    cin >> nthreads;
    omp_set_num_threads(nthreads);

    
    for (int i = 0; i < n; i++)
    {
        double x;
        cin >> x;
        signal[i] = Complex(x, 0.0);
    }

    auto start = chrono::high_resolution_clock::now();
    baileyFFTParallel(signal, n);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << fixed << setprecision(6) << elapsed.count() << " seconds\n";

    return 0;
}