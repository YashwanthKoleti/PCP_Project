// Bailey’s Parallel Algorithm

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <thread>
#include <algorithm>
#include <functional>
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

void transposeWorker(const vector<Complex> &a, vector<Complex> &temp, int n1, int n2, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            temp[j * n1 + i] = a[i * n2 + j];
        }
    }
}

void rowWorker1(vector<Complex> &a, int n2, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; ++i)
    {
        vector<Complex> row(n2);
        for (int j = 0; j < n2; ++j)
            row[j] = a[i * n2 + j];
        fftIterative(row);
        for (int j = 0; j < n2; ++j)
            a[i * n2 + j] = row[j];
    }
}

void twiddleWorker(vector<Complex> &a, int n, int n2, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            double angle = -2.0 * PI * i * j / n;
            a[i * n2 + j] *= Complex(cos(angle), sin(angle));
        }
    }
}

void rowWorker2(vector<Complex> &a, int n1, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; ++i)
    {
        vector<Complex> row(n1);
        for (int j = 0; j < n1; ++j)
            row[j] = a[i * n1 + j];
        fftIterative(row);
        for (int j = 0; j < n1; ++j)
            a[i * n1 + j] = row[j];
    }
}

void transposeParallel(vector<Complex> &a, int n1, int n2, int nthreads)
{
    vector<Complex> temp(a.size());
    vector<thread> threads;
    int rows_per_thread = n1 / nthreads;

    for (int t = 0; t < nthreads; ++t)
    {
        int start = t * rows_per_thread;
        int end = (t == nthreads - 1) ? n1 : (t + 1) * rows_per_thread;

        threads.emplace_back(transposeWorker, ref(a), ref(temp), n1, n2, start, end);
    }
    for (auto &th : threads)
        th.join();
    a = move(temp);
}

void baileyFFTParallel(vector<Complex> &a, int n, int nthreads)
{
    int logn = __builtin_ctz(n);
    int n1 = 1 << (logn / 2);
    int n2 = n / n1;

    transposeParallel(a, n1, n2, nthreads);

    vector<thread> threads;
    int rpt1 = n1 / nthreads;
    for (int t = 0; t < nthreads; ++t)
    {
        int start = t * rpt1;
        int end = (t == nthreads - 1) ? n1 : (t + 1) * rpt1;
        threads.emplace_back(rowWorker1, ref(a), n2, start, end);
    }
    for (auto &th : threads)
        th.join();
    threads.clear();

    for (int t = 0; t < nthreads; ++t)
    {
        int start = t * rpt1;
        int end = (t == nthreads - 1) ? n1 : (t + 1) * rpt1;
        threads.emplace_back(twiddleWorker, ref(a), n, n2, start, end);
    }
    for (auto &th : threads)
        th.join();
    threads.clear();

    transposeParallel(a, n1, n2, nthreads);

    int rpt2 = n2 / nthreads;
    for (int t = 0; t < nthreads; ++t)
    {
        int start = t * rpt2;
        int end = (t == nthreads - 1) ? n2 : (t + 1) * rpt2;
        threads.emplace_back(rowWorker2, ref(a), n1, start, end);
    }
    for (auto &th : threads)
        th.join();
    threads.clear();

    transposeParallel(a, n2, n1, nthreads);
}

int main()
{
    int n, nthreads;
    cin>>n;

    vector<Complex> signal(n);

    // for (int i = 0; i < n; i++) {
    //     double x; cin >> x;
    //     signal[i] = Complex(x, 0.0);
    // }

    cin >> nthreads;

    auto start = chrono::high_resolution_clock::now();
    baileyFFTParallel(signal, n, nthreads);
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