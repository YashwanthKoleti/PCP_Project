#include <omp.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

using Complex = complex<double>;
const double PI = acos(-1.0);

struct StageRoots
{
    int offset;
    int count;
};

vector<Complex> all_roots;
vector<StageRoots> stage_info;

void precomputeAllRoots(int n)
{
    all_roots.resize(n);
    stage_info.clear();

    int offset = 0;
    for (int len = 2; len <= n; len <<= 1)
    {
        int half = len / 2;
        StageRoots sr;
        sr.offset = offset;
        sr.count = half;
        for (int j = 0; j < half; j++)
        {
            double angle = -2.0 * PI * j / len;
            all_roots[offset + j] = Complex(cos(angle), sin(angle));
        }
        stage_info.push_back(sr);
        offset += half;
    }
}

void bitReverse(vector<Complex> &a)
{
    int n = a.size();
    int j = 0;
    for (int i = 1; i < n; i++)
    {
        int bit = n >> 1;
        while (j & bit)
        {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j)
            swap(a[i], a[j]);
    }
}

void fftIterative(vector<Complex> &a, int NUM_THREADS)
{
    int n = a.size();

    if ((n & (n - 1)) != 0)
    {
        cerr << "Error: Signal size N must be a power of 2." << endl;
        return;
    }

    bitReverse(a);
    precomputeAllRoots(n);

    int num_stages = stage_info.size();
    Complex *a_ptr = a.data();

#pragma omp parallel num_threads(NUM_THREADS)
    {

        int tid = omp_get_thread_num();
        int T = omp_get_num_threads();

        const int total_butterflies = n / 2;
        int bpt = total_butterflies / T;
        int rem = total_butterflies % T;
        int b_start, b_end;

        if (tid < rem)
        {
            b_start = tid * (bpt + 1);
            b_end = b_start + bpt + 1;
        }
        else
        {
            b_start = tid * bpt + rem;
            b_end = b_start + bpt;
        }

        for (int s = 0; s < num_stages; s++)
        {
            const int half = stage_info[s].count;
            const int len = half << 1;
            const int roff = stage_info[s].offset;

            int g = b_start / half;
            int j = b_start % half;

            for (int b = b_start; b < b_end; b++)
            {
                int base = g * len;
                int idx = base + j;

                Complex u = a_ptr[idx];
                Complex v = a_ptr[idx + half] * all_roots[roff + j];
                a_ptr[idx] = u + v;
                a_ptr[idx + half] = u - v;

                j++;
                if (j == half)
                {
                    j = 0;
                    g++;
                }
            }

#pragma omp barrier
        }
    }
}

int main()
{
    int n;
    cin >> n;

    int nthreads;

    vector<Complex> signal(n);

    // for (int i = 0; i < n; i++) {
    //     double x; cin >> x;
    //     signal[i] = Complex(x, 0.0);
    // }


    cin >> nthreads;

    mt19937 gen(42);
    uniform_real_distribution<double> dist(-100.0, 100.0);

    for (int i = 0; i < n; i++)
    {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = chrono::high_resolution_clock::now();
    fftIterative(signal, nthreads);
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time: " << fixed << setprecision(6) << elapsed.count() << " seconds\n";

    return 0;
}