// cooley tukey Iterative Parallel with pthreads

#include <pthread.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include <cmath>

using namespace std;

using Complex = complex<double>;
const double PI = acos(-1.0);

#ifdef __APPLE__
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int total;
    int generation;
} pthread_barrier_t;

int pthread_barrier_init(pthread_barrier_t *b, void *, int count)
{
    b->total = count;
    b->count = 0;
    b->generation = 0;
    pthread_mutex_init(&b->mutex, NULL);
    pthread_cond_init(&b->cond, NULL);
    return 0;
}

int pthread_barrier_wait(pthread_barrier_t *b)
{
    pthread_mutex_lock(&b->mutex);
    int gen = b->generation;
    b->count++;
    if (b->count == b->total)
    {
        b->count = 0;
        b->generation++;
        pthread_cond_broadcast(&b->cond);
        pthread_mutex_unlock(&b->mutex);
        return 1;
    }
    while (gen == b->generation)
        pthread_cond_wait(&b->cond, &b->mutex);
    pthread_mutex_unlock(&b->mutex);
    return 0;
}

int pthread_barrier_destroy(pthread_barrier_t *b)
{
    pthread_mutex_destroy(&b->mutex);
    pthread_cond_destroy(&b->cond);
    return 0;
}
#endif

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

struct ThreadArgs
{
    int tid;
    int num_threads;
    Complex *a;
    int n;
    int num_stages;
    pthread_barrier_t *barrier;
};

void *fft_worker(void *arg)
{
    ThreadArgs *ctx = (ThreadArgs *)arg;
    const int tid = ctx->tid;
    const int T = ctx->num_threads;
    Complex *a = ctx->a;
    const int n = ctx->n;
    pthread_barrier_t *barrier = ctx->barrier;

    for (int s = 0; s < ctx->num_stages; s++)
    {
        const int half = stage_info[s].count;
        const int len = half << 1;
        const int roff = stage_info[s].offset;
        const int num_groups = n / len;
        const int total_butterflies = num_groups * half;

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

        for (int b = b_start; b < b_end; b++)
        {
            int g = b / half;
            int j = b % half;
            int base = g * len;
            int idx = base + j;

            Complex u = a[idx];
            Complex v = a[idx + half] * all_roots[roff + j];
            a[idx] = u + v;
            a[idx + half] = u - v;
        }

        pthread_barrier_wait(barrier);
    }
    return nullptr;
}

void fftIterative(vector<Complex> &a, int NUM_THREADS)
{
    int n = a.size();
    bitReverse(a);
    precomputeAllRoots(n);

    int num_stages = stage_info.size();

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, NUM_THREADS);

    vector<pthread_t> threads(NUM_THREADS);
    vector<ThreadArgs> args(NUM_THREADS);

    for (int i = 0; i < NUM_THREADS; i++)
    {
        args[i] = {i, NUM_THREADS, a.data(), n, num_stages, &barrier};
        pthread_create(&threads[i], NULL, fft_worker, &args[i]);
    }

    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }

    pthread_barrier_destroy(&barrier);
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
    fftIterative(signal, nthreads);
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