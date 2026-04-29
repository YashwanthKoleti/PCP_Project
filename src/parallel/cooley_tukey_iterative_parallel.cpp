// #include <cmath>
// #include <complex>
// #include <iomanip>
// #include <iostream>
// #include <vector>

// using Complex = std::complex<double>;
// const double PI = std::acos(-1.0);

// struct Loop {
//   int start;
//   int end;
//   std::vector<Complex> *a;
//   int len;
//   int step;
// };

// void bitReverse(std::vector<Complex> &a) {
//   int n = a.size();
//   int j = 0;

//   for (int i = 1; i < n; i++) {
//     int bit = n >> 1;
//     while (j & bit) {
//       j ^= bit;
//       bit >>= 1;
//     }
//     j ^= bit;

//     if (i < j) {
//       std::swap(a[i], a[j]);
//     }
//   }
// }

// std::vector<Complex> roots;
// void precompute(int n) {
//   roots.resize(n);
//   for (int i = 0; i < n; i++) {
//     double angle = -2.0 * PI * i / n;
//     roots[i] = Complex(std::cos(angle), std::sin(angle));
//   }
// }

// void *opti(void *args) {
//   Loop *para = (Loop *)args;
//   std::vector<Complex> &a = *(para->a);
//   int len = para->len;
//   int step = para->step;
//   int start = para->start;
//   int end = para->end;

//   for (int i = start; i < end; i += len) {
//     for (int j = 0; j < len / 2; ++j) {
//       Complex u = a[i + j];
//       Complex v = a[i + j + len / 2] * roots[j * step];

//       a[i + j] = u + v;
//       a[i + j + len / 2] = u - v;
//     }
//   }

//   return nullptr;
// }

// void fftIterative(std::vector<Complex> &a, int NUM_THREAD) {
//   int n = a.size();

//   bitReverse(a);
//   precompute(n);

//   std::vector<pthread_t> loop_thread(NUM_THREAD);
//   std::vector<Loop> thread_data(NUM_THREAD);

//   for (int len = 2; len <= n; len <<= 1) {
//     int step = n / len;
//     int chunk = (step + NUM_THREAD - 1) / NUM_THREAD;

//     for (int i = 0; i < NUM_THREAD; i++) {
//       thread_data[i].start = i * chunk;
//       thread_data[i].end = std::min((i + 1) * chunk, n);
//       thread_data[i].a = &a;
//       thread_data[i].len = len;
//       thread_data[i].step = step;
//       pthread_create(&loop_thread[i], NULL, opti, &thread_data[i]);
//     }

//     for (int i = 0; i < NUM_THREAD; i++) {
//       pthread_join(loop_thread[i], NULL);
//     }
//   }
// }

// int main() {
//   int NUM_THREADS;
//   std::cin >> NUM_THREADS;

//   int signal_lenght;
//   std::cin >> signal_lenght;
//   std::vector<Complex> signal(signal_lenght);
//   for (int i = 0; i < signal_lenght; i++) {
//     double x;
//     std::cin >> x;
//     signal[i] = Complex(x, 0.0);
//   }

//   int n = signal.size();

//   roots.resize(n);
//   auto start = std::chrono::high_resolution_clock::now();
//   fftIterative(signal, NUM_THREADS);
//   auto end = std::chrono::high_resolution_clock::now();

//   std::chrono::duration<double> elapsed = end - start;
//   std::cout << "Time: " << elapsed.count() << " seconds\n";

//   for (const auto &val : signal) {

//     double real = (std::abs(val.real()) < 1e-10) ? 0.0 : val.real();
//     double imag = (std::abs(val.imag()) < 1e-10) ? 0.0 : val.imag();

//     std::cout << "(" << std::fixed << std::setprecision(4) << real << ", "
//               << imag << ")\n";
//   }

//   return 0;
// }
// 
#include <pthread.h>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

// --- Apple Compatibility for pthread_barrier ---
#ifdef __APPLE__
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int total;
    int generation;
} pthread_barrier_t;

int pthread_barrier_init(pthread_barrier_t *b, void *, int count) {
    b->total = count;
    b->count = 0;
    b->generation = 0;
    pthread_mutex_init(&b->mutex, NULL);
    pthread_cond_init(&b->cond, NULL);
    return 0;
}

int pthread_barrier_wait(pthread_barrier_t *b) {
    pthread_mutex_lock(&b->mutex);
    int gen = b->generation;
    b->count++;
    if (b->count == b->total) {
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

int pthread_barrier_destroy(pthread_barrier_t *b) {
    pthread_mutex_destroy(&b->mutex);
    pthread_cond_destroy(&b->cond);
    return 0;
}
#endif
// ---------------------------------------------

struct StageRoots {
    int offset;  // offset into flat array
    int count;   // len/2
};

std::vector<Complex> all_roots;
std::vector<StageRoots> stage_info;

void precomputeAllRoots(int n) {
    all_roots.resize(n);
    stage_info.clear();
    
    int offset = 0;
    for (int len = 2; len <= n; len <<= 1) {
        int half = len / 2;
        StageRoots sr;
        sr.offset = offset;
        sr.count = half;
        for (int j = 0; j < half; j++) {
            double angle = -2.0 * PI * j / len;
            all_roots[offset + j] = Complex(std::cos(angle), std::sin(angle));
        }
        stage_info.push_back(sr);
        offset += half;
    }
}

void bitReverse(std::vector<Complex> &a) {
    int n = a.size();
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        while (j & bit) { j ^= bit; bit >>= 1; }
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

struct ThreadArgs {
    int tid;
    int num_threads;
    Complex *a;          // raw pointer to avoid vector overhead
    int n;
    int num_stages;
    pthread_barrier_t *barrier;
};

void *fft_worker(void *arg) {
    ThreadArgs *ctx = (ThreadArgs *)arg;
    const int tid = ctx->tid;
    const int T   = ctx->num_threads;
    Complex *a    = ctx->a;
    const int n   = ctx->n;
    pthread_barrier_t *barrier = ctx->barrier;

    // Calculate this thread's workload once per thread, not per stage.
    // Total butterflies per stage is always n / 2
    const int total_butterflies = n / 2;  
    int bpt = total_butterflies / T;
    int rem = total_butterflies % T;
    int b_start, b_end;
    
    if (tid < rem) {
        b_start = tid * (bpt + 1);
        b_end   = b_start + bpt + 1;
    } else {
        b_start = tid * bpt + rem;
        b_end   = b_start + bpt;
    }

    for (int s = 0; s < ctx->num_stages; s++) {
        const int half   = stage_info[s].count;       
        const int len    = half << 1;
        const int roff   = stage_info[s].offset;      

        // OPTIMIZATION 1: Calculate 2D indices outside the loop
        // This completely eliminates integer division and modulo in the inner loop.
        int g = b_start / half;
        int j = b_start % half;

        for (int b = b_start; b < b_end; b++) {
            int base = g * len;
            int idx = base + j;
            
            Complex u = a[idx];
            Complex v = a[idx + half] * all_roots[roff + j];
            a[idx]        = u + v;
            a[idx + half] = u - v;

            // Increment local coordinates directly
            j++;
            if (j == half) {
                j = 0;
                g++;
            }
        }

        // Wait for all threads to finish the current stage
        pthread_barrier_wait(barrier);
    }
    return nullptr;
}

void fftIterative(std::vector<Complex> &a, int NUM_THREADS) {
    int n = a.size();
    
    // Check if N is a power of 2 (Required for Radix-2 DIT FFT)
    if ((n & (n - 1)) != 0) {
        std::cerr << "Error: Signal size N must be a power of 2." << std::endl;
        return;
    }

    bitReverse(a);
    precomputeAllRoots(n);

    int num_stages = stage_info.size();

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, NUM_THREADS);

    std::vector<pthread_t> threads(NUM_THREADS);
    std::vector<ThreadArgs> args(NUM_THREADS);

    for (int i = 0; i < NUM_THREADS; i++) {
        args[i] = {i, NUM_THREADS, a.data(), n, num_stages, &barrier};
        pthread_create(&threads[i], NULL, fft_worker, &args[i]);
    }

    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }

    pthread_barrier_destroy(&barrier);
}

int main() {
    int n;
    if (!(std::cin >> n)) return 0;
    
    int nthreads;
    if (!(std::cin >> nthreads)) return 0;

    std::vector<Complex> signal(n);
    std::mt19937 gen(42); 
    std::uniform_real_distribution<double> dist(-100.0, 100.0);
    
    for (int i = 0; i < n; i++) {
        signal[i] = Complex(dist(gen), 0.0);
    }

    auto start = std::chrono::high_resolution_clock::now();
    fftIterative(signal, nthreads);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time: " << std::fixed << std::setprecision(6) << elapsed.count() << " seconds\n";

    return 0;
}