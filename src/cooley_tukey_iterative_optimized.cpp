#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

struct Loop
{
    int start;
    int end;
    std::vector<Complex> *a;
    int len;
    int step;
};

void bitReverse(std::vector<Complex> &a)
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
        {
            std::swap(a[i], a[j]);
        }
    }
}

std::vector<Complex> roots;
void precompute(int n) {
    roots.resize(n);
    for (int i = 0; i < n; i++) {
        double angle = -2.0 * PI * i / n;
        roots[i] = Complex(std::cos(angle), std::sin(angle));
    }
}

void *opti(void *args)
{
    Loop *para = (Loop *)args;
    std::vector<Complex> &a = *(para->a);
    int len = para->len;
    int step = para->step;
    int start = para->start;
    int end = para->end;

    for (int i = start; i < end; i += len)
    {
        for (int j = 0; j < len / 2; ++j)
        {
            Complex u = a[i + j];
            Complex v = a[i + j + len / 2] * roots[j * step];

            a[i + j] = u + v;
            a[i + j + len / 2] = u - v;
        }
    }

    return nullptr;
}

void fftIterative(std::vector<Complex> &a, int NUM_THREAD)
{
    int n = a.size();

    bitReverse(a);
    precompute(n);

    std::vector<pthread_t> loop_thread(NUM_THREAD);
    std::vector<Loop> thread_data(NUM_THREAD);

    for (int len = 2; len <= n; len <<= 1)
    {
        int step = n / len;
        int chunk = (step + NUM_THREAD - 1) / NUM_THREAD;

        for (int i = 0; i < NUM_THREAD; i++)
        {
            thread_data[i].start = i * chunk;
            thread_data[i].end = std::min((i + 1) * chunk, n);
            thread_data[i].a = &a;
            thread_data[i].len = len;
            thread_data[i].step = step;
            pthread_create(&loop_thread[i], NULL, opti, &thread_data[i]);
        }

        for (int i = 0; i < NUM_THREAD; i++)
        {
            pthread_join(loop_thread[i], NULL);
        }
    }
}

int main()
{
    int NUM_THREADS;
    std::cin >> NUM_THREADS;

    int signal_lenght;
    std::cin>>signal_lenght;
    std::vector<Complex> signal(signal_lenght);
    for(int i = 0;i < signal_lenght;i++)
    {
        double x;
        std::cin >> x;
        signal[i] = Complex(x, 0.0);
    }

    int n = signal.size();
    

    roots.resize(n);
    fftIterative(signal, NUM_THREADS);

    std::cout << "FFT Result (Frequency Domain):\n";
    for (const auto &val : signal)
    {

        double real = (std::abs(val.real()) < 1e-10) ? 0.0 : val.real();
        double imag = (std::abs(val.imag()) < 1e-10) ? 0.0 : val.imag();

        std::cout << "(" << std::fixed << std::setprecision(4) << real
                  << ", " << imag << ")\n";
    }

    return 0;
}