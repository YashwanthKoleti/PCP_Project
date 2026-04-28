#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using Complex = std::complex<double>;
const double PI = std::acos(-1.0);

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

void bitReverse(std::vector<Complex> &a)
{
    int n = a.size();
    int logn = __builtin_ctz(n);

    for (int i = 0; i < n; i++)
    {
        int j = reverseBits(i, logn);

        if (i < j)
        {
            std::swap(a[i], a[j]);
        }
    }
}

void fftIterative(std::vector<Complex> &a)
{
    int n = a.size();

    bitReverse(a);

    for (int len = 2; len <= n; len <<= 1)
    {

        double angle = -2.0 * PI / len;
        Complex wlen(std::cos(angle), std::sin(angle));

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
    std::cin>>signal_lenght;
    std::vector<Complex> signal(signal_lenght);
    for(int i = 0;i < signal_lenght;i++)
    {
        double x;
        std::cin >> x;
        signal[i] = Complex(x, 0.0);
    }

    std::cout << "\n\n";

    auto start = std::chrono::high_resolution_clock::now();
    fftIterative(signal);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time: " << elapsed.count() << " seconds\n";

    // for (const auto &val : signal)
    // {

    //     double real = (std::abs(val.real()) < 1e-10) ? 0.0 : val.real();
    //     double imag = (std::abs(val.imag()) < 1e-10) ? 0.0 : val.imag();

    //     std::cout << "(" << std::fixed << std::setprecision(4) << real
    //               << ", " << imag << ")\n";
    // }

    return 0;
}