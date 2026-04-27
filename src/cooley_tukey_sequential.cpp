#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

using Complex = std::complex<double>;

const double PI = std::acos(-1.0);

void fft(std::vector<Complex> &a)
{
    int n = a.size();

    if (n <= 1)
    {
        return;
    }

    std::vector<Complex> even(n / 2);
    std::vector<Complex> odd(n / 2);
    for (int i = 0; i < n / 2; ++i)
    {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    fft(even);
    fft(odd);

    for (int k = 0; k < n / 2; ++k)
    {

        Complex t = std::polar(1.0, -2.0 * PI * k / n) * odd[k];

        a[k] = even[k] + t;
        a[k + n / 2] = even[k] - t;
    }
}

int main() {
    int n;
    std::cin >> n;

    std::vector<Complex> signal(n);
    for (int i = 0; i < n; i++) {
        double x;
        std::cin >> x;
        signal[i] = Complex(x, 0.0);
    }

    auto start = std::chrono::high_resolution_clock::now();
    fft(signal);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time: " << elapsed.count() << " seconds\n";

    return 0;
}