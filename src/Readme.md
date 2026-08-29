# Parallel FFT Project

This project contains sequential and parallel C++ implementations of FFT. Results, plots, and analysis are provided in `Rpt_Prjt-Co23btech11009.pdf`.

## Files

Sequential implementations:

- `Src_Prjt-Co23btech11009-file1.cpp` - Recursive Cooley-Tukey FFT
- `Src_Prjt-Co23btech11009-file2.cpp` - Iterative Cooley-Tukey FFT
- `Src_Prjt-Co23btech11009-file3.cpp` - Sequential Bailey FFT

Parallel implementations:

- `Src_Prjt-Co23btech11009-file4.cpp` - Parallel recursive Cooley-Tukey FFT
- `Src_Prjt-Co23btech11009-file5.cpp` - Parallel iterative Cooley-Tukey FFT
- `Src_Prjt-Co23btech11009-file6.cpp` - Parallel Bailey FFT

## Requirements

- `g++` with C++17 support
- POSIX threads support

Input size `N` should be a power of 2.

## Compilation

Run these commands from the directory containing the submitted files.

```bash
g++ Src_Prjt-Co23btech11009-file1.cpp -o file1 -O2 -std=c++17
g++ Src_Prjt-Co23btech11009-file2.cpp -o file2 -O2 -std=c++17
g++ Src_Prjt-Co23btech11009-file3.cpp -o file3 -O2 -std=c++17

g++ Src_Prjt-Co23btech11009-file4.cpp -o file4 -O2 -std=c++17 -pthread
g++ Src_Prjt-Co23btech11009-file5.cpp -o file5 -O2 -std=c++17 -pthread
g++ Src_Prjt-Co23btech11009-file6.cpp -o file6 -O2 -std=c++17 -pthread
```

## Execution

Each program reads from standard input and prints execution time.

Sequential programs `file1`, `file2`, and `file3` read only `N`.

```bash
./file1
```

Sample input:

```txt
8
```

All parallel programs `file4`, `file5` and `file6` read `N` and number of threads.

```bash
./file4
```

Sample input:

```txt
8
4
```

## Output

All programs print output in this format:

```txt
Time: <seconds> seconds
```

The FFT values are not printed because the programs are used for timing experiments.

## Signal Input Note

Most programs currently do not take the full signal values as input. They generate the signal internally using a fixed random seed or leave it zero initialized. If custom signal input is required, uncomment the signal input loop inside the `main()` function of the corresponding file:

```cpp
// for (int i = 0; i < n; i++) {
//     double x;
//     cin >> x;
//     signal[i] = Complex(x, 0.0);
// }
```

After uncommenting this block, provide `N` real values at the same position as the uncommented loop in `main()`. In `file4`,`file5` and `file6`, this means the input order becomes `N`, then `N` signal values, then the number of threads.
Also uncomment the below code in `main()` for printing the output signal
```cpp
    // for (const auto& val : signal) {

    //     double real = (abs(val.real()) < 1e-10) ? 0.0 : val.real();
    //     double imag = (abs(val.imag()) < 1e-10) ? 0.0 : val.imag();

    //     cout << "(" << fixed << setprecision(4) << real
    //               << ", " << imag << ")\n";
    // }
```

## Submission Contents

The submitted zip file contains the following files:

- `Src_Prjt-Co23btech11009-file1.cpp`
- `Src_Prjt-Co23btech11009-file2.cpp`
- `Src_Prjt-Co23btech11009-file3.cpp`
- `Src_Prjt-Co23btech11009-file4.cpp`
- `Src_Prjt-Co23btech11009-file5.cpp`
- `Src_Prjt-Co23btech11009-file6.cpp`
- `Readme.txt`
- `Rpt_Prjt-Co23btech11009.pdf`
- `Turitin_Prjt-Co23btech11009.pdf`
