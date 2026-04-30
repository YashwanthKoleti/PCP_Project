# Parallel FFT Project

This project contains sequential and parallel C++ implementations of FFT. Results, plots, and analysis are provided in `report.pdf`.

## Files

Sequential implementations:

- `src/sequential/Src_Prjt-Co23btech11009-file1.cpp` - Recursive Cooley-Tukey FFT
- `src/sequential/Src_Prjt-Co23btech11009-file2.cpp` - Iterative Cooley-Tukey FFT
- `src/sequential/Src_Prjt-Co23btech11009-file3.cpp` - Sequential Bailey FFT

Parallel implementations:

- `src/parallel/Src_Prjt-Co23btech11009-file4.cpp` - Parallel recursive Cooley-Tukey FFT
- `src/parallel/Src_Prjt-Co23btech11009-file5.cpp` - Parallel iterative Cooley-Tukey FFT
- `src/parallel/Src_Prjt-Co23btech11009-file6.cpp` - Parallel Bailey FFT

## Requirements

- `g++` with C++17 support
- POSIX threads support

Input size `N` should be a power of 2.

## Compilation

Run these commands from the project root.

```bash
g++ src/sequential/Src_Prjt-Co23btech11009-file1.cpp -o file1 -O2 -std=c++17
g++ src/sequential/Src_Prjt-Co23btech11009-file2.cpp -o file2 -O2 -std=c++17
g++ src/sequential/Src_Prjt-Co23btech11009-file3.cpp -o file3 -O2 -std=c++17

g++ src/parallel/Src_Prjt-Co23btech11009-file4.cpp -o file4 -O2 -std=c++17 -pthread
g++ src/parallel/Src_Prjt-Co23btech11009-file5.cpp -o file5 -O2 -std=c++17 -pthread
g++ src/parallel/Src_Prjt-Co23btech11009-file6.cpp -o file6 -O2 -std=c++17 -pthread
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

Parallel recursive and iterative programs `file4` and `file5` read `N` and number of threads.

```bash
./file4
```

Sample input:

```txt
8
4
```

Parallel Bailey program `file6` reads `N`, number of threads, and `N` real signal values.

```bash
./file6
```

Sample input:

```txt
8
4
1.0
2.0
3.0
4.0
5.0
6.0
7.0
8.0
```

## Output

All programs print output in this format:

```txt
Time: <seconds> seconds
```

The FFT values are not printed because the programs are used for timing experiments.