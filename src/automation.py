import subprocess
import csv
import random
import re
import statistics
import matplotlib.pyplot as plt

# ----------------------------
# CONFIG
# ----------------------------
CPP_FILE = "cooley_tukey_recursive_parallel.c++"
EXECUTABLE = "./parallel_fft"
CSV_FILE = "fft_results.csv"

thread_counts = [2**i for i in range(0, 7)]      # 1 to 64
signal_sizes = [2**i for i in range(1, 21)]      # 2 to 2^15
REPEATS = 20   # 🔥 number of runs per configuration

# ----------------------------
# STEP 1: Compile C++ program
# ----------------------------
print("Compiling C++ FFT program...")
compile_cmd = ["g++", "-O2", "-std=c++17", CPP_FILE, "-o", "parallel_fft", "-pthread"]
subprocess.run(compile_cmd, check=True)
print("Compilation successful.\n")

# ----------------------------
# STEP 2: Run experiments
# ----------------------------
results = []

for size in signal_sizes:
    # Generate ONE signal per size (important for fairness)
    signal = [str(random.randint(1, 100)) for _ in range(size)]

    for threads in thread_counts:
        print(f"Running: size={size}, threads={threads}")

        times = []

        for r in range(REPEATS):
            input_data = f"{size}\n" + "\n".join(signal) + f"\n{threads}\n"

            process = subprocess.run(
                [EXECUTABLE],
                input=input_data,
                text=True,
                capture_output=True
            )

            output = process.stdout

            match = re.search(r"Time:\s*([0-9.]+)", output)
            if match:
                times.append(float(match.group(1)))

        # Compute stats
        if times:
            avg_time = statistics.mean(times)
            std_dev = statistics.stdev(times) if len(times) > 1 else 0.0
        else:
            avg_time = None
            std_dev = None

        results.append([size, threads, avg_time, std_dev])

# ----------------------------
# STEP 3: Save results to CSV
# ----------------------------
with open(CSV_FILE, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Signal_Size", "Threads", "Avg_Time", "Std_Dev"])
    writer.writerows(results)

print(f"\nResults saved to {CSV_FILE}")

# ----------------------------
# STEP 4: Plot results
# ----------------------------
plt.figure(figsize=(10, 6))

for threads in thread_counts:
    sizes = [r[0] for r in results if r[1] == threads]
    times = [r[2] for r in results if r[1] == threads]

    plt.plot(sizes, times, marker='o', label=f"{threads} threads")

plt.xscale("log", base=2)
plt.yscale("log")
plt.xlabel("Signal Size (log2 scale)")
plt.ylabel("Execution Time (seconds, log scale)")
plt.title("Parallel FFT Performance (Averaged)")
plt.legend()
plt.grid(True)
plt.show()