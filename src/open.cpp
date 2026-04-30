#include <iostream>
#include <omp.h>

int main() {
    int thread_id, num_threads;

    // Optional: Set the number of threads you want to use. 
    // If you remove this, it will default to your CPU's max thread count.
    omp_set_num_threads(4);

    // This pragma tells the compiler to parallelize the following block
    #pragma omp parallel private(thread_id)
    {
        // Get the specific ID for the current thread
        thread_id = omp_get_thread_num();
        
        // This will print out of order because threads run simultaneously
        cout << "Hello from thread " << thread_id << "\n";

        // The 'master' block is only executed by the master thread (Thread 0)
        #pragma omp master
        {
            num_threads = omp_get_num_threads();
            cout << "--- Master Thread Reporting ---\n";
            cout << "Total threads in action: " << num_threads << "\n";
        }
    }

    return 0;
}