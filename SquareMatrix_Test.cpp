#include <vector>
#include <chrono>
#include <iomanip>
#include <iostream>
#include "SquareMatrix.h"

using namespace std;
int main() {
    vector<size_t> sizes = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512};

    // timing storage
    vector<double> bf_times, dc_times, tdc_times, str_times;

    for (size_t n : sizes) {
        cout << "=== " << n << "x" << n << " Test ===" << endl;

        SquareMatrix A{n}, B{n};

        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) {
                A.data[i][j] = (i * n + j) % 10 + 1;
                B.data[i][j] = (j * n + i) % 10 + 1;
            }

        // time BruteForce
        auto start = chrono::high_resolution_clock::now();
        SquareMatrix* bf = BruteForce(A, B);
        auto end = chrono::high_resolution_clock::now();
        bf_times.push_back(chrono::duration<double, milli>(end - start).count());

        // time DivideAndConquer
        start = chrono::high_resolution_clock::now();
        SquareMatrix* dc = DivideAndConquer(A, B);
        end = chrono::high_resolution_clock::now();
        dc_times.push_back(chrono::duration<double, milli>(end - start).count());

        // time ThreadedDivideAndConquer
        start = chrono::high_resolution_clock::now();
        SquareMatrix* tdc = ThreadedDivideAndConquer(A, B);
        end = chrono::high_resolution_clock::now();
        tdc_times.push_back(chrono::duration<double, milli>(end - start).count());

        // time Strassen
        start = chrono::high_resolution_clock::now();
        SquareMatrix* str = Strassen(A, B);
        end = chrono::high_resolution_clock::now();
        str_times.push_back(chrono::duration<double, milli>(end - start).count());

        // correctness check
        bool dc_pass = true, tdc_pass = true, str_pass = true;
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) {
                if (bf->data[i][j] != dc->data[i][j])  dc_pass  = false;
                if (bf->data[i][j] != tdc->data[i][j]) tdc_pass = false;
                if (bf->data[i][j] != str->data[i][j]) str_pass = false;
            }

        cout << "DivideAndConquer:         " << (dc_pass  ? "PASS" : "FAIL") << endl;
        cout << "ThreadedDivideAndConquer: " << (tdc_pass ? "PASS" : "FAIL") << endl;
        cout << "Strassen:                 " << (str_pass ? "PASS" : "FAIL") << endl;

        if (n <= 4) {
            cout << "BruteForce:"               << endl; bf->display(cout);
            cout << "DivideAndConquer:"         << endl; dc->display(cout);
            cout << "ThreadedDivideAndConquer:" << endl; tdc->display(cout);
            cout << "Strassen:"                 << endl; str->display(cout);
        }

        delete bf;
        delete dc;
        delete tdc;
        delete str;

        cout << endl;
    }

    // display timing table
    cout << endl << "=== Performance Results (ms) ===" << endl << endl;
    cout << setw(6)  << "n"
         << setw(15) << "BruteForce"
         << setw(15) << "D&C"
         << setw(15) << "Strassen"
         << setw(15) << "TD&C" << endl;
    cout << string(66, '-') << endl;

    for (size_t i = 0; i < sizes.size(); i++) {
        cout << setw(6)  << sizes[i]
             << setw(15) << fixed << setprecision(4) << bf_times[i]
             << setw(15) << fixed << setprecision(4) << dc_times[i]
             << setw(15) << fixed << setprecision(4) << str_times[i]
             << setw(15) << fixed << setprecision(4) << tdc_times[i]
             << endl;
    }

    return 0;
}