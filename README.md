# Program 2 - SquareMatrix Multiplication
> Derived from Dr. Bill Booth's Programming Project 2
>
> **Author:** Oluwatosin Ojo

## Files
- **SquareMatrix.h** — defines the `SquareMatrix` struct with constructor, destructor, and display method, along with all helper function declarations
- **SquareMatrix.cpp** — implements four matrix multiplication algorithms
- **SquareMatrix_Test.cpp** — test harness that verifies correctness of all four methods against BruteForce and measures runtime in milliseconds across matrix sizes 1, 2, 4, 8, 16, 32, 64, 128, 256, and 512

## Algorithms

### BruteForce
Standard triple-loop matrix multiplication. O(n³) with minimal overhead — simple and fast in practice due to tight loop structure with no function call overhead.

### DivideAndConquer
Recursively splits matrices into quadrants down to 1x1 elements. Also O(n³) by the Master Theorem (8 subproblems of size n/2), but carries significant recursion overhead — for a 512x512 matrix this means roughly 134 million function calls. Notably slower than BruteForce in practice despite identical theoretical complexity.

### ThreadedDivideAndConquer
Same single-level divide and conquer approach but computes the 4 output quadrants in parallel using 4 POSIX threads. Demonstrates a clear tradeoff — thread creation and joining overhead makes it slower than BruteForce for small matrices (1x1 to ~16x16), but at large sizes (256x256+) the parallelism pays off significantly.

### Strassen
Reduces multiplications from 8 to 7 per recursive level, achieving O(n^2.807). Uses a single level of Strassen decomposition before handing off to BruteForce for subproblems. Faster than BruteForce at large sizes, but the added complexity of 10 S-matrices and 7 P-matrices makes it slower for very small inputs.

## How to Run
Compile and run `SquareMatrix_Test.cpp` with `SquareMatrix.cpp`:
```bash
g++ -o test SquareMatrix.cpp SquareMatrix_Test.cpp -lpthread
./test
```

## Testing
For each matrix size the test harness:
1. Fills matrices A and B with values `(i*n + j) % 10 + 1` to keep values small and avoid overflow
2. Runs and times all four methods
3. Compares each result element-by-element against BruteForce and reports PASS or FAIL
4. Displays actual matrix values for sizes up to 4x4
5. Prints a final performance table showing runtime in milliseconds for all methods at all sizes

## Key Takeaways
| Algorithm | Strength | Weakness |
|---|---|---|
| BruteForce | Fast in practice, simple | Scales poorly at very large sizes |
| DivideAndConquer | Correct, demonstrates recursion | Slowest in practice due to recursion depth overhead |
| ThreadedDivideAndConquer | Best for large matrices | Thread overhead hurts at small sizes |
| Strassen | Best asymptotic complexity | Overhead not worth it for small matrices |
