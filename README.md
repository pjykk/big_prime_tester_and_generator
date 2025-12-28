# Big Prime Tester and Generator

A very simple prime tester and generator using Miller-Rabin test and Montgomery multiplication.

This is a final project of C-programming class.

## How to Use It

It's recommended to compile the file with

```text
gcc BigPrimeGenerator.c -o BigPrimeGenerator -O3
```

for higher efficiency (about 3x faster).

All inputs and outputs are in hexadecimal form.

Detailed Instructions are printed in stdout.

## Average Running Time

For each Miller-Rabin test, it takes about 0.2s (without -O3) or 0.06-0.07s (-O3 enabled).

It takes on average about 80-120 tries to succuessfully generate a prime.

This program measures its running time for each prime generation, detailed things are printed in stdout, like this:

```text
Total Time: 50.368s
Avg: 5.037s per prime; 10 prime(s) generated
Avg: 0.062s per try; 806 try/tries in total
```

Here a 'try' basically means a Miller-Rabin test.
