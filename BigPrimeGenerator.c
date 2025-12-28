/*
Big Prime Tester & Generator by pjykk

Randomly test primes in [2, 2^1024) & generate primes in (2^(BIT-1), 2^BIT) (BIT-bit), with BIT <= MAX_BIT = 1024
using Miller-Rabin test and Montgomery mul-mod.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_BIT 1024

int BIT = 1024, PRIME_CNT = 1;
// PRIME_CNT: number of primes generated

#define TEST_CNT 10
// TEST_CNT: count of test numbers for each prime

#define PRIME_TABLE_LENGTH 70
// number of primes in table used, may be within [1, 100]

#define BLOCK 32
// 2 ^ 32 a block in mul (unsigned int)

typedef unsigned int uint;
typedef unsigned long long ull;

// just list them to use the feature of static const (faster)
static const int prime_table[100] = 
{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 
31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 
467, 479, 487, 491, 499, 503, 509, 521, 523, 541};

typedef struct 
{
    int len; // from num[0] to num[len-1] are valid; x=0 <=> len=0
    int num[MAX_BIT * 2 + 10]; // 0/1 per digit; x=num[0]+2*num[1]+4*num[2]+...
}BigInt;

BigInt Trans(int x) // trans int into BigInt
{
    BigInt res;
    res.len = res.num[0] = 0;
    while(x)
    {
        res.num[res.len] = (x & 1);
        res.len++;
        x >>= 1;
    }
    return res;
}

void PrintHex(FILE* stream, BigInt x) // print x in hex form to stream (file or stdout)
{
    x.num[x.len] = x.num[x.len + 1] = x.num[x.len + 2] = 0;
    for(int i = (x.len + 3) / 4 - 1; i >= 0; i--)
    {
        int now = x.num[i * 4] + 2 * x.num[i * 4 + 1] + 4 * x.num[i * 4 + 2] + 8 * x.num[i * 4 + 3];
        if(now <= 9)
        {
            fprintf(stream, "%d", now);
        }
        else
        {
            fprintf(stream, "%c", now - 10 + 'A');
        }
    }
}

int Compare(BigInt a, BigInt b) // return sgn(a-b): -1/0/1
{
    if(a.len < b.len) // first len
    {
        return -1;
    }
    if(a.len > b.len)
    {
        return 1;
    }
    for(int i = a.len - 1; i >= 0; i--) // then number from high to low
    {
        if(a.num[i] < b.num[i])
        {
            return -1;
        }
        if(a.num[i] > b.num[i])
        {
            return 1;
        }
    }
    return 0;
}

BigInt Shift(BigInt a, int x) // return a >> x (x>0) or a << (-x) (x<0)
{
    if(x > 0) // a >> x
    {
        if(a.len <= x)
        {
            a.len = 0;
            return a;
        }
        a.len -= x;
        for(int i = 0; i < a.len; i++)
        {
            a.num[i] = a.num[i + x];
        }
        return a;
    }
    if(x < 0) // a << (-x)
    {
        x = -x;
        for(int i = a.len - 1; i >= 0; i--)
        {
            a.num[i + x] = a.num[i];
        }
        for(int i = 0; i < x; i++)
        {
            a.num[i] = 0;
        }
        a.len += x;
        return a;
    }
    return a; // nothing happen when x = 0
}

BigInt Cut(BigInt a, int x) // return the last x digits of a
{
    if(a.len <= x)
    {
        return a;
    }
    a.len = x;
    while(a.len > 0 && a.num[a.len - 1] == 0) // note that this operation can introduce prefix 0
    {
        a.len--;
    }
    return a;
}

void Swap(BigInt* pa, BigInt* pb) // swap *pa and *pb
{
    BigInt c;
    c = *pa;
    *pa = *pb;
    *pb = c;
}

BigInt Add(BigInt a, BigInt b) // return a + b
{
    if(a.len < b.len)
    {
        Swap(&a, &b);
    }
    BigInt c; 
    c.len = a.len;
    int add = 0;
    for(int i = 0; i < a.len; i++)
    {
        c.num[i] = a.num[i] + add;
        if(i < b.len)
        {
            c.num[i] += b.num[i]; // note that b[pos] (pos >= b.len) can be non-zero
        }
        if(c.num[i] > 1)
        {
            add = 1;
            c.num[i] -= 2;
        }
        else
        {
            add = 0;
        }
    }
    if(add) // extra 1 added
    {
        c.num[c.len] = add;
        c.len++;
    }
    return c;
}

BigInt Sub(BigInt a, BigInt b) // return a - b, only when a >= b (comparison not included)
{
    BigInt c; 
    c.len = a.len;
    int sub = 0;
    for(int i = 0; i < a.len; i++)
    {
        c.num[i] = a.num[i] - sub;
        if(i < b.len)
        {
            c.num[i] -= b.num[i];
        }
        if(c.num[i] < 0)
        {
            c.num[i] += 2;
            sub = 1;
        }
        else
        {
            sub = 0;
        }
    }
    while(c.len > 0 && c.num[c.len - 1] == 0) // note that this operation can introduce prefix 0
    {
        c.len--;
    }
    return c;
}

BigInt Mul(BigInt a, BigInt b)
{
    if(a.len < b.len)
    {
        Swap(&a, &b);
    }
    BigInt c;
    if(b.len == 0)
    {
        c.len = 0;
        return c;
    }
    
    // turn bit-operation into unsigned int-operation
    // calloc & free is much faster than global array & clear (by test)
    int ua_len = (a.len + BLOCK - 1) / BLOCK, ub_len = (b.len + BLOCK - 1) / BLOCK;
    int uc_len = ua_len + ub_len;
    uint* ua = (uint*)calloc(ua_len, sizeof(uint));
    uint* ub = (uint*)calloc(ub_len, sizeof(uint));
    uint* uc = (uint*)calloc(uc_len, sizeof(uint));
    
    for(int i = 0; i < a.len; i++) // turn bits into unsigned int
    {
        if(a.num[i])
        {
            ua[i / BLOCK] |= (1u << (i % BLOCK));
        }
    }
    for(int i = 0; i < b.len; i++)
    {
        if(b.num[i])
        {
            ub[i / BLOCK] |= (1u << (i % BLOCK));
        }
    }
    
    for(int i = 0; i < ua_len; i++) // mul in unsigned int
    {
        if(ua[i] == 0)
        {
            continue;
        }
        for(int j = 0; j < ub_len; j++)
        {
            if(ub[j] == 0)
            {
                continue;
            }
            ull prod = (ull)ua[i] * ub[j];
            uint low = (uint)prod, high = (uint)(prod >> BLOCK); // split the prod into high + low * 2^BLOCK
            ull now = low;
            int k = i + j;
            while(now)
            {
                now += uc[k];
                uc[k] = (uint)now;
                now >>= BLOCK;
                k++;
            }
            now = high;
            k = i + j + 1;
            while(now)
            {
                now += uc[k];
                uc[k] = (uint)now;
                now >>= BLOCK;
                k++;
            }
        }
    }

    c.len = uc_len * BLOCK; // calculate c.len
    while(c.len > 0 && ((uc[(c.len - 1) / BLOCK] >> ((c.len - 1) % BLOCK)) & 1) == 0)
    {
        c.len--;
    }
    for(int i = 0; i < c.len; i++) // turn unsigned int into bits
    {
        if((uc[i / BLOCK] >> (i % BLOCK)) & 1)
        {
            c.num[i] = 1;
        }
        else
        {
            c.num[i] = 0;
        }
    }

    free(ua);
    free(ub);
    free(uc);    

    return c;
}

BigInt N, R, r, M;
/*
some nums used in Montgomery
R = 2 ^ BIT
N: primelike in (R/2, R) to be tested
r = R^2 mod N
M satisfies NM = -1 (mod R)
*/

BigInt ONE; // 1, useful in Miller-Rabin

void Init() // init srand, R and ONE
{
    srand(time(0));
    R.len = BIT + 1;
    R.num[BIT] = 1;
    for(int i = 0; i < R.len - 1; i++) // must be cleared (!)
    {
        R.num[i] = 0;
    }
    ONE.len = 1;
    ONE.num[0] = 1;
}

void Mont_Init() // for given N, init r and M
{
    // r should be 2 ^ {2 * BIT} mod N. just times 2 and minus N each time
    r = Sub(R, N);
    for(int i = 0; i < BIT; i++)
    {
        r = Shift(r, -1);
        if(Compare(r, N) >= 0)
        {
            r = Sub(r, N);
        }
    }

    // from mod 2 ^ {k - 1} to 2 ^ k: inv remains or becomes inv + 2 ^ {k - 1}
    M = ONE; 
    for(int i = 2; i <= BIT; i++)
    {
        if(Compare(Cut(Mul(N, M), i), ONE)!=0)
        {
            M = Add(M, Shift(ONE, -(i-1)));
        }
    }
    M = Sub(R, M); // M = -inv
}

// use suffix -M to represent num in Montgomery form (not a*M)
// then aM = Mont(a*r), a = Mont(aM)
BigInt Mont(BigInt T) // calculate T/R mod N; T should be within [0, R*N)
{
    BigInt m = Cut(Mul(Cut(T, BIT), M), BIT);
    BigInt a = Shift(Add(T, Mul(m, N)), BIT);
    if(Compare(a, N) >= 0) // here, a must be within [0, 2*N)
    {
        a = Sub(a, N);
    }
    return a;
}

BigInt MulMod(BigInt a, BigInt b) // return a*b (mod N); a, b should be within [0, R)
{
    BigInt aM = Mont(Mul(a,r)), bM = Mont(Mul(b,r)); 
    return Mont(Mont(Mul(aM, bM)));
}

BigInt PowMod(BigInt a, BigInt b) // return a^b (mod N); a should be within [0, R)
{
    BigInt ans = Mont(r), now = Mont(Mul(a,r));
    for(int i = 0; i < b.len; i++)
    {
        if(b.num[i] == 1)
        {
            ans = Mont(Mul(ans, now));
        }
        now = Mont(Mul(now, now));
    }
    return Mont(ans);
}

BigInt RandTest() // random integer in [2, N) used for prime test
{
    BigInt x;
    // we randomly generate integer in [0,R) and it is expected to be in [2, N) in about 2 tries
    while(1)
    {
        x.len = BIT;
        for(int i = 0; i < x.len; i++)
        {
            x.num[i] = rand() & 1;
        }
        while(x.len > 0 && x.num[x.len - 1] == 0)
        {
            x.len--;
        }
        if(Compare(N, x) == 1 && Compare(x, ONE) == 1)
        {
            return x;
        }
    }
}

int Div_able(BigInt x, int k) // test if x % k == 0
{
    int sum = 0, base = 1;
    for(int i = 0; i < x.len; i++)
    {
        sum = (sum + base * x.num[i]) % k;
        base = base * 2 % k;
    }
    return (sum == 0);
}

BigInt RandPrimeLike() // random integer in (R/2, R), which is not multiple of small primes
{
    BigInt x;
    x.len = BIT;
    x.num[BIT - 1] = 1; // final x = R/2 + (random gen)
    x.num[0] = 1; // should be odd when BIT > 1
    while(1)
    {
        for(int i = 1; i < x.len - 1; i++)
        {
            x.num[i] = rand() & 1; // random gen
        }
        int flag = 1;
        for(int i = 1; i < PRIME_TABLE_LENGTH; i++) // prime_table[0] = 2, useless
        {
            if(Compare(x, Trans(prime_table[i])) == 0)
            {
                return x;
            }
            if(Div_able(x, prime_table[i]))
            {
                flag = 0;
                break;
            }
        }
        if(flag)
        {
            return x;
        }
    }
    
}

// test if N in (R/2, R) is prime; return [N is prime]
int MillerRabin(int op) // op = 1: print results for each try
{
    Mont_Init(); // necessary for the use of MulMod and PowMod

    int flag = 1; // need to deal with op = 1

    int pow2 = 0;
    BigInt now = Sub(N, ONE); // now = N-1, useful in the test below
    while(pow2 < now.len && now.num[pow2] == 0)
    {
        pow2++;
    }
    BigInt res = Shift(now, pow2); // N = 2^t * res, where res === 1 (mod 2)

    for (int i = 1; i <= TEST_CNT; i++)
    {
        BigInt a = RandTest();
        if(op)
        {
            fprintf(stdout, "Test #%d: ", i);
            PrintHex(stdout, a);
            fprintf(stdout, ": ");
        }
        BigInt num = PowMod(a, res); // init num = a^res
        if(Compare(num, ONE) == 0)
        {
            if(op)
            {
                fprintf(stdout, "Not Sure\n");
            }
            continue;
        }
        int j = 0;
        for( ; j < pow2; j++)
        {
            if(Compare(num, now) == 0) // num = N-1
            {
                break; 
            }
            num = MulMod(num, num); // num = a ^ {2^j * res}
        }
        if(j == pow2) // no 1 or N-1 => N is composite
        {
            flag = 0;
            if(op)
            {
                fprintf(stdout, "Composite\n");
                continue; // for op = 1 we need Miller-Rabin exactly 10 times
            }
            else
            {
                return 0;
            } 
        }
        if(op)
        {
            fprintf(stdout, "Not Sure\n");
        }
    }
    if(op)
    {
        if(flag == 0)
        {
            fprintf(stdout, "\nFinal result: Composite\n\n");
            return 0;
        }
        else
        {
            fprintf(stdout, "\nFinal Result: Very very very likely to be a prime\n\n");
            return 1;
        }
    }
    else
    {
        return 1;
    }
}

int main()
{
    //Instructions
    printf("\n\nBig Prime Tester & Generator.\n\n\n");
    printf("This program has so far two functions: primality testing and prime generating.\n");
    printf("Please carefully read the guideline below first.\n\n");

    printf("Primality testing is suitable for numbers in [0, 2^1024).\n");
    printf("If you want to test whether a hexadecimal number is a prime number, input -test [NUMBER]- (without slashes and brackets) and press enter.\n");
    printf("Both abcdef and ABCDEF in hexadecimal are accepted.\n");
    printf("The program will output the 10 random numbers (in hexadecimal form) used in primality testing and the corresponding results.\n");
    printf("In the last line the program gives the final result. Notice that Miller-Rabin algorithm used here has an error rate of 4^(-10).\n\n");
    
    printf("If you want to generate primes, it's expected to input -gen [BIT] [COUNT]- (without slashes and brackets) and press enter.\n");
    printf("The program will randomly generate COUNT prime(s) in [2^(BIT-1), 2^BIT) (BIT-bit) in hexadecimal form and output them both on the screen and in BigPrime.txt.\n");
    printf("BIT should be an integer within [2, 1024] and COUNT should be an integer for at least 1.\n");
    printf("It may take from less than 1 second to over 15 seconds to generate a single prime. Please wait patiently.\n");

    printf("If you want to end this program (not when producing primes), input -end- (without slashes).");

    printf("Please follow the input formats strictly. Other inputs may cause unexpected problems.\n\n");

    char input[6], test_num[270]; // 1024 / 4 = 256, better to be a bit larger  

    while(1)
    {
        printf("Your input:\n\n");
        scanf("%s", input);
        if(input[0] == 't' && input[1] == 'e' && input[2] == 's' && input[3] == 't' && input[4] == '\0')
        {
            scanf("%s", test_num);
            printf("\n");

            int str_len = 0;
            while(test_num[str_len] != '\0') // strlen by hand
            {
                str_len++;
            }
            if(str_len > 256)
            {
                printf("ERROR: NUMBER too large\n\n");
                continue;
            }

            // init N and BIT from s
            int pos = 0, now, error_flag = 0;
            for(int i = str_len - 1; i >= 0; i--)
            {
                if(test_num[i] >= '0' && test_num[i] <= '9')
                {
                    now = test_num[i] - '0';
                }
                else if(test_num[i] >= 'A' && test_num[i] <= 'F')
                {
                    now = test_num[i] - 'A' + 10;
                }
                else if(test_num[i] >= 'a' && test_num[i] <= 'f')
                {
                    now = test_num[i] - 'a' + 10;
                }
                else
                {
                    printf("ERROR: Invalid hexadecimal\n\n");
                    error_flag = 1;
                    break;
                }
                N.num[pos * 4] = (now & 1);
                N.num[pos * 4 + 1] = (now & 2) >> 1;
                N.num[pos * 4 + 2] = (now & 4) >> 2;
                N.num[pos * 4 + 3] = (now & 8) >> 3;
                N.len = pos * 4 + 4;
                pos++;
            }
            if(error_flag)
            {
                continue;
            }
            while(N.len > 0 && N.num[N.len - 1] == 0) // eliminate prefix 0
            {
                N.len--;
            }
            if(N.len == 0 || N.len == 1) // 0 and 1
            {
                printf("Not Prime\n\n");
                continue;
            }
            if(Compare(N, Trans(2)) == 0) // 2 can't be judged by Miller-Rabin, should be specially treated
            {
                printf("Prime\n\n");
                continue;
            }

            if(N.num[0] == 0) // even numbers can't do Montgomery 
            {
                printf("Composite because it's even\n\n");
                continue;
            }
            
            BIT = N.len;
            Init();
            MillerRabin(1);
        }
        else if(input[0] == 'g' && input[1] == 'e' && input[2] == 'n' && input[3] == '\0')
        {
            scanf("%d%d", &BIT, &PRIME_CNT);
            printf("\n");

            if(BIT <= 1 || BIT > MAX_BIT)
            {
                printf("ERROR: Invalid BIT\n\n");
                continue;
            }
            if(PRIME_CNT <= 0)
            {
                printf("ERROR: Invalid COUNT\n\n");
                continue;
            }

            FILE* fp;
            fp = fopen("BigPrime.txt", "w");
            if(fp == NULL)
            {
                printf("ERROR: Failed to open file BigPrime.txt\n\n");
                continue;
            }

            clock_t start, end;
            double total_time;
            start = clock(); 

            int total_tries = 0;

            if(BIT == 2) // 2 can't be judged by Miller-Rabin, should be specially treated
            {
                for(int prime_cnt = 0; prime_cnt < PRIME_CNT; )
                {
                    int op = rand() & 1;
                    if(op == 0)
                    {
                        N = Trans(2);
                    }
                    else
                    {
                        N = Trans(3);
                    }
                    prime_cnt++;
                    PrintHex(fp, N);
                    fprintf(fp, "\n");
                    fprintf(stdout, "Prime #%d: ", prime_cnt);
                    PrintHex(stdout, N);
                    fprintf(stdout, "\n");
                }
            }
            else
            {
                Init();
                for(int prime_cnt = 0; prime_cnt < PRIME_CNT; )
                {
                    total_tries++;
                    N = RandPrimeLike();
                    int now = MillerRabin(0);
                    if(now == 0) // failed to generate a prime
                    {
                        continue;
                    }
                    // successfully generated a prime
                    prime_cnt++;
                    PrintHex(fp, N);
                    fprintf(fp, "\n");
                    fprintf(stdout, "Prime #%d: ", prime_cnt);
                    PrintHex(stdout, N);
                    fprintf(stdout, "\n");
                }
            }

            end = clock(); 
            total_time = ((double)(end - start)) / CLOCKS_PER_SEC;

            printf("\nTotal Time: %.3fs\n", total_time);
            printf("Avg: %.3fs per prime; %d prime(s) generated\n", total_time / PRIME_CNT, PRIME_CNT);
            if(BIT != 2)
            {
                printf("Avg: %.3fs per try; %d try/tries in total\n\n\n", total_time / total_tries, total_tries);
            }
            else
            {
                printf("\n\n");
            }

            fclose(fp);
        }
        else if(input[0] == 'e' && input[1] == 'n' && input[2] == 'd' && input[3] == '\0')
        {
            printf("\nProgram ended. Thanks for your using.\n");
            break;
        }
        else
        {
            printf("ERROR: No such command\n\n");
        }
    }

    return 0;
}
