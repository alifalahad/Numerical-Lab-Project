# Numerical Lab Project
**CSE 2208 - Numerical Methods Implementation**

---

# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)

- [Newton Interpolation](#newton-interpolation)
  - [Forward Interpolation](#forward-interpolation)
    - [Theory](#forward-interpolation-theory)
    - [Code](#forward-interpolation-code)
    - [Input](#forward-interpolation-input)
    - [Output](#forward-interpolation-output)
  - [Backward Interpolation](#backward-interpolation)
    - [Theory](#backward-interpolation-theory)
    - [Code](#backward-interpolation-code)
    - [Input](#backward-interpolation-input)
    - [Output](#backward-interpolation-output)
  - [Divided Difference](#divided-difference)
    - [Theory](#divided-difference-theory)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)

- [Curve Fitting](#curve-fitting)
  - [Linear Regression](#linear-regression)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)
  - [Polynomial Regression](#polynomial-regression)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)
  - [Transcendental Regression](#transcendental-regression)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)

- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule](#simpsons-13-rule)
    - [Theory](#simpsons-13-theory)
    - [Code](#simpsons-13-code)
    - [Input](#simpsons-13-input)
    - [Output](#simpsons-13-output)
  - [Simpson's 3/8 Rule](#simpsons-38-rule)
    - [Theory](#simpsons-38-theory)
    - [Code](#simpsons-38-code)
    - [Input](#simpsons-38-input)
    - [Output](#simpsons-38-output)

- [Direct Differentiation](#direct-differentiation)
  - [Forward Differentiation](#forward-differentiation)
    - [Theory](#forward-differentiation-theory)
    - [Code](#forward-differentiation-code)
    - [Input](#forward-differentiation-input)
    - [Output](#forward-differentiation-output)
  - [Backward Differentiation](#backward-differentiation)
    - [Theory](#backward-differentiation-theory)
    - [Code](#backward-differentiation-code)
    - [Input](#backward-differentiation-input)
    - [Output](#backward-differentiation-output)

- [Solution of Differential Equations](#solution-of-differential-equations)
  - [Runge Kutta Method](#runge-kutta-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)

---

## Solution of Linear Equations

### Gauss Elimination Method

#### <a name="gauss-elimination-theory"></a>Theory
[Add theory content here]

#### <a name="gauss-elimination-code"></a>Code
```cpp
// View the code file here:
#include <bits/stdc++.h>
using namespace std;
#define ld long double
#define int long long
#define EPS (ld)(1e-12)
bool gaussianElimination(vector<vector<ld> > &A, vector<ld> &B, vector<ld> &x)
{
    int n = A.size();
    // [A | B]
    // Forward elimination
    // Echeleon Form
    for(int i = 0; i < n; i++)
    {
        // Pivoting : find row with max A[row][i]
        int pivot = i;
        for(int row = i + 1; row < n; row++)
        {
            if(fabs(A[row][i]) > fabs(A[pivot][i])) pivot = row;
        }
        if(pivot != i)
        {
            swap(A[i], A[pivot]);
            swap(B[i], B[pivot]);
        }
        // check for zero pivot
        if(fabs(A[i][i]) < EPS) return false;
        // eleminate rows below
        for(int row = i + 1; row < n; row++)
        {
            ld factor = A[row][i] / A[i][i];
            for(int col = i; col < n; col++) A[row][col] -= factor * A[i][col];
            B[row] -= factor * B[i];
        }
    }

    // finding solution
    x.resize(n, 0.0);
    for(int i = n - 1; i >= 0; i--)
    {
        ld sum = B[i];
        for(int j = i + 1; j < n; j++) sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
    }
    return true;
}
void fun()
{
    int n;

    cin >> n;
    vector<vector<ld> > A(n, vector<ld> (n));
    vector<ld> B(n);

    for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) cin >> A[i][j];

    for(int i = 0; i < n; i++) cin >> B[i];
    vector<ld> x;

    // bool check = gaussianElimination(A, B, x);
    bool check = gaussianElimination(A, B, x);

    if(check)
    {
        cout << fixed << setprecision(6) << "Solution is :\n";
        for(auto it : x) cout << it << '\n';
        cout << "Augmented Mstrix:\n";
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++) cout << A[i][j] << ' ';
            cout << B[i] << '\n';
        }
    }
    else
    {
        for(int k = 0; k < n; k++)
        {
            for(int i = k + 1; i < n; i++)
            {
                ld factor = B[k] / B[i];
                for(int j = 0; j < n; j++)
                {
                    if(factor != A[k][j] / A[i][j])
                    {
                        cout << "No solution\n";
                        return;
                    }
                }
            }
        }
        cout << "Infinitely many solutions\n";
    }
}
int32_t main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    int test;

    cin >> test;
    for(int i = 1; i <= test; i++)
    {
        cout << "\nTest case : " << i << "...\n\n";
        fun();
    }
}
```
[Open Gauss_Elimination.cpp](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_ELIMINATION/Gauss_Elimination.cpp)

#### <a name="gauss-elimination-input"></a>Input
```
[Add input format/example here]
3

2
1 2
2 4
3
6

5
1 2 3 4 5
5 6 7 8 9
6 5 4 3 2
3 4 5 43 2
0 9 2 1 3
1
2
3
4
5

3
1 2 3
4 5 6
7 8 8
1
2
4
```
[Open input.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_ELIMINATION/input.txt)

#### <a name="gauss-elimination-output"></a>Output
```
[Add output format/example here]

Test case : 1...

Infinitely many solutions

Test case : 2...

No solution

Test case : 3...

Solution is :
-1.333333
2.666667
-1.000000
Augmented Mstrix:
7.000000 8.000000 8.000000 4.000000
0.000000 0.857143 1.857143 0.428571
0.000000 0.000000 0.500000 -0.500000
```
[Open output.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_ELIMINATION/output.txt)

---

### Gauss Jordan Elimination Method

#### <a name="gauss-jordan-theory"></a>Theory
[Add theory content here]

#### <a name="gauss-jordan-code"></a>Code
```cpp
// View the code file here:  
```
[Open Gauss_Jordan.cpp](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/Gauss_Jordan.cpp)

#### <a name="gauss-jordan-input"></a>Input
```
[Add input format/example here]
```
[Open input.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/input.txt)

#### <a name="gauss-jordan-output"></a>Output
```
[Add output format/example here]
```
[Open output.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/output.txt)

---

### LU Decomposition Method

#### <a name="lu-decomposition-theory"></a>Theory
[Add theory content here]

#### <a name="lu-decomposition-code"></a>Code
```cpp
// View the code file here:   
```
[Open LU_DECOMPOSITION.cpp](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/LU_DECOMPOSITION/LU_DECOMPOSITION.cpp)

#### <a name="lu-decomposition-input"></a>Input
```
[Add input format/example here]
```
[Open input.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/input.txt)

#### <a name="lu-decomposition-output"></a>Output
```
[Add output format/example here]
```
[Open output.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/output. txt)

---

### Matrix Inversion

#### <a name="matrix-inversion-theory"></a>Theory
[Add theory content here]

#### <a name="matrix-inversion-code"></a>Code
```cpp
// View the code file here:   
```
[Open Matrix_Inversion.cpp](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/Matrix_Inversion. cpp)

#### <a name="matrix-inversion-input"></a>Input
```
[Add input format/example here]
```
[Open input.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/input.txt)

#### <a name="matrix-inversion-output"></a>Output
```
[Add output format/example here]
```
[Open output.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/output. txt)

---

## Solution of Non-Linear Equations

### Bisection Method

#### <a name="bisection-theory"></a>Theory

**Overview**

The Bisection Method is a reliable numerical technique for finding roots of continuous functions. It operates on the principle of the Intermediate Value Theorem:  if a continuous function f(x) changes sign over an interval [a, b], then there exists at least one root within that interval.

**Mathematical Foundation**

For a root to exist in interval [a, b]: 
- f(a) × f(b) < 0 (opposite signs)

**Algorithm Steps**

1. **Initialization**: Select interval [a, b] where f(a) × f(b) < 0
2. **Midpoint Calculation**: Compute c = (a + b) / 2
3. **Function Evaluation**: Calculate f(c)
4. **Interval Update**:
   - If f(a) × f(c) < 0 → New interval:  [a, c] (set b = c)
   - If f(c) × f(b) < 0 → New interval:  [c, b] (set a = c)
5. **Convergence Check**:  Repeat until |b - a| < tolerance
6. **Result**: Final midpoint c approximates the root

**Key Characteristics**

✅ **Advantages:**
- Guaranteed convergence for continuous functions
- Simple implementation
- Reliable and robust

❌ **Limitations:**
- Slow convergence (linear rate)
- Requires initial interval with sign change
- Inefficient for multiple or complex roots

#### <a name="bisection-code"></a>Code
```cpp
// View the code file here:  
```
[Open Bisection.cpp](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/Bisection.cpp)

#### <a name="bisection-input"></a>Input
```
[Add input format/example here]
```
[Open input.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/input.txt)

#### <a name="bisection-output"></a>Output
```
[Add output format/example here]
```
[Open output.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/output.txt)

---

### False Position Method

#### <a name="false-position-theory"></a>Theory

**Overview**

The False Position Method (Regula Falsi) is an enhanced bracketing technique that combines bisection with linear interpolation. It typically converges faster than bisection by using a weighted average based on function values.

**Mathematical Formula**

Root approximation: 
```
c = [a × f(b) - b × f(a)] / [f(b) - f(a)]
```

Where:
- [a, b] is the bracketing interval
- f(a) · f(b) < 0

**Algorithm**

1. **Input**: Degree, coefficients, initial guesses a and b, tolerance E
2. **Validation**:  Verify f(a) · f(b) < 0
3. **Compute**: c using linear interpolation formula
4. **Update Interval**:
   - If f(a) · f(c) < 0 → set b = c
   - Otherwise → set a = c
5. **Iterate**: Repeat until |f(c)| ≤ E
6. **Output**: Root value and iteration count

**Implementation Features**

- Uses `long double` for enhanced precision
- File-based I/O (input. txt, output.txt)
- Displays iteration table with a, b, c, and f(c) values
- Polynomial evaluation via coefficient vector

**Key Characteristics**

✅ **Advantages:**
- Faster than bisection method
- Guaranteed convergence with valid bracket
- Good accuracy for polynomial equations

❌ **Limitations:**
- Can be slow when one endpoint is fixed
- Requires sign change in initial interval
- Not ideal for multiple roots

#### <a name="false-position-code"></a>Code
```cpp
// View the code file here:
#include <bits/stdc++.h>
#define ld long double
using namespace std;
ld f(ld x, const vector<ld> &cof) {
    ld sum = 0.0;
    int n = cof.size();
    for (int i = 0; i < n; i++)
        sum += cof[i] * powl(x, n - 1 - i);
    return sum;
}
int main() {
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    int test;
    cin >> test;
    for (int tt = 1; tt <= test; tt++) {
        cout << "Test Case : " << tt << "...\n\n";
        int degree;
        cin >> degree;
        vector<ld> cof(degree + 1);
        for (int i = 0; i <= degree; i++)
            cin >> cof[i];
        ld a, b, E;
        cin >> a >> b >> E;
        ld fa = f(a, cof);
        ld fb = f(b, cof);
        if (fa * fb >= 0) {
            cout << "Invalid initial guesses. f(a) and f(b) must have opposite signs.\n\n";
            continue;
        }
        cout << setw(10) << "Iteration"<< setw(15) << "a"<< setw(15) << "b"<< setw(15) << "c"<< setw(20) << "f(c)" << endl;
        ld c, fc;
        int iter = 0;
        do {
            c = (a * fb - b * fa) / (fb - fa);
            fc = f(c, cof);
            cout << setw(10) << ++iter<< setw(15) << fixed << setprecision(6) << a<< setw(15) << b<< setw(15) << c<< setw(20) << fc << endl;
            if (fa * fc < 0) {
                b = c;
                fb = fc;
            } else {
                a = c;
                fa = fc;
            }
        } while (fabsl(fc) > E);
        cout << "\nRoot is approximately: " << c << endl;
        cout << "Total iterations: " << iter << "\n\n";
    }
    return 0;
}
```
[Open False_Position.cpp](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/FALSE_POSITION/False_Position.cpp)

#### <a name="false-position-input"></a>Input
```
[Add input format/example here]
3
2
1 -4 -10
1
3
0.0001

3
1 -6 11 -6
1
2
0.00001

4
1 0 -10 0 9
0
2
0.00001

```
[Open input.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/FALSE_POSITION/input.txt)

#### <a name="false-position-output"></a>Output
```
[Add output format/example here]
Test Case : 1...

Invalid initial guesses. f(a) and f(b) must have opposite signs.

Test Case : 2...

Invalid initial guesses. f(a) and f(b) must have opposite signs.

Test Case : 3...

 Iteration              a              b              c                f(c)
         1       0.000000       2.000000       0.750000            3.691406
         2       0.750000       2.000000       0.996865            0.050117
         3       0.996865       2.000000       1.000206           -0.003291
         4       0.996865       1.000206       1.000000            0.000003

Root is approximately: 1.000000
Total iterations: 4

```
[Open output.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/FALSE_POSITION/output.txt)

---

### Newton Raphson Method

#### <a name="newton-raphson-theory"></a>Theory

**Overview**

The Newton-Raphson Method is one of the most powerful and widely-used iterative techniques for solving nonlinear equations f(x) = 0. It leverages the tangent line at an initial guess to rapidly approximate the actual root, offering quadratic convergence under suitable conditions.

**Mathematical Principle**

The method uses linear approximation via the derivative: 

```
x_(n+1) = x_n - f(x_n) / f'(x_n)
```

Where:
- x_n is the current approximation
- f(x_n) is the function value
- f'(x_n) is the derivative at x_n

**Algorithm**

1. **Initialize**: Choose x₀ near the expected root
2. **Compute**:  Evaluate f(x₀) and f'(x₀)
3. **Update**: Calculate x₁ = x₀ - [f(x₀) / f'(x₀)]
4. **Check Convergence**: If |x₁ - x₀| < ε, stop
5. **Iterate**: Set x₀ = x₁ and repeat steps 2-4
6. **Result**: Final x_n is the approximate root

**Convergence Behavior**

- **Quadratic convergence** near the root (very fast)
- Requires good initial guess
- Derivative must be non-zero

**Key Characteristics**

✅ **Advantages:**
- Extremely fast convergence (quadratic)
- Efficient for well-behaved functions
- Fewer iterations than bracketing methods

❌ **Limitations:**
- Requires derivative calculation
- May diverge with poor initial guess
- Fails when f'(x) = 0

#### <a name="newton-raphson-code"></a>Code
```cpp
// View the code file here: 
```
[Open Newton_Raphson.cpp](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/Newton_Raphson.cpp)

#### <a name="newton-raphson-input"></a>Input
```
[Add input format/example here]
```
[Open input_nr.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/input_nr.txt)

#### <a name="newton-raphson-output"></a>Output
```
[Add output format/example here]
```
[Open output_nr.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/output_nr.txt)

---

### Secant Method

#### <a name="secant-theory"></a>Theory
[Add theory content here]

#### <a name="secant-code"></a>Code
```cpp
// View the code file here: 
```
[Open Secant.cpp](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/Secant.cpp)

#### <a name="secant-input"></a>Input
```
[Add input format/example here]
```
[Open input_sec.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/input_sec.txt)

#### <a name="secant-output"></a>Output
```
[Add output format/example here]
```
[Open output_sec.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/output_sec.txt)

---

## Newton Interpolation

### Forward Interpolation

#### <a name="forward-interpolation-theory"></a>Theory
[Add theory content here]

#### <a name="forward-interpolation-code"></a>Code
```cpp
// View the code file here:  
```
[Open newton_forward.cpp](./src/Newton%20interpolation/Forward/newton_forward.cpp)

#### <a name="forward-interpolation-input"></a>Input
```
[Add input format/example here]
```
[Open input_fw.txt](./src/Newton%20interpolation/Forward/input_fw.txt)

#### <a name="forward-interpolation-output"></a>Output
```
[Add output format/example here]
```
[Open output_fw.txt](./src/Newton%20interpolation/Forward/output_fw. txt)

---

### Backward Interpolation

#### <a name="backward-interpolation-theory"></a>Theory
[Add theory content here]

#### <a name="backward-interpolation-code"></a>Code
```cpp
// View the code file here:  
```
[Open newton_backward.cpp](./src/Newton%20interpolation/Backward/newton_backward.cpp)

#### <a name="backward-interpolation-input"></a>Input
```
[Add input format/example here]
```
[Open input_bw.txt](./src/Newton%20interpolation/Backward/input_bw.txt)

#### <a name="backward-interpolation-output"></a>Output
```
[Add output format/example here]
```
[Open output_bw. txt](./src/Newton%20interpolation/Backward/output_bw.txt)

---

### Divided Difference

#### <a name="divided-difference-theory"></a>Theory
[Add theory content here]

#### <a name="divided-difference-code"></a>Code
```cpp
// View the code file here: 
```
[Open div. cpp](./src/Newton%20interpolation/divided%20difference/div.cpp)

#### <a name="divided-difference-input"></a>Input
```
[Add input format/example here]
```
[Open input_div.txt](./src/Newton%20interpolation/divided%20difference/input_div.txt)

#### <a name="divided-difference-output"></a>Output
```
[Add output format/example here]
```
[Open output_div.txt](./src/Newton%20interpolation/divided%20difference/output_div.txt)

---

## Curve Fitting

### Linear Regression

#### <a name="linear-regression-theory"></a>Theory
[Add theory content here]

#### <a name="linear-regression-code"></a>Code
```cpp
// View the code file here:
```
[Open linear. cpp](./src/Curve%20Fitting/regression/Linear/linear.cpp)

#### <a name="linear-regression-input"></a>Input
```
[Add input format/example here]
```
[Open input_lin.txt](./src/Curve%20Fitting/regression/Linear/input_lin.txt)

#### <a name="linear-regression-output"></a>Output
```
[Add output format/example here]
```
[Open output_lin.txt](./src/Curve%20Fitting/regression/Linear/output_lin.txt)

---

### Polynomial Regression

#### <a name="polynomial-regression-theory"></a>Theory
[Add theory content here]

#### <a name="polynomial-regression-code"></a>Code
```cpp
// View the code file here:  
```
[Open polynomial.cpp](./src/Curve%20Fitting/regression/polynomial/polynomial.cpp)

#### <a name="polynomial-regression-input"></a>Input
```
[Add input format/example here]
```
[Open input_pol.txt](./src/Curve%20Fitting/regression/polynomial/input_pol.txt)

#### <a name="polynomial-regression-output"></a>Output
```
[Add output format/example here]
```
[Open output_pol.txt](./src/Curve%20Fitting/regression/polynomial/output_pol.txt)

---

### Transcendental Regression

#### <a name="transcendental-regression-theory"></a>Theory
[Add theory content here]

#### <a name="transcendental-regression-code"></a>Code
```cpp
// View the code file here:
```
[Open regression_transcendental.cpp](./src/Curve%20Fitting/regression/trascendental/regression_transcendental.cpp)

#### <a name="transcendental-regression-input"></a>Input
```
[Add input format/example here]
```
[Open input_trans.txt](./src/Curve%20Fitting/regression/trascendental/input_trans.txt)

#### <a name="transcendental-regression-output"></a>Output
```
[Add output format/example here]
```
[Open output_trans.txt](./src/Curve%20Fitting/regression/trascendental/output_trans.txt)

---

## Numerical Integration

### Simpson's 1/3 Rule

#### <a name="simpsons-13-theory"></a>Theory
[Add theory content here]

#### <a name="simpsons-13-code"></a>Code
```cpp
// View the code file here:  
```
[Open one_third.cpp](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/one%20_third. cpp)

#### <a name="simpsons-13-input"></a>Input
```
[Add input format/example here]
```
[Open input1_3.txt](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/input1_3.txt)

#### <a name="simpsons-13-output"></a>Output
```
[Add output format/example here]
```
[Open output1_3.txt](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/output1_3.txt)

---

### Simpson's 3/8 Rule

#### <a name="simpsons-38-theory"></a>Theory
[Add theory content here]

#### <a name="simpsons-38-code"></a>Code
```cpp
// View the code file here:  
```
[Open newton3_8.cpp](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/newton3_8.cpp)

#### <a name="simpsons-38-input"></a>Input
```
[Add input format/example here]
```
[Open input_3_8.txt](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/input_3_8.txt)

#### <a name="simpsons-38-output"></a>Output
```
[Add output format/example here]
```
[Open output_3_8.txt](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/output_3_8.txt)

---

## Direct Differentiation

### Forward Differentiation

#### <a name="forward-differentiation-theory"></a>Theory
[Add theory content here]

#### <a name="forward-differentiation-code"></a>Code
```cpp
// View the code file here:  
```
[Open Differentiation_FI.cpp](./src/Direct%20Differentiation/forward/Differentiation_FI.cpp)

#### <a name="forward-differentiation-input"></a>Input
```
[Add input format/example here]
```
[Open input_diff_fwd.txt](./src/Direct%20Differentiation/forward/input_diff_fwd.txt)

#### <a name="forward-differentiation-output"></a>Output
```
[Add output format/example here]
```
[Open output_diff_fwd. txt](./src/Direct%20Differentiation/forward/output_diff_fwd.txt)

---

### Backward Differentiation

#### <a name="backward-differentiation-theory"></a>Theory
[Add theory content here]

#### <a name="backward-differentiation-code"></a>Code
```cpp
// View the code file here: 
```
[Open differentiation_BI.cpp](./src/Direct%20Differentiation/backward/differentiation_BI.cpp)

#### <a name="backward-differentiation-input"></a>Input
```
[Add input format/example here]
```
[Open input_diff. txt](./src/Direct%20Differentiation/backward/input_diff.txt)

#### <a name="backward-differentiation-output"></a>Output
```
[Add output format/example here]
```
[Open output_diff.txt](./src/Direct%20Differentiation/backward/output_diff.txt)

---

## Solution of Differential Equations

### Runge Kutta Method

#### <a name="runge-kutta-theory"></a>Theory
[Add theory content here]

#### <a name="runge-kutta-code"></a>Code
```cpp
// View the code file here: 
```
[Open RungeKutta4th.cpp](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/RungeKutta4th.cpp)

#### <a name="runge-kutta-input"></a>Input
```
[Add input format/example here]
```
[Open input. txt](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/input.txt)

#### <a name="runge-kutta-output"></a>Output
```
[Add output format/example here]
```
[Open output.txt](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/output.txt)

---

## Contributing
This is a course project for CSE 2208.  Feel free to explore the implementations and learn from the numerical methods.  

## License
Educational purposes only. 
