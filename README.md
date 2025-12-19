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
Introduction
The Gauss Elimination method is a direct numerical technique used to solve systems of linear equations of the form Ax = b.
Theory: The method is theoretically based on the concept of row equivalence. It relies on the principle that applying elementary row operations (such as swapping rows, multiplying a row by a non-zero constant, or adding a multiple of one row to another) preserves the solution set of the system. By systematically applying these operations, the method transforms the dense coefficient matrix A into an upper triangular matrix (row-echelon form). Once in this simplified form, the unknowns can be solved sequentially without complex matrix inversion.
Algorithm
The process consists of two main stages: Forward Elimination and Back Substitution.
Stage 1: Forward Elimination
1.	Form the Augmented Matrix: Combine the coefficient matrix and the constant vector into a single matrix [A | b].
2.	Pivot Selection: Start with the first row. The first element is the "pivot."
3.	Eliminate Entries Below Pivot: Use row operations to create zeros in the column below the current pivot.
o	Formula: Row_i = Row_i – (factor) × Pivot_Row
o	Factor: (Element to eliminate) / (Pivot element)
4.	Iterate: Move to the next diagonal element (next pivot) and repeat the process until the matrix is in upper triangular form (all elements below the main diagonal are zero).
Stage 2: Back Substitution
1.	Solve for the Last Variable: The last row now contains an equation with only one variable. Solve for it directly.
2.	Substitute Upwards: Substitute the known value into the row above it to solve for the next variable.
3.	Repeat: Continue this process moving upwards until all variables (x₁, x₂, ..., xₙ) are found.
Advantages
•	Generality: It can solve any system of n linear equations with n unknowns, provided a unique solution exists.
•	Systematic Approach: The method is algorithmic and easy to program for computers.
•	Exactness: Theoretically, it produces the exact solution (ignoring computer round-off errors) unlike iterative methods that produce approximations.
Disadvantages
•	Computational Cost: It is computationally expensive for very large systems, with a time complexity of approximately O(n³).
•	Round-off Errors: In computer implementation, the repeated arithmetic operations can accumulate round-off errors, leading to inaccurate results for ill-conditioned matrices.
•	Division by Zero: If a pivot element is zero, the method fails unless row swapping (partial pivoting) is implemented.



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
#include <bits/stdc++.h>
using namespace std;
#define ld long double
#define int long long
#define EPS (ld)(1e-12)
bool gaussJordan(vector<vector<ld> > &A, vector<ld> &B, vector<ld> &x)
{
    int n = A.size();
    // [A | B]
    // Forward elimination
    // Canonical form
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
        // make pivot = 1
        ld pivotval = A[i][i];
        for(int j = 0; j < n; j++) A[i][j] /= pivotval;
        B[i] /= pivotval;
        // Eliminate all other entries in this column
        for(int row = 0; row < n; row++)
        {
            if(row == i) continue;
            ld factor = A[row][i];
            for(int col = 0; col < n; col++) A[row][col] -= factor * A[i][col];
            B[row] -= factor * B[i];
        }
    }

    // finding solution
    x = B;
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
    bool check = gaussJordan(A, B, x);

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
[Open Gauss_Jordan.cpp](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/Gauss_Jordan.cpp)

#### <a name="gauss-jordan-input"></a>Input
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
[Open input.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/input.txt)

#### <a name="gauss-jordan-output"></a>Output
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
1.000000 0.000000 0.000000 -1.333333
0.000000 1.000000 0.000000 2.666667
0.000000 0.000000 1.000000 -1.000000
```
[Open output.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/output.txt)

---

### LU Decomposition Method

#### <a name="lu-decomposition-theory"></a>Theory
[Add theory content here]

#### <a name="lu-decomposition-code"></a>Code
```cpp
// View the code file here:
#include <bits/stdc++.h>
using namespace std;
void displayMatrix(const vector<vector<double>> &matrix, string matrixName) {
    cout << "\nCurrent state of " << matrixName << " matrix:\n";
    for (auto &row : matrix) {
        for (auto val : row)
            cout << setw(10) << fixed << setprecision(4) << val << " ";
        cout << "\n";
    }
    cout << "---------------------------------------------\n";
}
int main() {
    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);
    int test;
    cin >> test;
    for(int tt = 1; tt <= test; tt++)
    {
        cout << "Test Case : " << tt << "...\n\n";
        int n;
        cout << "Enter number of equations: ";
        cin >> n;
        vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));
        cout << "Enter the augmented matrix (coefficients + constants):\n";
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= n; j++)
                cin >> augmentedMatrix[i][j];
        // Separate coefficient matrix and constants
        vector<vector<double>> coeffMatrix(n, vector<double>(n));
        vector<double> constants(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                coeffMatrix[i][j] = augmentedMatrix[i][j];
            constants[i] = augmentedMatrix[i][n];
        }
        vector<vector<double>> lowerTri(n, vector<double>(n, 0));
        vector<vector<double>> upperTri(n, vector<double>(n, 0));
        cout << "\nPerforming LU Decomposition...\n";
        for (int i = 0; i < n; i++) {
            // Compute upper triangular matrix
            for (int j = i; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++)
                    sum += lowerTri[i][k] * upperTri[k][j];
                upperTri[i][j] = coeffMatrix[i][j] - sum;
            }
            // Compute lower triangular matrix
            for (int j = i; j < n; j++) {
                if (i == j)
                    lowerTri[i][i] = 1; // diagonal element
                else {
                    double sum = 0;
                    for (int k = 0; k < i; k++)
                        sum += lowerTri[j][k] * upperTri[k][i];

                    if (fabs(upperTri[i][i]) < 1e-9) {
                        cout << "\nZero pivot encountered — system may have no or infinite solutions.\n";
                        return 0;
                    }
                    lowerTri[j][i] = (coeffMatrix[j][i] - sum) / upperTri[i][i];
                }
            }
            displayMatrix(lowerTri, "L");
            displayMatrix(upperTri, "U");
        }
        // Check determinant (product of diagonal elements of U)
        double determinant = 1;
        for (int i = 0; i < n; i++)
            determinant *= upperTri[i][i];
        if (fabs(determinant) < 1e-9) {
            cout << "\nDeterminant = 0 → No unique solution (infinite or none).\n";
            return 0;
        }
        cout << "\nLU Decomposition completed successfully!\n";
        displayMatrix(lowerTri, "Final L");
        displayMatrix(upperTri, "Final U");
        // Forward substitution: L*y = constants
        vector<double> intermediate(n);
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += lowerTri[i][j] * intermediate[j];
            intermediate[i] = constants[i] - sum;
        }
        cout << "\nForward Substitution (L*y = b):\n";
        for (int i = 0; i < n; i++)
            cout << "y" << i + 1 << " = " << intermediate[i] << "\n";
        // Back substitution: U*x = intermediate
        vector<double> solution(n);
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < n; j++)
                sum += upperTri[i][j] * solution[j];
            solution[i] = (intermediate[i] - sum) / upperTri[i][i];
        }
        cout << "\nBack Substitution (U*x = y):\n";
        cout << "Solution Vector (x):\n";
        for (int i = 0; i < n; i++)
            cout << "x" << i + 1 << " = " << solution[i] << "\n";
    }
    return 0;
}

```
[Open LU_DECOMPOSITION.cpp](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/LU_DECOMPOSITION/LU_DECOMPOSITION.cpp)

#### <a name="lu-decomposition-input"></a>Input
```
[Add input format/example here]
3

3
2 1 -1 8 
-3 -1 2 -11 
-2 1 2 -3

2
4 3 10 
2 1 4 

2
2 4 6
1 2 3



```
[Open input.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/input.txt)

#### <a name="lu-decomposition-output"></a>Output
```
[Add output format/example here]
Test Case : 1...

Enter number of equations: Enter the augmented matrix (coefficients + constants):

Performing LU Decomposition...

Current state of L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     0.0000     0.0000 
   -1.0000     0.0000     0.0000 
---------------------------------------------

Current state of U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.0000     0.0000 
    0.0000     0.0000     0.0000 
---------------------------------------------

Current state of L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     1.0000     0.0000 
   -1.0000     4.0000     0.0000 
---------------------------------------------

Current state of U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     0.0000     0.0000 
---------------------------------------------

Current state of L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     1.0000     0.0000 
   -1.0000     4.0000     1.0000 
---------------------------------------------

Current state of U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     0.0000    -1.0000 
---------------------------------------------

LU Decomposition completed successfully!

Current state of Final L matrix:
    1.0000     0.0000     0.0000 
   -1.5000     1.0000     0.0000 
   -1.0000     4.0000     1.0000 
---------------------------------------------

Current state of Final U matrix:
    2.0000     1.0000    -1.0000 
    0.0000     0.5000     0.5000 
    0.0000     0.0000    -1.0000 
---------------------------------------------

Forward Substitution (L*y = b):
y1 = 8.0000
y2 = 1.0000
y3 = 1.0000

Back Substitution (U*x = y):
Solution Vector (x):
x1 = 2.0000
x2 = 3.0000
x3 = -1.0000
Test Case : 2...

Enter number of equations: Enter the augmented matrix (coefficients + constants):

Performing LU Decomposition...

Current state of L matrix:
    1.0000     0.0000 
    0.5000     0.0000 
---------------------------------------------

Current state of U matrix:
    4.0000     3.0000 
    0.0000     0.0000 
---------------------------------------------

Current state of L matrix:
    1.0000     0.0000 
    0.5000     1.0000 
---------------------------------------------

Current state of U matrix:
    4.0000     3.0000 
    0.0000    -0.5000 
---------------------------------------------

LU Decomposition completed successfully!

Current state of Final L matrix:
    1.0000     0.0000 
    0.5000     1.0000 
---------------------------------------------

Current state of Final U matrix:
    4.0000     3.0000 
    0.0000    -0.5000 
---------------------------------------------

Forward Substitution (L*y = b):
y1 = 10.0000
y2 = -1.0000

Back Substitution (U*x = y):
Solution Vector (x):
x1 = 1.0000
x2 = 2.0000
Test Case : 3...

Enter number of equations: Enter the augmented matrix (coefficients + constants):

Performing LU Decomposition...

Current state of L matrix:
    1.0000     0.0000 
    0.5000     0.0000 
---------------------------------------------

Current state of U matrix:
    2.0000     4.0000 
    0.0000     0.0000 
---------------------------------------------

Current state of L matrix:
    1.0000     0.0000 
    0.5000     1.0000 
---------------------------------------------

Current state of U matrix:
    2.0000     4.0000 
    0.0000     0.0000 
---------------------------------------------

Determinant = 0 → No unique solution (infinite or none).
```
[Open output.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/output. txt)

---

### Matrix Inversion

#### <a name="matrix-inversion-theory"></a>Theory
[Add theory content here]

#### <a name="matrix-inversion-code"></a>Code
```cpp
// View the code file here:
#include <bits/stdc++.h>
using namespace std;
// Function to get cofactor of A[p][q]
void getCofactor(vector<vector<double>> &A, vector<vector<double>> &temp,int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = A[row][col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
// Determinant function
double determinant(vector<vector<double>> &A, int n) {
    if (n == 1)
        return A[0][0];
    double det = 0;
    vector<vector<double>> temp(n, vector<double>(n));
    int sign = 1;
    for (int f = 0; f < n; f++) {
        getCofactor(A, temp, 0, f, n);
        det += sign * A[0][f] * determinant(temp, n - 1);
        sign = -sign;
    }
    return det;
}
// Adjoint = transpose of cofactor matrix
void adjoint(vector<vector<double>> &A, vector<vector<double>> &adj, int n) {
    if (n == 1) {
        adj[0][0] = 1;
        return;
    }
    int sign;
    vector<vector<double>> temp(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            getCofactor(A, temp, i, j, n);
            sign = ((i + j) % 2 == 0) ? 1 : -1;
            // adjoint = transpose of cofactor matrix
            adj[j][i] = sign * determinant(temp, n - 1);
        }
    }
}
// Inverse = adj(A) / det(A)
bool inverse(vector<vector<double>> &A, vector<vector<double>> &inv, int n) {
    double det = determinant(A, n);
    if (det == 0) {
        cout << "Matrix is singular so No inverse exists.\n";
        return false;
    }
    vector<vector<double>> adj(n, vector<double>(n));
    adjoint(A, adj, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[i][j] / det;
    return true;
}
int main() {
    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);
    int test;
    cin >> test;
    for(int tt = 1; tt <= test; tt++)
    {
        cout << "Test Case : " << tt << "...\n\n";
        int n;
        cin >> n;
        vector<vector<double>> A(n, vector<double>(n));
        vector<vector<double>> inv(n, vector<double>(n));
        vector<double> B(n), X(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                cin >> A[i][j];
        for (int i = 0; i < n; i++)
            cin >> B[i];
        if (inverse(A, inv, n)) {
            cout << "\nInverse of matrix A:\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++)
                    cout << fixed << setprecision(4) << inv[i][j] << " ";
                cout << endl;
            }
            // Compute X = A^-1 * B
            for (int i = 0; i < n; i++) {
                X[i] = 0;
                for (int j = 0; j < n; j++)
                    X[i] += inv[i][j] * B[j];
            }
            cout << "\nRoots :\n";
            for (int i = 0; i < n; i++)
                cout << "x" << i + 1 << " = " << fixed << setprecision(4) << X[i] << endl;
        }
    }
    return 0;
}
```
[Open Matrix_Inversion.cpp](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/Matrix_Inversion. cpp)

#### <a name="matrix-inversion-input"></a>Input
```
[Add input format/example here]
3

3
1 1 1
2 5 3
1 2 4
6
23
17

2
2 1
1 1
5
3

2
1 1
2 2
2
4


```
[Open input.txt](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/input.txt)

#### <a name="matrix-inversion-output"></a>Output
```
[Add output format/example here]
Test Case : 1...


Inverse of matrix A:
1.7500 -0.2500 -0.2500 
-0.6250 0.3750 -0.1250 
-0.1250 -0.1250 0.3750 

Roots :
x1 = 0.5000
x2 = 2.7500
x3 = 2.7500
Test Case : 2...


Inverse of matrix A:
1.0000 -1.0000 
-1.0000 2.0000 

Roots :
x1 = 2.0000
x2 = 1.0000
Test Case : 3...

Matrix is singular so No inverse exists.
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
#include <bits/stdc++.h>
using namespace std;


int degree;
vector<double> coeffs;

double f(double x) {
    double result = coeffs[0];
    for (int i = 1; i <= degree; i++) {
        result = result * x + coeffs[i];
    }
    return result;
}

double f_prime(double x) {
    double result = degree * coeffs[0];
    for (int i = 1; i < degree; i++) {
        result = result * x + (degree - i) * coeffs[i];
    }
    return result;
}

int main() {
    freopen("input_nr.txt" , "r" , stdin);
    freopen("output_nr.txt" , "w" ,stdout);
    cout << "Newton-Raphson Method for Root Finding" << endl;

    
    //cout << "\nEnter the degree of polynomial: ";
    cin >> degree;
    
    if (degree < 1) {
        cout << "Error: Degree must be at least 1." << endl;
        return 1;
    }
    
    coeffs.resize(degree + 1);
    //cout << "Enter " << (degree + 1) << " coefficients (from a_" << degree << " to a_0):" << endl;
    for (int i = 0; i <= degree; i++) {
        //cout << "Coefficient a_" << (degree - i) << ": ";
        cin >> coeffs[i];
    }
    
    double x0;       
    double tol;     
    int max_iter;   
    //cout << "\nEnter initial guess: ";
    cin >> x0;
    //cout << "Enter tolerance: ";
    cin >> tol;
    //cout << "Enter maximum iterations: ";
    cin >> max_iter;

    double x1;
    int iter = 0;


    while (iter < max_iter) {
        double fx = f(x0);
        double fpx = f_prime(x0);
        
        if (fabs(fpx) < 1e-10) {
            cout << "Error: Derivative too close to zero. Method fails." << endl;
            return 1;
        }
        
        x1 = x0 - fx / fpx;

        cout << "Iteration " << (iter + 1) << ": x = " << x1 << endl;

        if (fabs(x1 - x0) < tol) {
            break;
        }
        x0 = x1;
        iter++;
    }

    cout << "Root = " << x1 << endl;
    cout << "Iterations = " << iter + 1 << endl;

    return 0;
}
```
[Open Newton_Raphson.cpp](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/Newton_Raphson.cpp)

#### <a name="newton-raphson-input"></a>Input
```
[Add input format/example here]
2
1 0 -4
6 0.0001 10
```
[Open input_nr.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/input_nr.txt)

#### <a name="newton-raphson-output"></a>Output
```
[Add output format/example here]
Newton-Raphson Method for Root Finding
Iteration 1: x = 3.33333
Iteration 2: x = 2.26667
Iteration 3: x = 2.01569
Iteration 4: x = 2.00006
Iteration 5: x = 2
Root = 2
Iterations = 5
```
[Open output_nr.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/output_nr.txt)

---

### Secant Method

#### <a name="secant-theory"></a>Theory
[Add theory content here]
Introduction
The Secant method is an iterative numerical technique used to find the roots of a function f(x) = 0. It is an open method that uses a secant line connecting two points on the curve to approximate the root.
Unlike the Newton-Raphson method, the Secant method does not require the calculation of the derivative f'(x). Instead, it approximates the derivative using finite differences based on two previous points (xₙ₋₁ and xₙ). This makes it particularly useful when differentiating the function is difficult or computationally expensive. The method requires two initial guesses but does not require the root to be bracketed.
Algorithm
1.	Initialize:
o	Choose two initial guesses, x₀ and x₁.
o	Set a tolerance error E and a maximum number of iterations.
2.	Iteration Step: For each iteration, calculate the next approximation xₙ₊₁ using the formula:
x_(n+1) = x_n - f(x_n) * [ (x_n - x_(n-1)) / (f(x_n) - f(x_(n-1))) ]
3.	Check Convergence:
o	Calculate the absolute error: |xₙ₊₁ – xₙ|.
o	If the error is ≤ E, the root is found.
4.	Update:
o	Set xₙ₋₁ = xₙ
o	Set xₙ = xₙ₊₁
o	Repeat Step 2 until convergence or maximum iterations are reached.
Advantages
•	No Derivative Needed: Unlike Newton-Raphson, it does not require analytical differentiation, making it suitable for complex functions.
•	Faster than Bisection: It has a superlinear convergence rate (approx. 1.618), which is significantly faster than the linear convergence of the Bisection method.
•	Efficiency: After the first step, it requires only one function evaluation per iteration (reusing the previous value).
Disadvantages
•	No Guaranteed Convergence: Unlike the Bisection method, the Secant method may diverge if the initial guesses are not close enough to the root.
•	Slower than Newton-Raphson: Its convergence rate is slightly slower than Newton’s method (which is quadratic, 2.0).
•	Numerical Instability: If f(xₙ) and f(xₙ₋₁) are nearly equal, the denominator approaches zero, potentially causing large errors or division by zero.



#### <a name="secant-code"></a>Code
```cpp
// View the code file here:
#include <bits/stdc++.h>
using namespace std;
int degree;
vector<double> coeffs;
double f(double x) {
    double result = coeffs[0];
    for (int i = 1; i <= degree; i++) {
        result = result * x + coeffs[i];
    }
    return result;
}
double secant(double x0, double x1, double tol = 0.0001, int maxIter = 100) {
    double x2;
    int iter = 0;
    while(abs(x2 - x1) > tol && iter < maxIter) {
        double f0 = f(x0);
        double f1 = f(x1);
        if(f1 - f0 == 0) return NAN;
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        x0 = x1;
        x1 = x2;
        iter++;
    }
    return x2;
}
int main() {
     freopen("input_sec.txt", "r", stdin);
    freopen("output_sec.txt", "w", stdout);
    cin >> degree;
    coeffs.resize(degree + 1);
    for (int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }
    double minRange, maxRange;
    cin >> minRange >> maxRange;
    double step = 0.1;
    double tol = 0.0001;
    vector<double> roots;
    for(double x = minRange; x < maxRange; x += step) {
        double x0 = x;
        double x1 = x + step;
        if(f(x0) * f(x1) <= 0) {
            double root = secant(x0, x1, tol);
            if(!isnan(root)) {
                roots.push_back(root);
            }
        }
    }
    cout << fixed << setprecision(6);
    for(double r : roots)
        cout << r << endl;
    return 0;
}
```
[Open Secant.cpp](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/Secant.cpp)

#### <a name="secant-input"></a>Input
```
[Add input format/example here]
2
1 0 -3
-2 2
```
[Open input_sec.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/input_sec.txt)

#### <a name="secant-output"></a>Output
```
[Add output format/example here]
-1.731429
1.731429
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
#include <bits/stdc++.h>
using namespace std;
int main() {
    freopen("input_lin.txt" , "r" , stdin);
    freopen("output_lin.txt" , "w" ,stdout);
    int n;
    //cout << "Enter number of data points: ";
    cin >> n;
    vector<double> x(n), y(n);
    //cout << "Enter x values:\n";
    for (int i = 0; i < n; i++) cin >> x[i];
    //cout << "Enter y values:\n";
    for (int i = 0; i < n; i++) cin >> y[i];
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    cout << "i\t x\t  y\t  xy\t   x^2\n";
    for (int i = 0; i < n; i++) {
        double xy = x[i] * y[i];
        double x2 = x[i] * x[i];
        cout << i+1 << "\t "
             << x[i] << "\t  "
             << y[i] << "\t  "
             << xy << "\t   "
             << x2 << endl;
        sumX += x[i];
        sumY += y[i];
        sumXY += xy;
        sumX2 += x2;
    }
    cout << "\nsum x   = " << sumX;
    cout << "\nsum y   = " << sumY;
    cout << "\nsum xy  = " << sumXY;
    cout << "\nsum x(square) = " << sumX2 << endl;

    double b = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    double a = (sumY - b * sumX) / n;
    cout << "\nRegression Line: y = " << a << " + " << b << "x\n";
    return 0;
}
```
[Open linear. cpp](./src/Curve%20Fitting/regression/Linear/linear.cpp)

#### <a name="linear-regression-input"></a>Input
```
[Add input format/example here]
4
1 2 3 4
3 5 7 9
```
[Open input_lin.txt](./src/Curve%20Fitting/regression/Linear/input_lin.txt)

#### <a name="linear-regression-output"></a>Output
```
[Add output format/example here]
i	 x	  y	  xy	   x^2
1	 1	  3	  3	   1
2	 2	  5	  10	   4
3	 3	  7	  21	   9
4	 4	  9	  36	   16

sum x   = 10
sum y   = 24
sum xy  = 70
sum x(square) = 30

Regression Line: y = 1 + 2x
```
[Open output_lin.txt](./src/Curve%20Fitting/regression/Linear/output_lin.txt)

---

### Polynomial Regression

#### <a name="polynomial-regression-theory"></a>Theory
[Add theory content here]

#### <a name="polynomial-regression-code"></a>Code
```cpp
// View the code file here:
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

vector<double> gaussianElimination(vector<vector<double>> a, int n) {
    for (int i = 0; i < n; i++) {
        
        for (int k = i + 1; k < n; k++)
            if (abs(a[k][i]) > abs(a[i][i]))
                swap(a[i], a[k]);

        for (int k = i + 1; k < n; k++) {
            double t = a[k][i] / a[i][i];
            for (int j = i; j <= n; j++)
                a[k][j] -= t * a[i][j];
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = a[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= a[i][j] * x[j];
        x[i] /= a[i][i];
    }
    return x;
}

int main() {
    freopen("input_pol.txt" , "r" , stdin);
    freopen("output_pol.txt" , "w" ,stdout);
    int n, m;
    //cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n), y(n);
   // cout << "Enter x values:\n";
    for (int i = 0; i < n; i++) cin >> x[i];
    //cout << "Enter y values:\n";
    for (int i = 0; i < n; i++) cin >> y[i];

    //cout << "Enter degree of polynomial (m): ";
    cin >> m;

   
    vector<double> X(2 * m + 1, 0);
    for (int i = 0; i < 2 * m + 1; i++)
        for (int j = 0; j < n; j++)
            X[i] += pow(x[j], i);

    vector<vector<double>> B(m + 1, vector<double>(m + 2, 0));
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= m; j++)
            B[i][j] = X[i + j];

        for (int j = 0; j < n; j++)
            B[i][m + 1] += pow(x[j], i) * y[j];
    }

    vector<double> a = gaussianElimination(B, m + 1);

    
    vector<double> sumCols(m + 3, 0); 
    
    cout << "\n i\t x\t y\t";
    for (int j = 1; j <= m; j++)
        cout << "x^" << j << "\t";
    cout << "y_pred\n";

    for (int i = 0; i < n; i++) {
        cout << i + 1 << "\t" << x[i] << "\t" << y[i] << "\t";

        sumCols[0] += x[i];
        sumCols[1] += y[i];

        double y_pred = a[0];
        for (int j = 1; j <= m; j++) {
            double x_pow = pow(x[i], j);
            cout << x_pow << "\t";
            y_pred += a[j] * x_pow;
            sumCols[j + 1] += x_pow;
        }
        cout << y_pred << endl;
        sumCols[m + 2] += y_pred;
    }

    cout << "sum\t" << sumCols[0] << "\t" << sumCols[1] << "\t";
    for (int j = 2; j <= m + 1; j++)
        cout << sumCols[j] << "\t";
    cout << sumCols[m + 2] << endl;

    cout << "\nPolynomial Regression Equation:\ny = ";
    for (int i = 0; i <= m; i++) {
        cout << fixed << setprecision(4) << a[i];
        if (i > 0) cout << "x^" << i;
        if (i != m) cout << " + ";
    }
    cout << endl;

    return 0;
}
```
[Open polynomial.cpp](./src/Curve%20Fitting/regression/polynomial/polynomial.cpp)

#### <a name="polynomial-regression-input"></a>Input
```
[Add input format/example here]
3
1 2 3
2 3 4
1
```
[Open input_pol.txt](./src/Curve%20Fitting/regression/polynomial/input_pol.txt)

#### <a name="polynomial-regression-output"></a>Output
```
[Add output format/example here]

 i	 x	 y	x^1	y_pred
1	1	2	1	2
2	2	3	2	3
3	3	4	3	4
sum	6	9	6	9

Polynomial Regression Equation:
y = 1.0000 + 1.0000x^1
```
[Open output_pol.txt](./src/Curve%20Fitting/regression/polynomial/output_pol.txt)

---

### Transcendental Regression

#### <a name="transcendental-regression-theory"></a>Theory
[Add theory content here]

#### <a name="transcendental-regression-code"></a>Code
```cpp
// View the code file here:
#include<bits/stdc++.h>
using namespace std;
 
pair<double , double> regression(vector<double>vx , vector<double>vy)
{
    double n = vx.size();
    double sx=0 , sy=0 , sxy=0 ,sx2=0;
    for(int i=0;i<n;i++)
    {
        double nx= exp(vx[i]/4);
        sx+= nx;
        sy+= vy[i];
        sxy+= nx*vy[i];
        sx2+= nx*nx;
    }
    double q = (n* sxy - sx*sy)/ (n *sx2 - sx*sx);
    double p = (sy - q*sx) /n;
    return {p , q};
 
}
 
int main(){
    freopen("input_trans.txt", "r", stdin);
    freopen("output_trans.txt", "w", stdout);
    int n;
    //cout << "Enter number of data points: ";
    cin >> n;
    
    vector<double>vz(n) , vy(n);
    //cout << "Enter z and y values (z y pairs):" << endl;
    for(int i=0;i<n;i++)
    {
        cin >> vz[i] >> vy[i];
    }

    auto [p , q] = regression(vz , vy);
    
    cout << "\nCalculated Parameters:" << endl;
    cout << "p = " << p << "  q = " << q << endl;
    cout << "Equation: y = " << p << " + " << q << "e^(z/4)" << endl;

    double z;
    cout << "\nEnter z to predict: ";
    cin >> z;
    
    double y = p + q * exp(z/4);
    cout << "Estimated y: " << y << endl;
 
}
 
```
[Open regression_transcendental.cpp](./src/Curve%20Fitting/regression/trascendental/regression_transcendental.cpp)

#### <a name="transcendental-regression-input"></a>Input
```
[Add input format/example here]
3
0 95
4 86.4
10 39.08
2
```
[Open input_trans.txt](./src/Curve%20Fitting/regression/trascendental/input_trans.txt)

#### <a name="transcendental-regression-output"></a>Output
```
[Add output format/example here]

Calculated Parameters:
p = 99.9968  q = -5.00041
Equation: y = 99.9968 + -5.00041e^(z/4)

Enter z to predict: Estimated y: 91.7525
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
#include <bits/stdc++.h>
using namespace std;
int degree;
vector<double> coeffs;
double f(double x)
{
    double result = coeffs[0];
    for (int i = 1; i <= degree; i++) {
        result = result * x + coeffs[i];
    }
    return result;
}
int main() {
    freopen("input1_3.txt", "r", stdin);
    freopen("output1_3.txt", "w", stdout);
    cin >> degree;
    coeffs.resize(degree + 1);
    for (int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }
    double a, b;
    cin >> a;
    cin >> b;
    int n;
    cin >> n;
    if (n % 2 != 0) {
        cout << "\nError: Number of intervals (n) must be even for Simpson's 1/3 Rule." << endl;
        cout << "You entered n = " << n << endl;
        return 1;
    }
    double h = (b - a) / n;
    double sumev = 0, sumodd = 0;
    vector<double> x(n+1), y(n+1);
    x[0] = a;
    y[0] = f(a);
    for(int i = 1; i <= n; i++)
    {
        x[i] = a + i * h;
        y[i] = f(x[i]);
    }
    for(int i = 1; i < n; i++)
    {
        if(i % 2 == 0)
            sumev += y[i];
        else
            sumodd += y[i];
    }
    double finalval = (h / 3.0) * ((y[0] + y[n]) + 4 * sumodd + 2 * sumev);
    cout << "Simpson's 1/3 Rule Integration" << endl;
    cout << "Integration value: " << fixed << setprecision(6) << finalval << endl;
    return 0;
}  
```
[Open one_third.cpp](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/one%20_third. cpp)

#### <a name="simpsons-13-input"></a>Input
```
[Add input format/example here]
2
1 0 0
0 2
4
```
[Open input1_3.txt](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/input1_3.txt)

#### <a name="simpsons-13-output"></a>Output
```
[Add output format/example here]
Simpson's 1/3 Rule Integration
Integration value: 2.666667
```
[Open output1_3.txt](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/output1_3.txt)

---

### Simpson's 3/8 Rule

#### <a name="simpsons-38-theory"></a>Theory
[Add theory content here]

#### <a name="simpsons-38-code"></a>Code
```cpp
// View the code file here:
#include<bits/stdc++.h>
using namespace std;
int degree;
vector<double> coeffs;
double func(double x) {
    double result = coeffs[0];
    for (int i = 1; i <= degree; i++) {
        result = result * x + coeffs[i];
    }
    return result;
}
int main() {
    freopen("input_3_8.txt" , "r" , stdin);
    freopen("output_3_8.txt" , "w" ,stdout);
    cout << "Newton's 3/8 Rule for Numerical Integration" << endl;
    cin >> degree;
    coeffs.resize(degree + 1);
    for (int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }
    double lower, upper;
    int n;
    cin >> lower;
    cin >> upper;
    cin >> n;
    if (n % 3 != 0) {
        cout << "Error: 'n' must be a multiple of 3 for Newton's 3/8 Rule." << endl;
        return 1;
    }
    double h = (upper - lower) / n;
    double sum = 0.0;
    sum = func(lower) + func(upper);
    for (int i = 1; i < n; i++) {
        double x = lower + i * h;
        if (i % 3 == 0) {
            // Multiplier is 2 for indices divisible by 3 (3, 6, 9...)
            sum += 2 * func(x);
        } else {
            // Multiplier is 3 for others (1, 2, 4, 5...)
            sum += 3 * func(x);
        }
    }
    // Apply the final 3h/8 factor
    double result = (3 * h / 8) * sum;
    cout << fixed << setprecision(6);
    cout << "The integral value is: " << result << endl;
    return 0;
} 
```
[Open newton3_8.cpp](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/newton3_8.cpp)

#### <a name="simpsons-38-input"></a>Input
```
[Add input format/example here]
3
1 0 0 1
0 3
6
```
[Open input_3_8.txt](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/input_3_8.txt)

#### <a name="simpsons-38-output"></a>Output
```
[Add output format/example here]
Newton's 3/8 Rule for Numerical Integration
The integral value is: 23.250000
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
