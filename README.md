# Numerical Lab Project
**CSE 2208 - Numerical Methods Implementation**
##  Contributors

| Name | Roll Number |
|------|-------------|
| MD KHORSHED SHEIKH | 2207013 |
| ALIF AL AHAD| 2207016 |
| MD MEHEDI HASAN SHUVO | 2207028 

---

# Table of Contents
<img width="1347" height="782" alt="image" src="https://github.com/user-attachments/assets/45ad1e0e-ae3a-4779-87f8-93d1da7aa9de" />


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
<details open>
<summary><b>📖 Theory</b></summary>

<br>

> ### Introduction
> 
> The Gauss Elimination method is a direct numerical technique used to solve systems of linear equations of the form **Ax = b**. 

<br>

### 📘 Theory

The method is theoretically based on the concept of **row equivalence**. It relies on the principle that applying **elementary row operations** (such as swapping rows, multiplying a row by a non-zero constant, or adding a multiple of one row to another) preserves the solution set of the system. 

By systematically applying these operations, the method transforms the dense coefficient matrix **A** into an **upper triangular matrix** (row-echelon form). Once in this simplified form, the unknowns can be solved sequentially without complex matrix inversion. 

<br>

### ⚙️ Algorithm

The process consists of two main stages: **Forward Elimination** and **Back Substitution**. 

<br>

<table>
<tr>
<td width="50%" valign="top">

#### 🔽 Stage 1: Forward Elimination

**Step 1:** Form the Augmented Matrix  
Combine the coefficient matrix and the constant vector into a single matrix **[A | b]**.

**Step 2:** Pivot Selection  
Start with the first row. The first element is the **"pivot"**.

**Step 3:** Eliminate Entries Below Pivot  
Use row operations to create zeros in the column below the current pivot. 

```
Formula: Row_i = Row_i – (factor) × Pivot_Row
Factor: (Element to eliminate) / (Pivot element)
```

**Step 4:** Iterate  
Move to the next diagonal element (next pivot) and repeat the process until the matrix is in **upper triangular form** (all elements below the main diagonal are zero).

</td>
<td width="50%" valign="top">

#### 🔼 Stage 2: Back Substitution

**Step 1:** Solve for the Last Variable  
The last row now contains an equation with only one variable. Solve for it directly.

**Step 2:** Substitute Upwards  
Substitute the known value into the row above it to solve for the next variable. 

**Step 3:** Repeat  
Continue this process moving upwards until all variables **(x₁, x₂, ..., xₙ)** are found.

</td>
</tr>
</table>

<br>

### 💪 Advantages

| Aspect | Description |
|--------|-------------|
| ✅ **Generality** | It can solve any system of n linear equations with n unknowns, provided a unique solution exists. |
| ✅ **Systematic Approach** | The method is algorithmic and easy to program for computers. |
| ✅ **Exactness** | Theoretically, it produces the exact solution (ignoring computer round-off errors) unlike iterative methods that produce approximations. |

<br>

### ⚠️ Disadvantages

| Issue | Description |
|-------|-------------|
| ❌ **Computational Cost** | It is computationally expensive for very large systems, with a time complexity of approximately **O(n³)**. |
| ❌ **Round-off Errors** | In computer implementation, the repeated arithmetic operations can accumulate round-off errors, leading to inaccurate results for ill-conditioned matrices. |
| ❌ **Division by Zero** | If a pivot element is zero, the method fails unless row swapping (partial pivoting) is implemented. |

</details>


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
<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> The Gauss-Jordan Elimination method is a variation of the Gauss Elimination method used to solve systems of linear equations or to find the inverse of a matrix. 

<br>

### Theory

While **Gauss Elimination** reduces the matrix to an **upper triangular form** (requiring back substitution), **Gauss-Jordan Elimination** continues the reduction process until the matrix reaches **Reduced Row Echelon Form (RREF)**. Theoretically, this transforms the coefficient matrix **A** into an **Identity Matrix (I)**. Since the operations are applied to the augmented matrix, the constant vector **b** is simultaneously transformed directly into the solution vector **x**, eliminating the need for back substitution entirely. 

<br>

### Algorithm

<table>
<tr>
<td width="5%"><b>1</b></td>
<td width="95%">

**Form the Augmented Matrix**

Create the matrix **[A | b]** containing coefficients and constants. 

</td>
</tr>

<tr>
<td><b>2</b></td>
<td>

**Pivot Selection & Normalization**

For each row **i** (from 1 to n):
- Select the diagonal element **A[i][i]** as the pivot. 
- **Normalize:** Divide the entire row **i** by the pivot element so that the pivot becomes **1**.

</td>
</tr>

<tr>
<td><b>3</b></td>
<td>

**Elimination (Above and Below)**

For every other row **j** (where **j ≠ i**):
- Eliminate the element in the pivot's column **(A[j][i])** by subtracting a multiple of the pivot row from row **j**. 
- **Formula:** `Row_j = Row_j – (A[j][i] × Row_i)`

</td>
</tr>

<tr>
<td><b>4</b></td>
<td>

**Repeat**

Continue this process for all diagonal elements. Upon completion, the matrix will look like **[I | x]**, where the right-hand column contains the solution values. 

</td>
</tr>
</table>

<br>

### Advantages

| Feature | Description |
|---------|-------------|
| **No Back Substitution** | The variables are solved for directly at the end of the process, simplifying the final extraction of results. |
| **Matrix Inversion** | It is the standard direct method for computing the inverse of a matrix (by augmenting A with the Identity Matrix I instead of vector b). |
| **Simplicity in Code** | The algorithm structure is slightly more uniform than Gauss Elimination since both upper and lower elements are treated similarly. |

<br>

### Disadvantages

| Issue | Description |
|-------|-------------|
| **Higher Computational Cost** | It requires approximately 50% more arithmetic operations than standard Gauss Elimination, making it less efficient for solving simple systems of equations. |
| **Round-off Errors** | Like Gauss Elimination, it is sensitive to error accumulation and requires pivoting strategies to handle small or zero diagonal elements. |
| **Memory Usage** | It modifies the entire matrix structure, which can be inefficient for sparse matrices.  |

</details>

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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> LU Decomposition (also known as LU Factorization) is a numerical method used to solve systems of linear equations, calculate determinants, and invert matrices.   

<br>

### Theory

The fundamental principle of this method is to factorize a square matrix **A** into the product of two matrices: 

- **L (Lower Triangular Matrix):** A matrix where all elements above the main diagonal are zero. 
- **U (Upper Triangular Matrix):** A matrix where all elements below the main diagonal are zero.

**Mathematically:** `A = L × U`

Given a system **Ax = b**, we can substitute **A** to get **(LU)x = b**. By introducing an intermediate vector **y**, the problem is split into two simpler steps: 

- **Ly = b** (solved by forward substitution)
- **Ux = y** (solved by back substitution)

<br>

### Algorithm

The solution process involves three distinct stages: 

<br>

**Stage 1: Decomposition (Factorization)**

Decompose matrix **A** into **L** and **U**. A common approach is **Doolittle's Method**, where the diagonal elements of **L** are set to 1.

- Compute elements of **U** (rows).
- Compute elements of **L** (columns) using the remaining values. 

<br>

**Stage 2: Forward Substitution (Solving for y)**

Solve the lower triangular system **Ly = b** to find the intermediate vector **y**. 

- Since **L** is lower triangular, calculate **y₁** first, then substitute to find **y₂**, and so on. 

<br>

**Stage 3: Back Substitution (Solving for x)**

Solve the upper triangular system **Ux = y** to find the solution vector **x**.

- Since **U** is upper triangular, calculate **xₙ** (the last variable) first, then substitute upwards. 

<br>

### Advantages

- **Efficiency for Multiple Inputs:** Once A is decomposed into L and U, you can solve Ax = b for many different constant vectors b without re-doing the computationally expensive decomposition step. 

- **Determinant Calculation:** The determinant of A is simply the product of the diagonal elements of L and U (often just the product of U's diagonal if L has 1s).

- **Memory Management:** The L and U matrices can often be stored in the same memory space as the original matrix A (in-place decomposition), saving memory.

<br>

### Disadvantages

- **Pivoting Required:** If a pivot element (diagonal element) becomes zero or very small, the decomposition fails or becomes unstable. Partial pivoting (PA = LU) is often required to handle this.

- **Square Matrices Only:** The standard method is strictly defined for square matrices. 

- **Complexity:** For a single system Ax = b, it requires approximately the same number of arithmetic operations (O(n³/3)) as Gauss Elimination, so it offers no speed advantage if the system is solved only once.

</details>

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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> The Matrix Inversion method using determinants and cofactors is a classic analytical approach to solving systems of linear equations **(Ax = b)** or finding the inverse of a matrix.  

<br>

### Theory

For a square matrix **A**, the inverse **A⁻¹** exists if and only if the determinant of **A** (denoted as **|A|**) is non-zero. The inverse is calculated using the formula: 

```
A⁻¹ = (1 / |A|) × adj(A)
```

Where **adj(A)** is the **Adjoint matrix**, which is the transpose of the **Cofactor matrix**. Once the inverse is found, the system solution is **x = A⁻¹b**.

<br>

### Algorithm

The process involves four major steps:  

<br>

#### **Step 1: Calculate the Determinant (|A|)**

Compute the determinant of the matrix.   

- If **|A| = 0**, the matrix is singular (no inverse exists), and the method terminates. 

<br>

#### **Step 2: Find the Cofactor Matrix**

For every element **A[i][j]** in the matrix, calculate its cofactor **C[i][j]**.

- **Minor (M):** The determinant of the sub-matrix left after removing row **i** and column **j**. 
- **Cofactor formula:** `C[i][j] = (-1)^(i+j) × M[i][j]`

<br>

#### **Step 3: Find the Adjoint Matrix (adj(A))**

Take the transpose of the Cofactor matrix.  

- Swap rows with columns: **adj[i][j] = C[j][i]**

<br>

#### **Step 4: Compute Inverse and Solve**

- Multiply every element of the Adjoint matrix by **(1 / |A|)** to get **A⁻¹**.
- Finally, multiply **A⁻¹** by the constant vector **b** to find the solution vector **x**. 

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Theoretical Clarity** | This method provides an exact, closed-form algebraic formula for the solution, which is useful for theoretical proofs and symbolic derivations. |
| **Small Systems** | It is straightforward and easy to implement by hand for small matrices (e.g., 2×2 or 3×3). |
| **Element Independence** | Unlike elimination methods, each element of the inverse can be calculated independently if needed. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Computational Complexity** | This is one of the most computationally expensive methods for large matrices.  Calculating determinants and cofactors for an n × n matrix has a factorial time complexity (O(n!)) if done recursively, making it impractical for n > 4. |
| **Memory Intensive** | It requires storing the cofactor matrix and the adjoint matrix. |
| **Rounding Errors** | For floating-point numbers, the calculation of determinants can become numerically unstable, leading to significant errors in the final result. |

</details>

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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> The Bisection Method is a reliable numerical technique for finding roots of continuous functions. It operates on the principle of the Intermediate Value Theorem:  if a continuous function f(x) changes sign over an interval [a, b], then there exists at least one root within that interval.  

<br>

### Theory

For a root to exist in interval **[a, b]**: 

```
f(a) × f(b) < 0  (opposite signs)
```

<br>

### Algorithm

The Bisection Method follows this iterative process:  

<br>

#### **Step 1: Initialization**

Select interval **[a, b]** where **f(a) × f(b) < 0**.

<br>

#### **Step 2: Midpoint Calculation**

Compute the midpoint: 

```
c = (a + b) / 2
```

<br>

#### **Step 3: Function Evaluation**

Calculate **f(c)**.

<br>

#### **Step 4: Interval Update**

- If **f(a) × f(c) < 0** → New interval: **[a, c]** (set **b = c**)
- If **f(c) × f(b) < 0** → New interval:  **[c, b]** (set **a = c**)

<br>

#### **Step 5: Convergence Check**

Repeat until **|b - a| < tolerance**. 

<br>

#### **Step 6: Result**

Final midpoint **c** approximates the root.

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Guaranteed Convergence** | Ensures convergence for continuous functions with sign change. |
| **Simple Implementation** | Straightforward algorithm with minimal computational requirements. |
| **Reliable and Robust** | Works consistently regardless of function complexity. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Slow Convergence** | Linear convergence rate makes it slower than methods like Newton-Raphson. |
| **Requires Initial Interval** | Must find an interval where the function changes sign. |
| **Inefficient for Multiple Roots** | Not well-suited for handling multiple or complex roots efficiently. |

</details>


#### <a name="bisection-code"></a>Code
```cpp
// View the code file here:
#include<bits/stdc++.h>
using namespace std;

int degree;
vector<double> coeffs;

//polynomial
float find(float x){
    double result = coeffs[0];
    for (int i = 1; i <= degree; i++) {
        result = result * x + coeffs[i];
    }
    return result;
}

float biSection(float x, float y){
    if (find(x) == 0) return x;
    if (find(y) == 0) return y;

    if(find(x) * find(y) > 0) return -1; 

    float p;
    while(abs(y - x) > 0.0001){ 
        p = (x + y) / 2;
        float val = find(p);
        
        if (val == 0) return p;

        if(val * find(x) < 0) 
            y = p;
        else 
            x = p;
    }
    return (x + y) / 2;
}

int main() {

    freopen("input_bs.txt" , "r" , stdin);
    freopen("output_bs.txt" , "w" ,stdout);

    cout << "Bisection Method for Root Finding " << endl;

    cin >> degree;
    if (degree < 1) return 1;

    coeffs.resize(degree + 1);
    for (int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }

    float min_range, max_range;
    cin >> min_range >> max_range;

    set<float> roots;
    float step = 0.1;

    for (float x = min_range; x < max_range; x += step) {
        float x1 = x;
        float x2 = x + step;

        if (find(x1) * find(x2) <= 0) {
            float result = biSection(x1, x2);
            if (result != -1) {
                float rounded = round(result * 10000) / 10000;
                roots.insert(rounded);
            }
        }
    }

    if (roots.empty()) {
        cout << "No roots found." << endl;
    } else {
        cout << "Found " << roots.size() << " root(s):" << endl;
        for (auto r : roots) {
            cout << fixed << setprecision(6) << r << endl;
        }
    }
}  
```
[Open Bisection.cpp](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/Bisection.cpp)

#### <a name="bisection-input"></a>Input
```
[Add input format/example here]
3
1 -4 -1 4
-2 5
```
[Open input.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/input.txt)

#### <a name="bisection-output"></a>Output
```
[Add output format/example here]
Bisection Method for Root Finding 
Found 3 root(s):
-1.000000
1.000000
4.000000

```
[Open output.txt](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/output.txt)

---

### False Position Method

#### <a name="false-position-theory"></a>Theory

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> The False Position Method (Regula Falsi) is an enhanced bracketing technique that combines bisection with linear interpolation.  It typically converges faster than bisection by using a weighted average based on function values.  

<br>

### Theory

**Root approximation formula:**

```
c = [a × f(b) - b × f(a)] / [f(b) - f(a)]
```

Where:
- **[a, b]** is the bracketing interval
- **f(a) · f(b) < 0**

<br>

### Algorithm

The False Position Method follows this iterative process:   

<br>

#### **Step 1: Input**

Degree, coefficients, initial guesses **a** and **b**, tolerance **E**. 

<br>

#### **Step 2: Validation**

Verify **f(a) · f(b) < 0**.

<br>

#### **Step 3: Compute**

Calculate **c** using linear interpolation formula. 

<br>

#### **Step 4: Update Interval**

- If **f(a) · f(c) < 0** → set **b = c**
- Otherwise → set **a = c**

<br>

#### **Step 5: Iterate**

Repeat until **|f(c)| ≤ E**. 

<br>

#### **Step 6: Output**

Root value and iteration count.

<br>

---

### Implementation Features

- Uses `long double` for enhanced precision
- File-based I/O (input. txt, output.txt)
- Displays iteration table with **a**, **b**, **c**, and **f(c)** values
- Polynomial evaluation via coefficient vector

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Faster than Bisection** | Converges more quickly by using linear interpolation instead of simple midpoint. |
| **Guaranteed Convergence** | Ensures convergence with a valid bracketing interval. |
| **Good Accuracy** | Provides reliable results for polynomial equations. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Slow Convergence (Edge Cases)** | Can be slow when one endpoint remains fixed for many iterations. |
| **Requires Sign Change** | Must find an initial interval where the function changes sign. |
| **Not Ideal for Multiple Roots** | Performance degrades when dealing with multiple or closely-spaced roots. |

</details>
-
#### <a name="false-position-code"></a>Code

```cpp
#include <bits/stdc++. h>
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
    freopen("output. txt", "w", stdout);
    
    int test;
    cin >> test;
    
    for (int tt = 1; tt <= test; tt++) {
        cout << "Test Case :  " << tt << "...\n\n";
        
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
        
        cout << setw(10) << "Iteration"
             << setw(15) << "a"
             << setw(15) << "b"
             << setw(15) << "c"
             << setw(20) << "f(c)" << endl;
        
        ld c, fc;
        int iter = 0;
        
        do {
            c = (a * fb - b * fa) / (fb - fa);
            fc = f(c, cof);
            
            cout << setw(10) << ++iter
                 << setw(15) << fixed << setprecision(6) << a
                 << setw(15) << b
                 << setw(15) << c
                 << setw(20) << fc << endl;
            
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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> The Newton-Raphson Method is one of the most powerful and widely-used iterative techniques for solving nonlinear equations f(x) = 0. It leverages the tangent line at an initial guess to rapidly approximate the actual root, offering quadratic convergence under suitable conditions.  

<br>

### Theory

The method uses linear approximation via the derivative:   

```
x_(n+1) = x_n - f(x_n) / f'(x_n)
```

Where:
- **x_n** is the current approximation
- **f(x_n)** is the function value
- **f'(x_n)** is the derivative at x_n

<br>

### Algorithm

The Newton-Raphson Method follows this iterative process:   

<br>

#### **Step 1: Initialize**

Choose **x₀** near the expected root. 

<br>

#### **Step 2: Compute**

Evaluate **f(x₀)** and **f'(x₀)**.

<br>

#### **Step 3: Update**

Calculate the next approximation:  

```
x₁ = x₀ - [f(x₀) / f'(x₀)]
```

<br>

#### **Step 4: Check Convergence**

If **|x₁ - x₀| < ε**, stop. 

<br>

#### **Step 5: Iterate**

Set **x₀ = x₁** and repeat steps 2-4. 

<br>

#### **Step 6: Result**

Final **x_n** is the approximate root.

<br>

---

### Convergence Behavior

- **Quadratic convergence** near the root (very fast)
- Requires good initial guess
- Derivative must be non-zero

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Extremely Fast Convergence** | Quadratic convergence rate makes it one of the fastest methods available. |
| **Efficient for Well-Behaved Functions** | Performs exceptionally well when the function is smooth and the initial guess is close to the root. |
| **Fewer Iterations** | Typically requires significantly fewer iterations than bracketing methods. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Requires Derivative Calculation** | Must compute or approximate f'(x), which can be difficult or computationally expensive. |
| **May Diverge** | Poor initial guess can cause divergence or oscillation away from the root. |
| **Fails at Critical Points** | Method breaks down when f'(x) = 0 (horizontal tangent). |

</details>
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
<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> The Secant method is an iterative numerical technique used to find the roots of a function f(x) = 0. It is an open method that uses a secant line connecting two points on the curve to approximate the root.   

<br>

### Theory

Unlike the **Newton-Raphson method**, the Secant method does not require the calculation of the derivative **f'(x)**. Instead, it approximates the derivative using finite differences based on two previous points **(xₙ₋₁ and xₙ)**. This makes it particularly useful when differentiating the function is difficult or computationally expensive. The method requires two initial guesses but does not require the root to be bracketed.  

<br>

### Algorithm

The Secant method follows this iterative process:   

<br>

#### **Step 1: Initialize**

- Choose two initial guesses, **x₀** and **x₁**.   
- Set a tolerance error **E** and a maximum number of iterations.   

<br>

#### **Step 2: Iteration Step**

For each iteration, calculate the next approximation **xₙ₊₁** using the formula:   

```
x_(n+1) = x_n - f(x_n) * [ (x_n - x_(n-1)) / (f(x_n) - f(x_(n-1))) ]
```

<br>

#### **Step 3: Check Convergence**

- Calculate the absolute error: **|xₙ₊₁ – xₙ|**.   
- If the error is **≤ E**, the root is found.  

<br>

#### **Step 4: Update**

- Set **xₙ₋₁ = xₙ**
- Set **xₙ = xₙ₊₁**
- Repeat Step 2 until convergence or maximum iterations are reached.  

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **No Derivative Needed** | Unlike Newton-Raphson, it does not require analytical differentiation, making it suitable for complex functions.   |
| **Faster than Bisection** | It has a superlinear convergence rate (approx. 1.618), which is significantly faster than the linear convergence of the Bisection method.  |
| **Efficiency** | After the first step, it requires only one function evaluation per iteration (reusing the previous value). |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **No Guaranteed Convergence** | Unlike the Bisection method, the Secant method may diverge if the initial guesses are not close enough to the root.  |
| **Slower than Newton-Raphson** | Its convergence rate is slightly slower than Newton's method (which is quadratic, 2.0). |
| **Numerical Instability** | If f(xₙ) and f(xₙ₋₁) are nearly equal, the denominator approaches zero, potentially causing large errors or division by zero.  |

</details>



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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Newton's Forward Difference Interpolation is a numerical method used to approximate the value of a function at a point within the range of a given dataset.   

<br>

### Theory

This method is specifically designed for datasets where the input values **(x)** are **equally spaced** (i.e., the step size **h** is constant). It utilizes **forward differences** (denoted by **Δ**) and is most accurate when the point to be interpolated lies near the **beginning (top)** of the data table.   

<br>

### Algorithm

The Newton's Forward Difference Interpolation follows this process:   

<br>

#### **Step 1: Construct the Forward Difference Table**

Calculate the first differences **(Δy)**, second differences **(Δ²y)**, and so on, until the differences become constant or zero.   

- **Δy₀ = y₁ – y₀**
- **Δ²y₀ = Δy₁ – Δy₀**

<br>

#### **Step 2: Calculate u**

Determine the normalized variable **u** using the formula:   

```
u = (x – x₀) / h
```

Where:
- **x** is the point to find
- **x₀** is the first value
- **h** is the step size

<br>

#### **Step 3: Apply Formula**

Use the formula to find **y(x)**:

```
y(x) = y₀ + u(Δy₀) + [u(u-1) / 2!] × Δ²y₀ + [u(u-1)(u-2) / 3!] × Δ³y₀ + ...  
```

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Simplicity** | The formula is straightforward and easy to compute for equally spaced data.   |
| **Accuracy at Start** | It provides high accuracy for values close to the starting value (x₀). |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Restriction** | It is strictly valid only for equally spaced intervals.  |
| **Location Specific** | It becomes less accurate as the point moves towards the end of the table (where Backward Interpolation is preferred). |

</details>



#### <a name="forward-interpolation-code"></a>Code
```cpp
// View the code file here:
 #include <bits/stdc++.h>
using namespace std;

long long fact(int n) {
    long long f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

int main() {
     freopen("input_fw.txt" , "r" , stdin);
    freopen("output_fw.txt" , "w" ,stdout);
    int n;
    //cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n), y(n);

    //cout << "Enter x values (equally spaced):\n";
    for (int i = 0; i < n; i++) cin >> x[i];

    //cout << "Enter y values:\n";
    for (int i = 0; i < n; i++) cin >> y[i];

    // Difference table
    vector<vector<double>> diff(n, vector<double>(n));
    for (int i = 0; i < n; i++) diff[i][0] = y[i];

    // forward difference table
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            diff[i][j] = diff[i+1][j-1] - diff[i][j-1];
        }
    }

    double xp;
    //cout << "Enter the value of x to interpolate: ";
    cin >> xp;

    double h = x[1] - x[0];      
    double s = (xp - x[0]) / h; 

    //  Newton forward formula
    double yp = y[0];
    double term = 1;

    for (int j = 1; j < n; j++) {
        term *= (s - (j - 1)); 
        yp += (term * diff[0][j]) / fact(j);
    }

    cout << "Interpolated value at x = " << xp << " is: " << yp << endl;
    return 0;
}

```
[Open newton_forward.cpp](./src/Newton%20interpolation/Forward/newton_forward.cpp)

#### <a name="forward-interpolation-input"></a>Input
```
[Add input format/example here]
5
1 2 3 4 5
1 4 9 16 25
1.5
```
[Open input_fw.txt](./src/Newton%20interpolation/Forward/input_fw.txt)

#### <a name="forward-interpolation-output"></a>Output
```
[Add output format/example here]
Interpolated value at x = 1.5 is: 2.25
```
[Open output_fw.txt](./src/Newton%20interpolation/Forward/output_fw. txt)

---

### Backward Interpolation

#### <a name="backward-interpolation-theory"></a>Theory
[Add theory content here]

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Newton's Backward Difference Interpolation is the counterpart to the Forward method, used for approximating values near the end of a dataset.   

<br>

### Theory

Like the Forward method, this requires **equally spaced intervals**. However, it uses **backward differences** (denoted by **∇**) and references values starting from the last point **(xₙ)** and moving backward. It is best used when the value to be interpolated lies near the **end (bottom)** of the table.   

<br>

### Algorithm

The Newton's Backward Difference Interpolation follows this process:   

<br>

#### **Step 1: Construct the Backward Difference Table**

Calculate differences similar to the forward method but indexed from the end.   

- **∇yₙ = yₙ – yₙ₋₁**

<br>

#### **Step 2: Calculate u**

Determine the normalized variable **u** using the formula:   

```
u = (x – xₙ) / h
```

Where **xₙ** is the last value in the dataset.

<br>

#### **Step 3: Apply Formula**

Use the formula to find **y(x)**:

```
y(x) = yₙ + u(∇yₙ) + [u(u+1) / 2!] × ∇²yₙ + [u(u+1)(u+2) / 3!] × ∇³yₙ + ...  
```

*Note the use of u+1 instead of u-1.*

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **End-of-Data Accuracy** | It is the most accurate method for estimating values near the end of the dataset.   |
| **Complementary** | Works perfectly in conjunction with the Forward method to cover the whole range.   |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Restriction** | Requires equal spacing of x values.    |
| **Not for Middle** | Less accurate for points located in the exact center of the table (Stirling's or Bessel's formulas are better there). |

</details>

</details>

#### <a name="backward-interpolation-code"></a>Code
```cpp
// View the code file here:
#include <bits/stdc++.h>
using namespace std;

// factorial function
long long fact(int n) {
    long long f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

int main() {

     freopen("input_bw.txt" , "r" , stdin);
    freopen("output_bw.txt" , "w" ,stdout);
    int n;
    //cout << "Enter number of data points: ";
    cin >> n;

    vector<double> x(n), y(n);

   // cout << "Enter x values (equally spaced):\n";
    for (int i = 0; i < n; i++) cin >> x[i];

    //cout << "Enter y values:\n";
    for (int i = 0; i < n; i++) cin >> y[i];

    // Creating backward difference table
    vector<vector<double>> diff(n, vector<double>(n));
    for (int i = 0; i < n; i++) diff[i][0] = y[i];

    // Construct backward difference table
    for (int j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            diff[i][j] = diff[i][j-1] - diff[i-1][j-1];
        }
    }

    double xp;
    //cout << "Enter value of x to interpolate: ";
    cin >> xp;

    double h = x[1] - x[0];
    double s = (xp - x[n-1]) / h;  

    // Newton backward formula
    double yp = y[n-1];  
    double term = 1;

    for (int j = 1; j < n; j++) {
        term *= (s + (j - 1)); 
        yp += (term * diff[n-1][j]) / fact(j);
    }

    cout << "Interpolated value at x = " << xp << " is: " << yp << endl;

    return 0;
}
```
[Open newton_backward.cpp](./src/Newton%20interpolation/Backward/newton_backward.cpp)

#### <a name="backward-interpolation-input"></a>Input
```
[Add input format/example here]
5
1 2 3 4 5
1 4 9 16 25
4.5
```
[Open input_bw.txt](./src/Newton%20interpolation/Backward/input_bw.txt)

#### <a name="backward-interpolation-output"></a>Output
```
[Add output format/example here]
Interpolated value at x = 4.5 is: 20.25
```
[Open output_bw. txt](./src/Newton%20interpolation/Backward/output_bw.txt)

---

### Divided Difference

#### <a name="divided-difference-theory"></a>Theory
[Add theory content here]
<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Newton's Divided Difference Interpolation is a versatile general formula used for interpolation when the interval step sizes are unequal.    

<br>

### Theory

Unlike the Forward/Backward methods, this does not assume the input values **(x)** are equally spaced. It uses **divided differences**, which are differences in function values divided by the differences in their corresponding x values. The resulting polynomial is constructed progressively, meaning new data points can be added without recalculating the entire set.    

<br>

### Algorithm

The method follows a systematic approach to construct an interpolating polynomial:    

<br>

#### **Step 1: Calculate Divided Differences**

Compute the first, second, and higher-order divided differences.    

**First Order:**
```
f[x₀, x₁] = (f(x₁) – f(x₀)) / (x₁ – x₀)
```

**Second Order:**
```
f[x₀, x₁, x₂] = (f[x₁, x₂] – f[x₀, x₁]) / (x₂ – x₀)
```

<br>

#### **Step 2: Apply Formula**

Construct the interpolating polynomial: 

```
f(x) = f(x₀) + (x – x₀) f[x₀, x₁] + (x – x₀)(x – x₁) f[x₀, x₁, x₂] + ...
```

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Unequal Intervals** | It works for both equally and unequally spaced data, making it universally applicable. |
| **Flexibility** | Adding a new data point only adds one new term to the polynomial; it does not require recalculating previous coefficients (unlike Lagrange Interpolation). |
| **Polynomial Degree** | The degree of the polynomial is evident immediately from the number of terms.    |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Calculation Complexity** | The manual construction of the divided difference table is more arithmetic-heavy than simple difference tables. |
| **Sensitivity** | Like all polynomial interpolation, it can be sensitive to errors in data points (oscillations). |

</details>

</details>

#### <a name="divided-difference-code"></a>Code
```cpp
// View the code file here:
#include<bits/stdc++.h>
using namespace std;

vector<vector<double>> calculateDifferenceTable(int pointCount, const vector<double>& xData, const vector<double>& yData) {
    vector<vector<double>> differenceTable(pointCount, vector<double>(pointCount, 0.0));

    for (int i = 0; i < pointCount; ++i) {
        differenceTable[i][0] = yData[i];
    }

    for (int j = 1; j < pointCount; ++j) {
        for (int i = 0; i < pointCount - j; ++i) {
            double diffTop = differenceTable[i + 1][j - 1] - differenceTable[i][j - 1];
            double diffBottom = xData[i + j] - xData[i];
            differenceTable[i][j] = diffTop / diffBottom;
        }
    }
    return differenceTable;
}

double getInterpolatedValue(double targetValue, int pointCount, const vector<double>& xData, const vector<vector<double>>& differenceTable) {
    double finalValue = differenceTable[0][0];
    double termProduct = 1.0;

    for (int i = 1; i < pointCount; ++i) {
        termProduct *= (targetValue - xData[i - 1]);
        finalValue += differenceTable[0][i] * termProduct;
    }
    return finalValue;
}

double estimateError(double targetValue, int pointCount, const vector<double>& xData, const vector<vector<double>>& differenceTable) {
    int lastIndex = pointCount - 1;
    double errorCoefficient = differenceTable[0][lastIndex];
    
    double termProduct = 1.0;
    for (int i = 0; i < lastIndex; ++i) {
        termProduct *= (targetValue - xData[i]);
    }
    
    return fabs(errorCoefficient * termProduct);
}

void displayTable(int pointCount, const vector<vector<double>>& differenceTable) {
    cout << "\nDifferenceTable:\n\n";
    for (int i = 0; i < pointCount; ++i) {
        for (int j = 0; j < pointCount; ++j) {
            cout << differenceTable[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void processTestCase() {
    int polynomialDegree;
    cin >> polynomialDegree;
    
    int totalPoints = polynomialDegree + 1;
    vector<double> xValues(totalPoints), yValues(totalPoints);

    for (int i = 0; i < totalPoints; ++i) cin >> xValues[i];
    for (int i = 0; i < totalPoints; ++i) cin >> yValues[i];

    double queryX;
    cin >> queryX;

    vector<vector<double>> differenceTable = calculateDifferenceTable(totalPoints, xValues, yValues);
    double result = getInterpolatedValue(queryX, totalPoints, xValues, differenceTable);
    double error = estimateError(queryX, totalPoints, xValues, differenceTable);

    cout << fixed << setprecision(4);
    cout << "n: " << polynomialDegree << "\n";
    
    cout << "x: "; 
    for (auto val : xValues) cout << val << " ";
    cout << "\ny: "; 
    for (auto val : yValues) cout << val << " ";
    cout << "\nX: " << queryX << endl;

    displayTable(totalPoints, differenceTable);

    cout << "Interpolation: " << result << endl;
    cout << "Truncation error: " << error << endl << endl;
}

int main() {
    freopen("input_div.txt", "r", stdin);
    freopen("output_div.txt", "w", stdout);

    processTestCase();
    return 0;
}
```
[Open div. cpp](./src/Newton%20interpolation/divided%20difference/div.cpp)

#### <a name="divided-difference-input"></a>Input
```
[Add input format/example here]
2
0 1 2
0 1 4
1.5
```
[Open input_div.txt](./src/Newton%20interpolation/divided%20difference/input_div.txt)

#### <a name="divided-difference-output"></a>Output
```
[Add output format/example here]
n: 2
x: 0.0000 1.0000 2.0000 
y: 0.0000 1.0000 4.0000 
X: 1.5000

DifferenceTable:

0.0000 1.0000 1.0000 
1.0000 3.0000 0.0000 
4.0000 0.0000 0.0000 

Interpolation: 2.2500
Truncation error: 0.7500

```
[Open output_div.txt](./src/Newton%20interpolation/divided%20difference/output_div.txt)

---

## Curve Fitting

### Linear Regression

#### <a name="linear-regression-theory"></a>Theory
[Add theory content here]

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Linear Regression is a statistical method used to model the relationship between a dependent variable **(y)** and an independent variable **(x)** by fitting a straight line to the observed data. 

<br>

### Theory

The goal is to find the equation of a straight line **y = a + bx** that best fits the data points. The **method of Least Squares** is used, which minimizes the sum of the squares of the vertical differences (residuals) between the observed data points and the fitted line.

<br>

### Algorithm

The Linear Regression process follows these systematic steps: 

<br>

#### **Step 1: Calculate Sums**

Compute the necessary summations from the data set of **n** points: 

- **Σx** (sum of all x values)
- **Σy** (sum of all y values)
- **Σxy** (sum of products of x and y)
- **Σx²** (sum of squares of x)

<br>

#### **Step 2: Solve Normal Equations**

The constants **a** (intercept) and **b** (slope) are found by solving the system of linear equations:

```
Σy = n·a + b·Σx
Σxy = a·Σx + b·Σx²
```

<br>

#### **Step 3: Compute Coefficients**

Alternatively, calculate **a** and **b** directly using the formulas derived from the system above. 

<br>

#### **Step 4: Formulate Equation**

Substitute **a** and **b** into **y = a + bx** to get the final regression line.

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Simplicity** | It is the simplest form of regression and easy to interpret. |
| **Best Fit** | Mathematically guarantees the minimum possible error (in the squared sense) for linear relationships. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Linearity Assumption** | It assumes a linear relationship; if the data is curved, this method yields poor results. |
| **Sensitive to Outliers** | A single extreme value (outlier) can significantly skew the line. |

</details>


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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Polynomial Regression is an extension of linear regression used to model relationships where the data points follow a curve rather than a straight line.  

<br>

### Theory

It fits a polynomial equation of degree **m** to the data: 

```
y = a₀ + a₁x + a₂x² + ... + aₘxᵐ
```

Like linear regression, it uses the **Least Squares method** to minimize the error. For a 2nd-degree (parabolic) curve, the system involves solving for three unknowns **(a₀, a₁, a₂)**.

<br>

### Algorithm

The Polynomial Regression process is structured as follows: 

<br>

#### **Step 1: Choose Degree**

Decide the degree **m** of the polynomial (e.g., **m=2** for a parabola). 

<br>

#### **Step 2: Calculate Sums**

Compute summations for higher powers of x:  

- **Σx, Σx², Σx³, Σx⁴, ..., Σx²ᵐ**
- **Σy, Σxy, Σx²y, .. .**

<br>

#### **Step 3: Form Matrix System**

Construct the system of linear equations (Normal Equations) in matrix form. For a 2nd degree polynomial:  

```
n·a₀ + a₁Σx + a₂Σx² = Σy
a₀Σx + a₁Σx² + a₂Σx³ = Σxy
a₀Σx² + a₁Σx³ + a₂Σx⁴ = Σx²y
```

<br>

#### **Step 4: Solve**

Use **Gaussian Elimination** or **LU Decomposition** to solve for coefficients **a₀, a₁, a₂, .. .**

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Flexibility** | Can model complex, non-linear relationships and curves. |
| **Better Fit** | Provides a much lower error than linear regression for curved data. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Overfitting** | Using a polynomial of too high a degree can cause the curve to oscillate wildly between points, fitting the "noise" rather than the trend. |
| **Complexity** | Requires solving a larger system of equations. |

</details>


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
<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Transcendental regression involves fitting non-polynomial functions—such as exponential, logarithmic, or power functions—to a dataset.   

<br>

### Theory

These equations are **non-linear in their parameters**.  To solve them using standard regression techniques, the equation is usually **linearized using logarithms**.  

**Common Models:**

- **Exponential Model:** `y = ae^(bx)`
- **Power Model:** `y = ax^b`
- **Growth Model:** `y = a + b ln(x)`

<br>

### Algorithm (Linearization Method)

The transcendental regression process uses linearization to transform non-linear equations:   

<br>

#### **Step 1: Identify Model**

Choose the form that best represents the data (e.g., Exponential: **y = ae^(bx)**).  

<br>

#### **Step 2: Linearize the Equation**

Apply logarithms to transform the non-linear equation into a linear form **(Y = A + BX)**.  

**Example (Exponential):**
```
ln(y) = ln(a) + bx
```

Let: 
- **Y = ln(y)**
- **A = ln(a)**
- **B = b**
- **X = x**

<br>

#### **Step 3: Apply Linear Regression**

Perform standard linear regression on the transformed data **(X, Y)** to find **A** and **B**.  

<br>

#### **Step 4: Transform Back**

Convert **A** and **B** back to the original parameters:   

```
a = e^A
b = B
```

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Real-world Modeling** | Many natural phenomena (population growth, radioactive decay, cooling rates) follow transcendental laws, not polynomial ones. |
| **Simplicity** | Linearization allows solving complex curves using simple linear regression tools. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Data Restrictions** | Logarithms cannot be taken of negative numbers or zero; data must be positive. |
| **Distorted Error** | Minimizing the error of the transformed (log) data is not mathematically identical to minimizing the error of the original data, though it is usually close enough. |

</details>



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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Simpson's 1/3 Rule is a numerical integration method used to approximate the definite integral of a function f(x) within a given range [a, b]. It is generally more accurate than the Trapezoidal Rule for the same number of intervals.   

<br>

### Theory

This method approximates the area under the curve by connecting groups of three points using **parabolas (second-degree polynomials)**. Because it uses pairs of sub-intervals to form these parabolas, the total number of sub-intervals **n** must be an **even number**.  

<br>

### Algorithm

The Simpson's 1/3 Rule follows a structured integration approach:   

<br>

#### **Step 1: Discretize the Domain**

- Choose an even number of intervals, **n**.  
- Calculate the step size: **h = (b – a) / n**.  

<br>

#### **Step 2: Generate Points**

Compute the function values **(y₀, y₁, ..., yₙ)** at each step **x**.  

<br>

#### **Step 3: Apply Formula**

Use the Simpson's 1/3 formula:   

```
Area ≈ (h / 3) × [ (First + Last) + 4 × (Sum of Odds) + 2 × (Sum of Evens) ]
```

**Mathematically:**
```
I = (h / 3) [ (y₀ + yₙ) + 4(y₁ + y₃ + ... + yₙ₋₁) + 2(y₂ + y₄ + ... + yₙ₋₂) ]
```

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Higher Accuracy** | It provides much better accuracy (order O(h⁴)) compared to the Trapezoidal Rule (O(h²)) because it accounts for the curvature of the function.  |
| **Exact for Quadratics** | It gives the exact area if the function being integrated is a polynomial of degree 2 or less. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Constraint** | It strictly requires the number of sub-intervals n to be even. If n is odd, the method cannot be applied directly over the whole range. |

<br>






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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> Simpson's 3/8 Rule is another numerical integration technique, similar to the 1/3 rule but based on a higher-order interpolation.    

<br>

### Theory

This method approximates the function using **cubic interpolation (third-degree polynomials)**. It connects **four points at a time**. Consequently, the total number of sub-intervals **n** must be a **multiple of 3** (e. g., 3, 6, 9, 12).  

<br>

### Algorithm

The Simpson's 3/8 Rule follows this integration approach:   

<br>

#### **Step 1: Discretize the Domain**

- Choose a number of intervals **n** that is a multiple of 3.  
- Calculate the step size: **h = (b – a) / n**.   

<br>

#### **Step 2: Generate Points**

Compute the function values **(y₀, y₁, .. ., yₙ)**.  

<br>

#### **Step 3: Apply Formula**

Use the Simpson's 3/8 formula:   

```
Area ≈ (3h / 8) × [ (First + Last) + 3 × (Sum of Non-Multiples of 3) + 2 × (Sum of Multiples of 3) ]
```

**Mathematically:**
```
I = (3h / 8) [ (y₀ + yₙ) + 3(y₁ + y₂ + y₄ + y₅ + ...) + 2(y₃ + y₆ + y₉ + .. .) ]
```

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **Smooth Functions** | It is generally more accurate than the 1/3 rule for functions that are smoother or have higher-order curvature.    |
| **Exact for Cubics** | It gives the exact area if the function is a polynomial of degree 3 or less.   |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Stricter Constraint** | The requirement that n must be a multiple of 3 makes it less flexible than the 1/3 rule or the Trapezoidal rule.  |
| **Computation** | It requires slightly more arithmetic operations per interval than the 1/3 rule. |

</details>


#### <a name="simpsons-38-code"></a>Code

```cpp
#include<bits/stdc++. h>
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
    freopen("input_3_8.txt", "r", stdin);
    freopen("output_3_8.txt", "w", stdout);
    
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

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> This method is used to compute the derivative (rate of change) of a function given in a tabular form. It is specifically applicable when the value of x at which the derivative is required is near the beginning of the table and the data is equally spaced. 

<br>

### Theory

Instead of differentiating the function directly (which may be unknown), we approximate the function using **Newton's Forward Difference Interpolation formula** and then differentiate that polynomial with respect to **x**.

By applying the **chain rule:**

```
dy/dx = (dy/du) × (du/dx)
```

Since **u = (x – x₀) / h**, we know that **du/dx = 1/h**.

<br>

### Algorithm

The forward differentiation process follows these steps: 

<br>

#### **Step 1: Construct Table**

Create the forward difference table up to the required order.

<br>

#### **Step 2: Calculate u**

Find **u = (x – x₀) / h**.

<br>

#### **Step 3: Apply Formula (General)**

**First Derivative:**

```
dy/dx = (1/h) [ Δy₀ + ((2u – 1)/2) Δ²y₀ + ((3u² – 6u + 2)/6) Δ³y₀ + ... ]
```

<br>

#### **Simplified Formula (at x = x₀)**

If calculating exactly at the starting point **(u = 0)**, the formula simplifies to:

```
dy/dx = (1/h) [ Δy₀ – (1/2)Δ²y₀ + (1/3)Δ³y₀ – (1/4)Δ⁴y₀ + ... ]
```

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **No Analytical Function Needed** | Can find the derivative purely from data points. |
| **Simple Calculation** | Very efficient for finding derivatives at tabular points near the start. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Error Amplification** | Differentiation is an "error-magnifying" process; small errors in the data (noise) can lead to large errors in the derivative. |
| **Restriction** | Strictly requires equally spaced intervals. |

</details>

#### <a name="forward-differentiation-code"></a>Code
```cpp
// View the code file here:
#include<bits/stdc++.h>
using namespace std;

double evaluatePolynomial(int degree, vector<double>& coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i <= degree; i++) {
        result += coeffs[i] * pow(x, degree - i);
    }
    return result;
}

int main() {
    freopen("input_diff_fwd.txt", "r", stdin);
    freopen("output_diff_fwd.txt", "w", stdout);

    int degree, n;
    double lower, upper, p;
    
    // Input
    cin >> degree;
    vector<double> coeffs(degree + 1);
    for(int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }
    cin >> lower >> upper >> n >> p;

    double h = (upper - lower) / n;
    vector<double> x_val(n + 1);
    vector<vector<double>> diff(n + 1, vector<double>(n + 1));

    for (int i = 0; i <= n; i++) {
        x_val[i] = lower + i * h;
        diff[i][0] = evaluatePolynomial(degree, coeffs, x_val[i]);
    }

    //Forward Difference Table
    for (int j = 1; j <= n; j++) {
        for (int i = 0; i <= n - j; i++) {
            diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];
        }
    }

    int p_index = -1;
    double min_dist = 1e10;
    for(int i = 0; i <= n; i++) {
        double dist = abs(x_val[i] - p);
        if(dist < min_dist) {
            min_dist = dist;
            p_index = i;
        }
    }

    double first_derivative = 0;

    // Derivative 
    int max_order = min(degree, n - p_index);
    if (max_order >= 1) {
        double term_sum = 0;
        for(int i = 1; i <= max_order; i++) {
            double term = diff[p_index][i] / i;
            if (i % 2 == 0) term_sum -= term; 
            else term_sum += term;            
        }
        first_derivative = (1.0 / h) * term_sum;
    }

    cout << "Point p: " << p << endl;
    cout << "First Derivative: " << fixed << setprecision(6) << first_derivative << endl;

    return 0;
}
```
[Open Differentiation_FI.cpp](./src/Direct%20Differentiation/forward/Differentiation_FI.cpp)

#### <a name="forward-differentiation-input"></a>Input
```
[Add input format/example here]
4
1 0 0 0 0
0 2 100 0.5
```
[Open input_diff_fwd.txt](./src/Direct%20Differentiation/forward/input_diff_fwd.txt)

#### <a name="forward-differentiation-output"></a>Output
```
[Add output format/example here]
Point p: 0.5
First Derivative: 0.500000
```
[Open output_diff_fwd. txt](./src/Direct%20Differentiation/forward/output_diff_fwd.txt)

---

### Backward Differentiation

#### <a name="backward-differentiation-theory"></a>Theory
[Add theory content here]
<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> This method is the counterpart to the forward method and is used to compute the derivative when the point of interest lies near the end of the table. 

<br>

### Theory

We differentiate **Newton's Backward Difference Interpolation formula**. Here, **u = (x – xₙ) / h**, so **du/dx = 1/h**. 

<br>

### Algorithm

The backward differentiation process follows these steps:  

<br>

#### **Step 1: Construct Table**

Create the backward difference table. 

<br>

#### **Step 2: Calculate u**

Find **u = (x – xₙ) / h**.

<br>

#### **Step 3: Apply Formula (General)**

**First Derivative:**

```
dy/dx = (1/h) [ ∇yₙ + ((2u + 1)/2) ∇²yₙ + ((3u² + 6u + 2)/6) ∇³yₙ + ... ]
```

<br>

#### **Simplified Formula (at x = xₙ)**

If calculating exactly at the end point **(u = 0)**, the formula simplifies to:

```
dy/dx = (1/h) [ ∇yₙ + (1/2)∇²yₙ + (1/3)∇³yₙ + (1/4)∇⁴yₙ + ... ]
```

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **End-Point Accuracy** | Best suited for analyzing rates of change for the most recent data (bottom of the table). |
| **Consistent** | Uses the same backward difference table constructed for interpolation. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Truncation Error** | Higher-order differences are often neglected, which can reduce accuracy if the step size h is not small enough. |
| **Stability** | Like all numerical differentiation, it is sensitive to round-off errors.  |

</details>

#### <a name="backward-differentiation-code"></a>Code
```cpp
// View the code file here:
#include<bits/stdc++.h>
using namespace std;

double evaluatePolynomial(int degree, vector<double>& coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i <= degree; i++) {
        result += coeffs[i] * pow(x, degree - i);
    }
    return result;
}

int main() {
    freopen("input_diff.txt", "r", stdin);
    freopen("output_diff.txt", "w", stdout);

    int degree, n;
    double lower, upper, p;
    
    // Input
    cin >> degree;
    vector<double> coeffs(degree + 1);
    for(int i = 0; i <= degree; i++) {
        cin >> coeffs[i];
    }
    cin >> lower >> upper >> n >> p;

    double h = (upper - lower) / n;
    vector<double> x_val(n + 1);
    vector<vector<double>> diff(n + 1, vector<double>(n + 1));

    for (int i = 0; i <= n; i++) {
        x_val[i] = lower + i * h;
        diff[i][0] = evaluatePolynomial(degree, coeffs, x_val[i]);
    }

    //  backward difference table
    for (int j = 1; j <= n; j++) {
        for (int i = j; i <= n; i++) {
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];
        }
    }

    int p_index = -1;
    double min_dist = 1e10;
    for(int i = 0; i <= n; i++) {
        double dist = abs(x_val[i] - p);
        if(dist < min_dist) {
            min_dist = dist;
            p_index = i;
        }
    }

    double first_derivative = 0;
    
    // Derivative
    int max_order = min(degree, p_index);
    if (max_order >= 1) {
        double term_sum = 0;
        
        for(int i = 1; i <= max_order; i++) {
            double term = diff[p_index][i] / i;
            term_sum += term;
        }
        
        first_derivative = (1.0 / h) * term_sum;
    }

    cout << "Point p: " << p << endl;
    cout << "First Derivative: " << fixed << setprecision(6) << first_derivative << endl;

    return 0;
}
```
[Open differentiation_BI.cpp](./src/Direct%20Differentiation/backward/differentiation_BI.cpp)

#### <a name="backward-differentiation-input"></a>Input
```
[Add input format/example here]
6
1 0 0 0 0 0 0
0 2 20 1.8
```
[Open input_diff. txt](./src/Direct%20Differentiation/backward/input_diff.txt)

#### <a name="backward-differentiation-output"></a>Output
```
[Add output format/example here]
Point p: 1.8
First Derivative: 113.374080
```
[Open output_diff.txt](./src/Direct%20Differentiation/backward/output_diff.txt)

---

## Solution of Differential Equations

### Runge Kutta Method

#### <a name="runge-kutta-theory"></a>Theory
[Add theory content here]

<details open>
<summary><b>Theory</b></summary>

<br>

> ### Introduction
> 
> The Runge-Kutta 4th Order method (often referred to simply as RK4) is the most widely used numerical technique for solving Ordinary Differential Equations (ODEs). It provides a solution to initial value problems of the form **dy/dx = f(x, y)**. 

<br>

### Theory

The RK4 method calculates a weighted average of **four different slopes** within the interval.  This results in significantly higher accuracy without requiring the calculation of higher-order derivatives. 

**The four slopes (k) are:**

- **k₁:** Slope at the beginning of the interval. 
- **k₂:** Slope at the midpoint (using k₁ to step halfway).
- **k₃:** Slope at the midpoint (using k₂ to correct the estimate).
- **k₄:** Slope at the end of the interval. 

<br>

### Algorithm

The RK4 method follows this iterative process:  

<br>

#### **Step 1: Define the Problem**

Given the function **f(x, y)**, step size **h**, and initial condition **(x₀, y₀)**.  

<br>

#### **Step 2: Iterate**

For each step **n**, calculate the four intermediate values:

```
k₁ = h · f(xₙ, yₙ)
k₂ = h · f(xₙ + h/2, yₙ + k₁/2)
k₃ = h · f(xₙ + h/2, yₙ + k₂/2)
k₄ = h · f(xₙ + h, yₙ + k₃)
```

<br>

#### **Step 3: Update y**

Calculate the next value of **y** using the weighted average formula: 

```
yₙ₊₁ = yₙ + (1/6) · (k₁ + 2k₂ + 2k₃ + k₄)
```

<br>

#### **Step 4: Update x**

```
xₙ₊₁ = xₙ + h
```

<br>

#### **Step 5: Repeat**

Continue the process until the desired target **x** is reached.  

<br>

---

### Advantages

| Aspect | Benefit |
|--------|---------|
| **High Accuracy** | It has a global error of order O(h⁴), making it much more accurate than Euler's method or the Modified Euler method. |
| **No Higher Derivatives** | Unlike Taylor Series methods, it does not require analytically calculating second or third derivatives of the function, which can be difficult or impossible for complex equations. |
| **Stability** | It is generally stable for a wide range of functions and step sizes. |

<br>

---

### Disadvantages

| Limitation | Impact |
|------------|--------|
| **Computational Cost** | It requires four function evaluations per step (calculating k₁, k₂, k₃, and k₄), whereas Euler's method requires only one. This can be slower for very complex functions. |
| **Complexity** | The implementation is slightly more complex than basic linear methods. |

</details>



#### <a name="runge-kutta-code"></a>Code
```cpp
// View the code file here:
#include <bits/stdc++.h>
using namespace std;
using ld = long double;

// Polynomial derivative function: dy/dx = f(x)
ld dydx(ld x, const vector<ld> &cof) {
    ld sum = 0.0;
    int n = cof.size();
    for (int i = 0; i < n - 1; i++) { // last term is constant
        sum += cof[i] * (n - 1 - i) * powl(x, n - 2 - i);
    }
    return sum;
}

// Runge-Kutta 4th order method (dy/dx = f(x))
ld rungeKutta(ld x0, ld y0, ld xEnd, ld h, const vector<ld> &cof) {
    int n = (int)((xEnd - x0) / h);
    ld y = y0;

    cout << fixed << setprecision(6);
    cout << "Step\t x\t y\n";

    for (int i = 0; i <= n; i++) {
        cout << i << "\t" << x0 << "\t" << y << "\n";
        
        if (i < n) {  // Only update if not at the last step
            ld k1 = h * dydx(x0, cof);
            ld k2 = h * dydx(x0 + h/2.0, cof);
            ld k3 = h * dydx(x0 + h/2.0, cof);
            ld k4 = h * dydx(x0 + h, cof);

            y = y + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
            x0 += h;
        }
    }

    return y;
}

int main() {
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    int test;
    cin >> test;
    for(int tt = 1; tt <= test; tt++)
    {
        
        cout << "\nTest Case " << tt << "...\n\n";
        int degree;
        //    cout << "Enter degree of polynomial: ";
        cin >> degree;
        
        vector<ld> cof(degree + 1);
        //    cout << "Enter coefficients from highest degree to constant:\n";
        for (int i = 0; i <= degree; i++) cin >> cof[i];
        
        ld x0, y0, xEnd, h;
        //    cout << "Enter initial x0, initial y0, final x, and step h:\n";
        cin >> x0 >> y0 >> xEnd >> h;
        
        ld yFinal = rungeKutta(x0, y0, xEnd, h, cof);
        cout << "\nFinal value: y(" << xEnd << ") = " << yFinal << "\n";
    }

    return 0;
}
```
[Open RungeKutta4th.cpp](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/RungeKutta4th.cpp)

#### <a name="runge-kutta-input"></a>Input
```
[Add input format/example here]
2
3
1 0 0 0
0 0 2 0.5
2
4 -3 0
1 1 3 0.4
```
[Open input. txt](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/input.txt)

#### <a name="runge-kutta-output"></a>Output
```
[Add output format/example here]

Test Case 1...

Step	 x	 y
0	0.000000	0.000000
1	0.500000	0.125000
2	1.000000	1.000000
3	1.500000	3.375000
4	2.000000	8.000000

Final value: y(2.000000) = 8.000000

Test Case 2...

Step	 x	 y
0	1.000000	1.000000
1	1.400000	3.640000
2	1.800000	7.560000
3	2.200000	12.760000
4	2.600000	19.240000
5	3.000000	27.000000

Final value: y(3.000000) = 27.000000
```
[Open output.txt](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/output.txt)

---

|

---
## License
Educational purposes only. 
