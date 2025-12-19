# Numerical Lab Project
**CSE 2208 - Numerical Methods Implementation**

---

# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#theory)
    - [Code](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_ELIMINATION/Gauss_Elimination.cpp)
    - [Input](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_ELIMINATION/input.txt)
    - [Output](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_ELIMINATION/output.txt)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#theory-1)
    - [Code](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/Gauss_Jordan.cpp)
    - [Input](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/input. txt)
    - [Output](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/GAUSS_JORDAN/output.txt)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#theory-2)
    - [Code](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/LU_DECOMPOSITION/LU_DECOMPOSITION.cpp)
    - [Input](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/input.txt)
    - [Output](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/LU_DECOMPOSITION/output.txt)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#theory-3)
    - [Code](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/Matrix_Inversion.cpp)
    - [Input](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/input.txt)
    - [Output](./src/SOLUTION%20OF%20LINEAR%20EQUATIONS/MATRIX_INVERSION/output.txt)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#theory-4)
    - [Code](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/Bisection.cpp)
    - [Input](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/input.txt)
    - [Output](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/BISECTION/output.txt)
  - [False Position Method](#false-position-method)
    - [Theory](#theory-5)
    - [Code](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/FALSE_POSITION/False_Position.cpp)
    - [Input](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/FALSE_POSITION/input.txt)
    - [Output](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/FALSE_POSITION/output.txt)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#theory-6)
    - [Code](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/Newton_Raphson.cpp)
    - [Input](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/input_nr.txt)
    - [Output](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/Newto-Raphson/output_nr.txt)
  - [Secant Method](#secant-method)
    - [Theory](#theory-7)
    - [Code](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/Secant. cpp)
    - [Input](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/input_sec.txt)
    - [Output](./src/SOLUTION%20OF%20NON%20LINEAR%20EQUATIONS/SECANT_METHOD/output_sec. txt)

- [Newton Interpolation](#newton-interpolation)
  - [Forward Interpolation](#forward-interpolation)
    - [Theory](#theory-8)
    - [Code](./src/Newton%20interpolation/Forward/newton_forward.cpp)
    - [Input](./src/Newton%20interpolation/Forward/input_fw.txt)
    - [Output](./src/Newton%20interpolation/Forward/output_fw. txt)
  - [Backward Interpolation](#backward-interpolation)
    - [Theory](#theory-9)
    - [Code](./src/Newton%20interpolation/Backward/newton_backward.cpp)
    - [Input](./src/Newton%20interpolation/Backward/input_bw.txt)
    - [Output](./src/Newton%20interpolation/Backward/output_bw.txt)
  - [Divided Difference](#divided-difference)
    - [Theory](#theory-10)
    - [Code](./src/Newton%20interpolation/divided%20difference/div. cpp)
    - [Input](./src/Newton%20interpolation/divided%20difference/input_div.txt)
    - [Output](./src/Newton%20interpolation/divided%20difference/output_div.txt)

- [Curve Fitting](#curve-fitting)
  - [Linear Regression](#linear-regression)
    - [Theory](#theory-11)
    - [Code](./src/Curve%20Fitting/regression/Linear/linear. cpp)
    - [Input](./src/Curve%20Fitting/regression/Linear/input_lin.txt)
    - [Output](./src/Curve%20Fitting/regression/Linear/output_lin.txt)
  - [Polynomial Regression](#polynomial-regression)
    - [Theory](#theory-12)
    - [Code](./src/Curve%20Fitting/regression/polynomial/polynomial.cpp)
    - [Input](./src/Curve%20Fitting/regression/polynomial/input_pol.txt)
    - [Output](./src/Curve%20Fitting/regression/polynomial/output_pol.txt)
  - [Transcendental Regression](#transcendental-regression)
    - [Theory](#theory-13)
    - [Code](./src/Curve%20Fitting/regression/trascendental/regression_transcendental.cpp)
    - [Input](./src/Curve%20Fitting/regression/trascendental/input_trans.txt)
    - [Output](./src/Curve%20Fitting/regression/trascendental/output_trans.txt)

- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3 Rule](#simpsons-13-rule)
    - [Theory](#theory-14)
    - [Code](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/one%20_third. cpp)
    - [Input](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/input1_3.txt)
    - [Output](./src/NUMERICAL_INTEGRATION/ONE_THIRD_RULE/output1_3.txt)
  - [Simpson's 3/8 Rule](#simpsons-38-rule)
    - [Theory](#theory-15)
    - [Code](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/newton3_8.cpp)
    - [Input](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/input_3_8.txt)
    - [Output](./src/NUMERICAL_INTEGRATION/THREE_EIGHT_RULE/output_3_8.txt)

- [Direct Differentiation](#direct-differentiation)
  - [Forward Differentiation](#forward-differentiation)
    - [Theory](#theory-16)
    - [Code](./src/Direct%20Differentiation/forward/Differentiation_FI.cpp)
    - [Input](./src/Direct%20Differentiation/forward/input_diff_fwd.txt)
    - [Output](./src/Direct%20Differentiation/forward/output_diff_fwd.txt)
  - [Backward Differentiation](#backward-differentiation)
    - [Theory](#theory-17)
    - [Code](./src/Direct%20Differentiation/backward/differentiation_BI.cpp)
    - [Input](./src/Direct%20Differentiation/backward/input_diff. txt)
    - [Output](./src/Direct%20Differentiation/backward/output_diff. txt)

- [Solution of Differential Equations](#solution-of-differential-equations)
  - [Euler Method](#euler-method)
    - [Theory](#theory-18)
    - [Code](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/EULER/Euler.cpp)
    - [Input](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/EULER/input.txt)
    - [Output](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/EULER/output. txt)
  - [Runge Kutta Method](#runge-kutta-method)
    - [Theory](#theory-19)
    - [Code](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/RungeKutta4th.cpp)
    - [Input](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/input.txt)
    - [Output](./src/SOLUTION%20OF%20DIFFERENTIAL%20EQUATIONS/Runge%20Kutta/output.txt)

---

## Solution of Linear Equations

### Gauss Elimination Method

#### Theory
[Add theory content here]

---

### Gauss Jordan Elimination Method

#### Theory
[Add theory content here]

---

### LU Decomposition Method

#### Theory
[Add theory content here]

---

### Matrix Inversion

#### Theory
[Add theory content here]

---

## Solution of Non-Linear Equations

### Bisection Method

#### Theory
[Add theory content here]

---

### False Position Method

#### Theory
[Add theory content here]

---

### Newton Raphson Method

#### Theory
[Add theory content here]

---

### Secant Method

#### Theory
[Add theory content here]

---

## Newton Interpolation

### Forward Interpolation

#### Theory
[Add theory content here]

---

### Backward Interpolation

#### Theory
[Add theory content here]

---

### Divided Difference

#### Theory
[Add theory content here]

---

## Curve Fitting

### Linear Regression

#### Theory
[Add theory content here]

---

### Polynomial Regression

#### Theory
[Add theory content here]

---

### Transcendental Regression

#### Theory
[Add theory content here]

---

## Numerical Integration

### Simpson's 1/3 Rule

#### Theory
[Add theory content here]

---

### Simpson's 3/8 Rule

#### Theory
[Add theory content here]

---

## Direct Differentiation

### Forward Differentiation

#### Theory
[Add theory content here]

---

### Backward Differentiation

#### Theory
[Add theory content here]

---

## Solution of Differential Equations

### Euler Method

#### Theory
[Add theory content here]

---

### Runge Kutta Method

#### Theory
[Add theory content here]

---

## Contributing
This is a course project for CSE 2208. Feel free to explore the implementations and learn from the numerical methods. 

## License
Educational purposes only. 
