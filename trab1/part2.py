from part1 import (
    grad_descent,
    exponential_func,
    alfa,
) 
from sympy import (
    diff,
    lambdify,
    symbols,
    exp,
)
import numpy as np

alfa = 0.01

def main():
    print(f"Próximo vetor de entrada         {grad_descent(exponential_func,[1])}")

    x_symbol = symbols('x')
    expr = exp(-x_symbol) * x_symbol*(x_symbol**2-x_symbol-1)
    # print(expr)
    df = diff(expr, x_symbol)
    # print(df)
    f = lambdify(x_symbol, df, 'numpy')
    a = np.array([1.0])
    # print(f(a))
    a -= alfa * f(a)
    print(f"Próximo vetor de entrada (Sympy) {a}")


if __name__ == "__main__":
    main()
