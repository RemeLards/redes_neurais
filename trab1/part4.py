from part1 import (
    grad_descent,
    alfa,
)
from sympy import (
    diff,
    lambdify,
    symbols,
)
import numpy as np

def bi_dim_func(x: list[float]) -> float:
    return ( (4 - 2.1 * x[0]**2 + x[0]**2 / 3 ) * x[0]**3 + x[0] * x[1] + (-4 + 4 * x[1]**2) *x[1]**2 )


def main():
    # a)
    print(f"Próximo vetor de entrada         {grad_descent(bi_dim_func,[1,1])}")

    # b)
    x1,x2 = symbols('x1,x2')
    expr = ((4 - 2.1 * x1**2 + x1**2 / 3 ) * x1**3 + x1 * x2 + (-4 + 4 * x2**2) * x2**2)
    # print(expr)
    df_x1 = diff(expr, x1)
    df_x2 = diff(expr, x2)
    # print(df)
    f1 = lambdify([x1,x2], df_x1, 'numpy')
    f2 = lambdify([x1,x2], df_x2, 'numpy')
    a = np.array([1.0,1.0])
    # print(f(a))
    grad = np.array([
        f1(a[0], a[1]),
        f2(a[0], a[1])
    ])
    a -= alfa * grad
    print(f"Próximo vetor de entrada (Sympy) {a.tolist()}")



if __name__ == "__main__":
    main()

