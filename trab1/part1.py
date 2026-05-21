import math
from typing import Callable

alfa = 0.1

# Desisti de deixar genérico para a sáida de uma função
# Ao invés de ser apenas um valor ser um vetor de valores
# Estava dando trabalho e creio que para esse trabalho foge do escopo
def exponential_func(x: list[float]) -> float: 
    return ( math.exp(-x[0]) * x[0] * (x[0]**2-x[0]-1) )

def numerical_grad_op(
    func: Callable,
    args: list[float],
    **kwargs: dict,
) -> list[float]:
    
    dh = kwargs.get("dh") if kwargs.get("dh") else 0.01
    grad = []

    for i in range(len(args)):

        args_plus = args.copy()
        args_minus = args.copy()

        args_plus[i] += dh
        args_minus[i] -= dh

        f1 = func(args_plus)
        f2 = func(args_minus)

        grad.append((f1-f2)/(2*dh))
        
        #print(f"Gradiente : {grad}")
        
    return grad

def grad_descent(
    func: Callable,
    grad_func: Callable, 
    args: list[float],
    **kwargs: dict,
) -> list[float]:
    
    alfa = kwargs.get("alfa") if kwargs.get("alfa") else 0.1

    grad = grad_func(
        func,
        args,
        **kwargs,
    )
    next_x = []
    for dy,x in zip(grad,args):
        next_x.append(x-alfa*dy)
    return next_x


def main():
    print(f"Próximo vetor de entrada {grad_descent(exponential_func,numerical_grad_op,[1])}")


if __name__ == "__main__":
    main()
