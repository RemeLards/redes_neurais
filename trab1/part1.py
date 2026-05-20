import math
from typing import Callable

alfa = 0.1

def exponential_func(x: list[float]) -> list[float]:
    return [( math.exp(-x[0]) * x[0] * (x[0]**2-x[0]-1) )]

def numerical_grad_op(
    func: Callable,
    args: list[float],
    dh: float = 0.01
) -> list[float]:

    args1 = [x+dh for x in args]
    args2 = [x-dh for x in args]

    f1 = func(args1)
    f2 = func(args2)

    grad = []
    for y1,y2 in zip(f1,f2):
        grad.append((y1-y2)/(2*dh))
    
    #print(f"Gradiente : {grad}")
    
    return grad

def grad_descent(
    func: Callable, 
    args: list[float],
    dh: float = 0.01,
    alfa: float = 0.1,
) -> list[float]:
    grad = numerical_grad_op(func,args,dh)
    next_x = []
    for dy,x in zip(grad,args):
        next_x.append(x-alfa*dy)
    return next_x


def main():
    print(f"Próximo vetor de entrada {grad_descent(exponential_func,[1])}")


if __name__ == "__main__":
    main()
