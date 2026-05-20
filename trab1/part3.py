from part1 import (
    grad_descent,
    exponential_func,
    numerical_grad_op,
) 
from typing import (
    Callable,
    Union,
)
from pydantic import (
    BaseModel,
)

class GradDescMetrics(BaseModel):
    x: list[list[float]] = []
    f: list[float] = []
    df: list[float] = []

def grad_descent_algorithm(
    func: Callable,
    args: list[float],
    dh: float = 0.01,
    alfa: float = 0.1,
    gmin: float = 0.1,
    k_max: int = 100,
    enable_metrics: bool = True,
) -> Union[list[float],GradDescMetrics]:
    
    metrics = GradDescMetrics()
    next_args = args
    k = 0

    cont=True
    while(cont):
        if enable_metrics:
            metrics.x.append(next_args)
            metrics.f.append(func(next_args))
            metrics.df.append(
                numerical_grad_op(func=exponential_func,
                    args=next_args,
                    dh=dh,
                )
            )
        next_args = grad_descent(
            func=exponential_func,
            args=next_args,
            alfa=alfa,
            dh=dh,
        )
        grad = numerical_grad_op(
            func=exponential_func,
            args=next_args,
            dh=dh,
        )

        for dy in grad:
            if dy <= gmin:
                cont = False

        k += 1
        if k >= k_max:
            cont = False
    
    return next_args,metrics
        

def main():
    result,metrics = grad_descent_algorithm(
        exponential_func,
        [3],
        alfa=0.1,
        gmin=0.1,
        k_max=100
    )
    print(f"Algorítmo de Descida de Gradiente : {result}")
    for i, (x, f, df) in enumerate(zip(metrics.x, metrics.f, metrics.df)):
        # Sei que é unidimensional
        print(f"Iteração {i:4d} | x={x[0]:+8.6f} | f(x)={f:+8.6f} | f'(x)={df[0]:+8.6f}")



if __name__ == "__main__":
    main()
