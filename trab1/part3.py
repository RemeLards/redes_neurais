from part1 import (
    grad_descent,
    exponential_func,
    numerical_grad_op,
) 
from typing import (
    Callable,
    Tuple,
    Optional,
)
from pydantic import (
    BaseModel,
)
import math
import matplotlib.pyplot as plt
import numpy as np

class GradDescMetrics(BaseModel):
    x: list[list[float]] = []
    f: list[float] = []
    df: list[float] = []


def exponential_func_np(x: np.ndarray) -> np.ndarray: 
    return ( np.exp(-x) * x * (x**2-x-1) )

def grad_descent_algorithm(
    func: Callable,
    grad_func: Callable,
    grad_desc_func: Callable,
    args: list[float],
    gmin: float = 0.1,
    k_max: int = 100,
    enable_metrics: bool = True,
    **kwargs: dict,
) -> Tuple[list[float],Optional[GradDescMetrics]]:
    
    metrics = GradDescMetrics() if enable_metrics else None
    next_args = args
    k = 0

    cont=True
    while(cont):
        if enable_metrics:
            metrics.x.append(next_args)
            metrics.f.append(func(next_args))
            metrics.df.append(
                grad_func(
                    func=func,
                    args=next_args,
                    **kwargs,
                )
            )
        next_args = grad_desc_func(
            func=func,
            grad_func=grad_func,
            args=next_args,
            **kwargs,
        )
        grad = grad_func(
            func=func,
            args=next_args,
            **kwargs
        )

        grad_norm = 0
        for dydx in grad:
            grad_norm+= dydx*dydx
        grad_norm = math.sqrt(grad_norm)
        if grad_norm <= gmin:
            cont=False
   

        k += 1
        if k >= k_max:
            cont = False
    
    return next_args,metrics
        

def main():
    result,metrics = grad_descent_algorithm(
        exponential_func,
        numerical_grad_op,
        grad_descent,
        [3],
        gmin=0.001,
        k_max=100,
        alfa=0.1,
        dh=0.01,
    )

    for i, (x, f, df) in enumerate(zip(metrics.x, metrics.f, metrics.df)):
        # Sei que é unidimensional
        print(f"Iteração {i:4d} | x={x[0]:+8.6f} | f(x)={f:+8.6f} | f'(x)={df[0]:+8.6f}")
    
    print(f"Algorítmo de Descida de Gradiente : {result}")
    print(f"Valor da função:                    {exponential_func(result)}")
    
    fig, ax = plt.subplots(3, 1, figsize=(8, 6))

    x = [x for x in range(len(metrics.x))]

    ax[0].scatter(x,metrics.f)
    ax[0].set_title("f(x) x Iteração")
    ax[0].grid()
    ax[0].set_xlabel("Iteração")
    ax[0].set_ylabel("f(x)")


    # gráfico da derivada
    ax[1].scatter(x,metrics.df,color="red")
    ax[1].set_title("f'(x) x Iteração")
    ax[1].grid()
    ax[1].set_xlabel("Iteração")
    ax[1].set_ylabel("f(x)'")

    x = np.linspace(-1, 5, 100)
    y = exponential_func_np(x)
    ax[2].set_title("Trajetória da descida de gradiente")
    ax[2].set_ylabel("f(x)")
    ax[2].set_xlabel("x")
    ax[2].plot(x, y, 'b', label='f(x)')
    ax[2].grid()
    step = 5
    xs = metrics.x[::step]
    ys = metrics.f[::step]
    for i in range(len(xs)-1):
        ax[2].annotate(
            "",
            xy=(xs[i+1][0], ys[i+1]),
            xytext=(xs[i][0], ys[i]),
            arrowprops=dict(arrowstyle="->", lw=1, mutation_scale=20),
        )
    ax[2].scatter(
        metrics.x[::step],
        metrics.f[::step],
        s=50,
        color="r"
    )

    plt.tight_layout()
    plt.savefig('trab1_part3.png', dpi=300, bbox_inches='tight')
    plt.show()



if __name__ == "__main__":
    main()
