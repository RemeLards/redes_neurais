from part1 import (
    numerical_grad_op,
    grad_descent,
)
from part3 import(
    grad_descent_algorithm,
    GradDescMetrics
)
from sympy import (
    diff,
    lambdify,
    symbols,
)
from typing import (
    Callable,
    Optional,
)
import numpy as np
import matplotlib.pyplot as plt

def bi_dim_func(x: list[float]) -> float:
    return ( (4 - 2.1 * x[0]**2 + x[0]**2 / 3 ) * x[0]**3 + x[0] * x[1] + (-4 + 4 * x[1]**2) *x[1]**2 )


def bi_dim_func_np(x,y) -> np.ndarray:
    return ( (4 - 2.1 * x**2 + x**2 / 3 ) * x**3 + x * y + (-4 + 4 * y**2) *y**2 )


def sympy_func(args: Optional[list[float]] = None):
    if not args:
        if not hasattr(sympy_func, "df1") or not hasattr(sympy_func, "df2"):
            print("Calculando Gradiente Simbólico...")
            x1,x2 = symbols('x1,x2')
            expr = ((4 - 2.1 * x1**2 + x1**2 / 3 ) * x1**3 + x1 * x2 + (-4 + 4 * x2**2) * x2**2)
            # print(expr)
            df_x1 = diff(expr, x1)
            df_x2 = diff(expr, x2)
            # print(df)
            df1 = lambdify([x1,x2], df_x1, 'numpy')
            df2 = lambdify([x1,x2], df_x2, 'numpy')
            sympy_func.df1 = df1
            sympy_func.df2 = df2
        return sympy_func.df1,sympy_func.df2
    return bi_dim_func(args)


def sympy_grad(
    func: Callable,
    args: list[float],
    **kwargs: dict,
) -> list[float]:
    df1,df2 = func()
    grad = np.array([
        df1(args[0], args[1]),
        df2(args[0], args[1])
    ])

    return grad.tolist()


def sympy_grad_desc(
    func: Callable,
    grad_func: Callable,
    args: list[float],
    **kwargs: dict
) -> list[float]:
    alfa = kwargs.get("alfa") if kwargs.get("alfa") else 0.1
    a = np.array(args)
    grad = np.array(grad_func(func=func,args=args))
    a -= alfa * grad
    # print(f"Próximo vetor de entrada (Sympy) {a.tolist()}")
    return a.tolist()


def plot3d(
    metrics: GradDescMetrics, 
    spacing: int,
    title: str,
    ax=None,
):
    if ax is None:
        ax = plt.axes(projection='3d')

    # Facilita a criação de trajetória (setas) ->
    x1t = np.array([p[0] for p in metrics.x])
    x2t = np.array([p[1] for p in metrics.x])
    ft = np.array(metrics.f)

    dx1 = x1t[1:] - x1t[:-1:]
    dx2 = x2t[1:] - x2t[:-1]
    df = ft[1:] - ft[:-1]

    # diminuindo densidade das setas
    x1t_s = x1t[:-1:spacing]
    x2t_s = x2t[:-1:spacing]
    ft_s = ft[:-1:spacing]

    dx1_s = dx1[::spacing]
    dx2_s = dx2[::spacing]
    df_s = df[::spacing]

    # pontos
    idx = np.arange(0, len(dx1), spacing)
    dot_x1 = x1t[:-1][idx]
    dot_x2 = x2t[:-1][idx]
    dot_f = ft[:-1][idx]


    # domínio
    x = np.linspace(-2, 2, 25)
    y = np.linspace(-2, 2, 25)

    # cria grade 2D
    x1_mesh, x2_mesh = np.meshgrid(x, y)

    # função
    f = bi_dim_func_np(x1_mesh,x2_mesh)

    # superfície
    ax.plot_wireframe(x1_mesh, x2_mesh, f, color='green')

    # setas 3D
    ax.quiver(
        x1t_s, x2t_s, ft_s,
        dx1_s, dx2_s, df_s,
        color='black', length=0.5, normalize=True
    )
    # ponto inicial e final
    ax.scatter(metrics.x[0][0], metrics.x[0][1], metrics.f[0], color='blue', s=100)
    ax.scatter(dot_x1[1:-1],dot_x2[1:-1],dot_f[1:-1], color='cyan', s=100) # Cor diferentes para a transição
    ax.scatter(metrics.x[-1][0], metrics.x[-1][1], metrics.f[-1], color='red', s=100)

    ax.set_xlabel("x1")
    ax.set_ylabel("x2")
    ax.set_zlabel("f(x1,x2)")
    ax.set_zlim(-5, 5)

    ax.set_title(title + f" | It:{len(metrics.x)}")


def main():
    # a)
    result_min1,metrics_min1 = grad_descent_algorithm(
        bi_dim_func, # -> Função explode, tive que olhar ela para ver onde estão os minimos e tentar guiar para cair entre um dos dois mínimos sem a função explodir
        numerical_grad_op,
        grad_descent,
        [1.0,1.0], # [1,1],[-0.8,-0.8]
        gmin=0.01, # 0.01,0.01
        k_max=10000,
        alfa=0.01, # 0.01,0.001
        dh=0.01,
    )
    print("\n\n")
    print(f"Mínimo 1 ->    Ponto:                                      x_1={metrics_min1.x[-1][0]}    x_2={metrics_min1.x[-1][1]}    f(x_1,x_2)={metrics_min1.f[-1]}")
    print(f"Mínimo 1 ->    Algorítmo de Descida de Gradiente:          {result_min1}")
    print(f"Mínimo 1 ->    Valor da função:                            {bi_dim_func(result_min1)}")
    print("\n\n")
    result_min2,metrics_min2 = grad_descent_algorithm(
        bi_dim_func, # -> Função explode, tive que olhar ela para ver onde estão os minimos e tentar guiar para cair entre um dos dois mínimos sem a função explodir
        numerical_grad_op,
        grad_descent,
        [-0.8,-0.8], # [1,1],[-0.8,-0.8]
        gmin=0.01, # 0.01,0.01
        k_max=10000,
        alfa=0.001, # 0.01,0.001
        dh=0.01,
    )
    print(f"Mínimo 2 ->    Ponto:                                      x_1={metrics_min2.x[-1][0]}    x_2={metrics_min2.x[-1][1]}   f(x_1,x_2)={metrics_min2.f[-1]}")
    # if metrics:
    #     for i, (x, f, df) in enumerate(zip(metrics.x, metrics.f, metrics.df)):
    #         x_str = "[" + ", ".join(f"{v:+8.12f}" for v in x) + "]"
    #         df_str = "[" + ", ".join(f"{v:+8.12f}" for v in df) + "]"

    #         print(
    #             f"Iteração {i+1:6d} | "
    #             f"x={x_str} | "
    #             f"f(x)={f:+8.12f} | "
    #             f"f'(x)={df_str}"
    #         )
    print(f"Mínimo 2 ->    Algorítmo de Descida de Gradiente:          {result_min2}")
    print(f"Mínimo 2 ->    Valor da função:                            {bi_dim_func(result_min2)}")
    print("\n\n")

    # b)
    result_sympy,_ = grad_descent_algorithm(
        sympy_func,
        sympy_grad,
        sympy_grad_desc,
        [1.0,1.0], # [1,1],[-0.8,-0.8]
        gmin=0.01, # 0.01,0.01
        k_max=10000,
        alfa=0.01, # 0.01,0.001
        dh=0.01,
        # enable_metrics=False,
    )
    print(f"Mínimo 1 ->    Algorítmo de Descida de Gradiente (Sympy) : {result_sympy}")
    print(f"Mínimo 1 ->    Valor da função (Sympy):                    {bi_dim_func(result_sympy)}")
    print("\n\n")
    result_sympy,_ = grad_descent_algorithm(
        sympy_func,
        sympy_grad,
        sympy_grad_desc,
        [-0.8,-0.8], # [1,1],[-0.8,-0.8]
        gmin=0.01, # 0.01,0.01
        k_max=10000,
        alfa=0.001, # 0.01,0.001
        dh=0.01,
        # enable_metrics=False,
    )
    print(f"Mínimo 2 ->    Algorítmo de Descida de Gradiente (Sympy) : {result_sympy}")
    print(f"Mínimo 2 ->    Valor da função (Sympy):                    {bi_dim_func(result_sympy)}")
    print("\n\n")

    # c) Vou desenhar com os valores obtidos do gradiente aproximado!

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121, projection='3d')
    plot3d(metrics_min1,8,"Trajetória do Primeiro Mínimo",ax1)

    ax2 = fig.add_subplot(122, projection='3d')
    plot3d(metrics_min2,30,"Trajetória do Segundo Mínimo",ax2)

    #plt.savefig('trab1_part4.png', dpi=300, bbox_inches='tight')
    plt.show()


    

if __name__ == "__main__":
    main()

