import numpy as np
import pandas as pd
import math
from copy import deepcopy

from CSTRSIM.python.CSTR import(
    CSTR,
    Fault,
)
# from CSTRSIM.python.CSTR_plot import (
#     dataframe2sklearn,
# )
from copy import deepcopy
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score,
    recall_score,
    f1_score,
    precision_score,
    confusion_matrix,
    ConfusionMatrixDisplay,
) 
import matplotlib.pyplot as plt
from typing import Callable

FAULT_CHOSEN=7
SIMULATION_RUNS=3


def run_cst(
    cstr_sim_config: dict
) -> CSTR:

    cstr = CSTR(**cstr_sim_config)
    cstr.open()
    cstr.run(log_run=False)
    cstr.close()
    return cstr


def generate_runs(
    cstr_base_sim_config: dict,
    amount: int
) -> list[CSTR]:
    if not cstr_base_sim_config['faults']:
        cstr_base_sim_config['id'] = "normal"
    else:
        fault_n = cstr_base_sim_config['faults'][0].id
        cstr_base_sim_config['id'] = f"fault_{fault_n}"

    cstr_normal_list = []
    for i in range(amount):
        cstr_sim_config = deepcopy(cstr_base_sim_config)
        cstr_sim_config["id"]+=f"_{i}"
        cstr_sim_config["randseed"]+=i
        cstr_normal_list.append(run_cst(cstr_sim_config))
    return cstr_normal_list


def data_standartization(
    data: pd.DataFrame,
) -> pd.DataFrame:
    x = data.drop(columns=['CLASS'])
    y = data['CLASS']

    # Inicializando Scaler no Scikit para fazer a Estandartização
    scaler = StandardScaler()

    x_std = scaler.fit_transform(x)
    x_std_df = pd.DataFrame(x_std, columns=x.columns)
    df_std = x_std_df.copy()
    df_std['CLASS'] = y.values

    # Remove linhas com valor Nulo
    df_std = df_std.dropna()

    return df_std


def data_stratification_split(
    data: pd.DataFrame        
) -> dict:
    x = data.drop(columns=['CLASS'])
    y = data['CLASS']


    x_train, x_temp, y_train, y_temp = train_test_split(
        x, y, 
        test_size=0.30,
        stratify=y, 
        random_state=100,
        shuffle=True
    )

    x_val, x_test, y_val, y_test = train_test_split(
        x_temp, y_temp, 
        test_size=0.50,
        stratify=y_temp, 
        random_state=100,
        shuffle=True
    )

    datasets = {
        "train" : [x_train.to_numpy(),y_train.to_numpy().reshape(-1, 1)],
        "validation" : [x_val.to_numpy(),y_val.to_numpy().reshape(-1, 1)],
        "test" : [x_test.to_numpy(),y_test.to_numpy().reshape(-1, 1)]
    }

    return datasets


def concat_pds(
  cstr_normal_list: list[CSTR] | None,
  cstr_faulty_list:  list[CSTR] | None,    
) -> pd.DataFrame:
    df_final = pd.DataFrame()

    if cstr_normal_list:
        for normal in cstr_normal_list:
            df = pd.read_csv(normal.datafn, sep=';')
            df_final = pd.concat([df_final,df],ignore_index=True)
    if cstr_faulty_list:
        for faulty in cstr_faulty_list:
            df = pd.read_csv(faulty.datafn, sep=';')
            df_final = pd.concat([df_final,df],ignore_index=True)

    # Varíaveis de restrição utilizadas no simulador
    # não são importante para o classificador pois não são medidas sensoriais e sim medidas do processo (restrições)

    COLUNMS_TO_DROP = ["r1","r2","r3","r4"]
    df_final.columns = [col.strip().replace('{','').replace('}','').replace('$','').replace('_','') for col in df_final.columns]
    df_final = df_final.drop(columns=COLUNMS_TO_DROP)

    # Troca as strings por números binários
    class_mapping = {
        "normal": 0.0,
        f"{FAULT_CHOSEN}": 1.0,
    }
    df_final['CLASS'] = df_final['CLASS'].map(class_mapping)

    # Salvando CSV
    print("Saving dataset on ./output/X.csv ...")
    df_final.to_csv('./output/X.csv', index=False)
    print("Dataset saved on ./output/X.csv")


    return df_final


def loss_calc(
    w: np.ndarray,
    x: np.ndarray,
    y: np.ndarray
) -> tuple[np.ndarray,float]:
    # print(f"W shape {w.shape}")
    # print(f"X shape {x.shape}")
    y_pred = x @ w
    loss = np.sum((y - y_pred) ** 2,axis=0)/x.shape[0]
    # print(f"Loss: {loss}")
    loss_norm = np.linalg.norm(loss)
    return y_pred,loss_norm


def train_classifier(
    dataset: dict,
    learning_rate: float,
    epochs: int,
) -> tuple[np.ndarray,list,list]:
    x_train = dataset["train"][0]
    y_train = dataset["train"][1]
    y_val = dataset["validation"][1]
    x_val = dataset["validation"][0]

    bias_train = np.ones((x_train.shape[0], 1))
    x_train_with_bias = np.concatenate((bias_train, x_train), axis=1) 
    bias_val = np.ones((x_val.shape[0], 1))
    x_val_with_bias = np.concatenate((bias_val, x_val), axis=1) 

    epoch_n = 0
    w = np.zeros((x_train_with_bias.shape[1],y_train.shape[1])) # bias included
    train_loss_seq = []
    val_loss_seq = []
    train_size = x_train_with_bias.shape[0]
    loss = 0
    while epoch_n <= epochs:
        # print(f"y_pred = {y_pred.shape}")
        y_pred,loss = loss_calc(w,x_train_with_bias,y_train)
        # print(f"loss = {loss.shape}")
        dLw = -2*(y_train - y_pred).T @ x_train_with_bias
        dLw /= train_size # Normalizar pelo quantidade do lote
        # print(f"Lw = {Lw.shape}")
        w -= learning_rate*dLw.T

        train_loss_seq.append(loss)
        _,loss_val = loss_calc(w,x_val_with_bias,y_val)
        val_loss_seq.append(loss_val)

        if epoch_n % 100 == 0:
            print(f"Epoch {epoch_n} - Loss: {loss:.4f}")
        epoch_n +=1

        if epoch_n % 100 == 0:
            med = sum(train_loss_seq[-100:])/100
            diff_sum = 0
            for loss in train_loss_seq[-100:]:
                diff_sum+=(loss-med)**2
            
            if diff_sum <= 5e-9:
                print(f"Stopping Training, last 100 iterations had no better results!")
                break
    
    return w,train_loss_seq,val_loss_seq


def plot_loss(
    ax,
    train_loss: np.ndarray,
    val_loss: np.ndarray,
    text: str,
) -> None:
    epochs_range = range(1, len(train_loss) + 1)
    ax.plot(epochs_range, train_loss, label='Treino', color='blue', linewidth=2)
    ax.plot(epochs_range, val_loss, label='Validação', color='red', linewidth=2)

    ax.set_title(f'Loss durante o Treinamento ({text})', fontsize=14, fontweight='bold')
    ax.set_xlabel('Épocas', fontsize=12)
    ax.set_ylabel('Erro Quadrático Médio (MSE)', fontsize=12)

    ax.legend(fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    

def binary_pred(
    y: np.ndarray,
) -> np.ndarray:
    return (y >= 0.5).astype(int) # Já que  0 = Normal, 1 = Faulty -> score


def one_hot_pred(
    y: np.ndarray,
) -> np.ndarray:
    softmax = []
    for y_array in y:
        den = 0
        score = []
        for y_val in y_array:
            den+=math.e**y_val
        for y_val in y_array:
            score.append((math.e**y_val)/den)
        softmax.append(np.array(score))

    pred = []
    for score in softmax:
        pred.append(np.argmax(score)) 
    
    return np.array(pred)

    
    

def test_classifier(
    ax,
    dataset: dict,
    model: np.ndarray,
    pred: Callable,
    text: str,
) -> tuple[np.ndarray,list,list]:
    x_test = dataset["test"][0]
    y_test = dataset["test"][1]

    bias_test = np.ones((x_test.shape[0], 1))
    x_test_with_bias = np.concatenate((bias_test, x_test), axis=1) 
    y_pred,_ = loss_calc(model,x_test_with_bias,y_test)
    
    y_pred_classes = pred(y_pred)
    y_test = pred(y_test) # its easier to just convert both into one variable classification

    acuracia = accuracy_score(y_test, y_pred_classes)
    recall = recall_score(y_test, y_pred_classes)
    precision = precision_score(y_test, y_pred_classes)
    f1 = f1_score(y_test, y_pred_classes)
    cm = confusion_matrix(y_test, y_pred_classes)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['Normal', 'Falha'])
    disp.plot(
        ax=ax,
        cmap='Blues',
        values_format='d',
        colorbar=True
    )
    ax.set_title(f'Matriz de Confusão ({text})', fontsize=12, fontweight='bold')
    ax.set_ylabel("")
    ax.set_xlabel(
        f"Acurácia={acuracia:.3f} | "
        f"Precisão={precision:.3f} | "
        f"Recall={recall:.3f} | "
        f"F1={f1:.3f}"
    )

    print(f"Acurácia:       {acuracia:.4f}")
    print(f"Recall:         {recall:.4f}")
    print(f"F1-Score:       {f1:.4f}")
    print(f"Precision:      {precision:.4f}")


def to_one_hot(
    dataset: dict
) -> None:
    n_classes = 2
    for dataset_type in dataset:
        y_array = dataset[dataset_type][1]
        y_array = y_array.astype(np.int8)
        dataset[dataset_type][1] = np.eye(n_classes)[y_array]
        dataset[dataset_type][1] = dataset[dataset_type][1].astype(np.float32).reshape((y_array.shape[0],n_classes))
        # print(f"Dataset {dataset[dataset_type][1]} ")


  
def main():
    exp =  {
        'theta': 1,
        'randseed': 1240,
        'fortran_rand': False,
        'timehoriz': 100,
        'faults': ()
    }
    cstr_normal_list = generate_runs(exp,SIMULATION_RUNS)
    # Escolhi a falha 7 para classificar se há falha no sistema
    # EXTENT0 -> Valor nominal da falha
    # DELAY -> tempo em minutos após o início da operação em que a falha ocorre
    # TC (Time Constant) -> Controla o quão rápido é a crescimento/descrimento da exponencial até o ponto de estabelecimento
    #

    exp['faults'] = (
        Fault(
            id=FAULT_CHOSEN,
            EXTENT0=10000,
            DELAY=1.0,
            TC=0.01
            ),
        )
    cstr_faulty_list = generate_runs(exp,SIMULATION_RUNS)
    pd_data = concat_pds(
        cstr_normal_list=cstr_normal_list,
        cstr_faulty_list=cstr_faulty_list,
    )
    data_std = data_standartization(data=pd_data)
    dataset_bin = data_stratification_split(data=data_std)
    model_bin,train_seq,val_seq = train_classifier(dataset=dataset_bin,learning_rate=0.01,epochs=2000)
    _, ax = plt.subplots(3, 2, figsize=(14, 6))
    plot_loss(
        ax[0,0],
        train_loss=train_seq,
        val_loss=val_seq,
        text="Binary Classifier"
    )
    test_classifier(
        ax[0,1],
        dataset=dataset_bin,
        model=model_bin,
        pred=binary_pred,
        text="Binary Classifier"
    )

    dataset_onehot = deepcopy(dataset_bin)
    to_one_hot(dataset=dataset_onehot)
    model_onehot,train_seq,val_seq = train_classifier(dataset=dataset_onehot,learning_rate=0.01,epochs=2000)
    plot_loss(
        ax[1,0],
        train_loss=train_seq,
        val_loss=val_seq,
        text="One Hot Classifier"
    )
    test_classifier(
        ax[1,1],
        dataset=dataset_onehot,
        model=model_onehot,
        pred=one_hot_pred,
        text="One Hot Classifier"
    )

    ### Testando com falhas mais severas com variáveis diferentes
    exp['faults'] = (
        Fault(
            id=FAULT_CHOSEN,
            EXTENT0=20000,
            DELAY=10.0,
            TC=0.01
            ),
        )
    cstr_faulty_list = generate_runs(exp,SIMULATION_RUNS*5)
    pd_data = concat_pds(
        cstr_normal_list=cstr_normal_list,
        cstr_faulty_list=cstr_faulty_list,
    )
    data_std = data_standartization(data=pd_data)
    x = data_std.drop(columns=['CLASS'])
    y = data_std['CLASS']
    dataset_bin = {
        "test" : [x.to_numpy(),y.to_numpy().reshape(-1, 1)]
    }
    dataset_bin = data_stratification_split(data=data_std)
    dataset_onehot = deepcopy(dataset_bin)
    to_one_hot(dataset=dataset_onehot)
    test_classifier(
        ax[2,0],
        dataset=dataset_bin,
        model=model_bin,
        pred=binary_pred,
        text="Binary Classifier [Different Test Set]"
    )
    test_classifier(
        ax[2,1],
        dataset=dataset_onehot,
        model=model_onehot,
        pred=one_hot_pred,
        text="One Hot Classifier [Different Test Set]"
    )



    plt.tight_layout()
    #plt.savefig('resultado_modelo.png', dpi=300, bbox_inches='tight')
    plt.show()







if __name__ == "__main__":
    main()