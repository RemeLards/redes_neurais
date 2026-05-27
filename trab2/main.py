import numpy as np
import pandas as pd

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
    confusion_matrix,
    ConfusionMatrixDisplay,
) 
import matplotlib.pyplot as plt


FAULT_CHOSEN=7
SIMULATION_RUNS=20


def run_cst(
    cstr_sim_config: dict
) -> CSTR:

    cstr = CSTR(**cstr_sim_config)
    cstr.open()
    cstr.run()
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
  cstr_normal_list: list[CSTR],
  cstr_faulty_list:  list[CSTR],    
) -> pd.DataFrame:
    df_final = pd.DataFrame()

    for normal in cstr_normal_list:
        df = pd.read_csv(normal.datafn, sep=';')
        df_final = pd.concat([df_final,df],ignore_index=True)
    
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
    y_pred = x @ w
    loss = np.sum((y - y_pred) ** 2)/x.shape[0]
    return y_pred,loss


def train_classifier(
    dataset: dict,
    learning_rate: float,
    epochs: int
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
    w = np.zeros((x_train_with_bias.shape[1],1)) # bias included
    train_loss_seq = []
    val_loss_seq = []
    train_size = x_train_with_bias.shape[0]
    loss = 0
    while epoch_n <= epochs:
        # print(f"y_pred = {y_pred.shape}")
        y_pred,loss = loss_calc(w,x_train_with_bias,y_train)
        # print(f"loss = {loss.shape}")
        Lw = -2*(y_train - y_pred).T @ x_train_with_bias
        Lw /= train_size # Normalizar pelo quantidade do lote
        # print(f"Lw = {Lw.shape}")
        w -= learning_rate*Lw.T

        train_loss_seq.append(loss)
        _,loss_val = loss_calc(w,x_val_with_bias,y_val)
        val_loss_seq.append(loss_val)

        if epoch_n % 100 == 0:
            print(f"Época {epoch_n} - Loss: {loss:.4f}")
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
) -> None:
    epochs_range = range(1, len(train_loss) + 1)
    ax.plot(epochs_range, train_loss, label='Treino', color='blue', linewidth=2)
    ax.plot(epochs_range, val_loss, label='Validação', color='red', linewidth=2)

    ax.set_title('Loss durante o Treinamento', fontsize=14, fontweight='bold')
    ax.set_xlabel('Épocas', fontsize=12)
    ax.set_ylabel('Erro Quadrático Médio (MSE)', fontsize=12)

    ax.legend(fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    

def test_classifier(
    ax,
    dataset: dict,
    model: np.ndarray
) -> tuple[np.ndarray,list,list]:
    x_test = dataset["test"][0]
    y_test = dataset["test"][1]

    bias_test = np.ones((x_test.shape[0], 1))
    x_test_with_bias = np.concatenate((bias_test, x_test), axis=1) 
    y_pred,_ = loss_calc(model,x_test_with_bias,y_test)
    
    y_pred_classes = (y_pred >= 0.5).astype(int) # Já que  0 = Normal, 1 = Faulty -> score

    acuracia = accuracy_score(y_test, y_pred_classes)
    recall = recall_score(y_test, y_pred_classes)
    f1 = f1_score(y_test, y_pred_classes)
    cm = confusion_matrix(y_test, y_pred_classes)
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['Normal', 'Falha'])
    disp.plot(ax=ax, cmap='Blues', values_format='d', colorbar=False)
    ax.set_title('Matriz de Confusão', fontsize=12, fontweight='bold')

    print(f"Acurácia:       {acuracia:.4f}")
    print(f"Recall:         {recall:.4f}")
    print(f"F1-Score:       {f1:.4f}")




    
def main():
    exp =  {
        'theta': 1,
        'randseed': 1234,
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
            DELAY=1,
            TC=0.01
            ),
        )
    cstr_faulty_list = generate_runs(exp,SIMULATION_RUNS)
    pd_data = concat_pds(
        cstr_normal_list=cstr_normal_list,
        cstr_faulty_list=cstr_faulty_list,
    )
    data_std = data_standartization(data=pd_data)
    dataset = data_stratification_split(data=data_std)
    model,train_seq,val_seq = train_classifier(dataset=dataset,learning_rate=0.01,epochs=2000)

    _, ax = plt.subplots(1, 2, figsize=(14, 6))
    plot_loss(ax[0],train_loss=train_seq,val_loss=val_seq)
    test_classifier(ax[1],dataset=dataset,model=model)

    plt.tight_layout()
    #plt.savefig('resultado_modelo.png', dpi=300, bbox_inches='tight')
    plt.show()







if __name__ == "__main__":
    main()