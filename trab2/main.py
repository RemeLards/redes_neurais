from CSTRSIM.python.CSTR import(
    CSTR,
    Fault,
)
from CSTRSIM.python.CSTR_plot import (
    read_X,
    dataframe2sklearn,
)
from copy import deepcopy

def run_cst(
    cstr_sim_config: dict
) -> CSTR:

    cstr = CSTR(**cstr_sim_config)
    cstr.open()
    cstr.run()
    cstr.close()
    return cstr

def generate_normal_op(
    cstr_base_sim_config: dict,
    amount: int
) -> list[CSTR]:
    cstr_normal_list = []
    for i in range(amount):
        cstr_sim_config = deepcopy(cstr_base_sim_config)
        cstr_sim_config["id"]+=f"_{i}"
        cstr_sim_config["randseed"]+=i
        cstr_normal_list.append(run_cst(cstr_sim_config))
    return cstr_normal_list

    
    

def main():
    exp_normal =  {
        'id': 'normal100',
        'theta': 1,
        'randseed': 1234,
        'fortran_rand': False,
        'timehoriz': 100,
        #'plotmask': None,
        'faults': ()
    }
    cstr_list = generate_normal_op(exp_normal,2)
    for cstr in cstr_list:
        datafn = cstr.datafn
        print('Opening input data file %s.' % datafn)
        df = read_X(datafn=datafn)
        X, y, featname = dataframe2sklearn(df)
        labels = y



if __name__ == "__main__":
    main()