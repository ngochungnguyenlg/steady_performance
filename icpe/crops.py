import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ruptures import Pelt
from typing import Callable, Tuple, List, Union
from functools import lru_cache
import sys

def pelt_for_crops(beta: float, algo : Pelt) -> Tuple[float, List[float]]:
    result = algo.predict(pen=beta)
    sum_cost = 0
    result_ = [0]+result.copy()
    for idx, seg in enumerate(result_, start=1):
        if idx < len(result_):
            sum_cost += algo.cost.error(result_[idx-1], result_[idx]-1)
    
    return sum_cost, result

class CropsClass:
    def __init__(self, method: Callable[[float], Tuple[float, List[float]]], betas: List[float], algo: Pelt):
        self.method = method
        self.betas = betas
        self.data = algo

def segmentations(crops: CropsClass) -> Union[pd.DataFrame, None]:
    print("Analysis Crops object starting ...")
    segs = []
    betas = []
    for beta in crops.betas:
        #[beta, (Qm, segs)]
        segs.append([beta[0], crops.method(beta[0], crops.data)])
        segs.append([beta[1], crops.method(beta[1], crops.data)])
        
    valid_segs = [seg for seg in segs if len(seg[1][1]) > 1]

    if not valid_segs:
        return None

    n = max(len(seg[1][1]) for seg in valid_segs)
    mat = np.array([seg[1][1] + [np.nan] * (n - len(seg[1][1])) for seg in valid_segs])
    col_names = [f"cpt.{i+1}" for i in range(n)]

    df = pd.DataFrame({
        'beta': [seg[0] for seg in valid_segs],
        'Qm': [seg[1][0] for seg in valid_segs],
        'Q': [seg[1][0] + seg[0] * len(seg[1][1]) for seg in valid_segs],
        'm': [len(seg[1]) for _, seg in valid_segs]
    })

    df = pd.concat([df, pd.DataFrame(mat, columns=col_names)], axis=1)
    print("Analysis Crops object ending ...")
    return df

def print_summary(crops: CropsClass):
    print("Crops analysis")
    print(f"Minimum penalty value = {min(crops.betas)} : Maximum penalty value = {max(crops.betas)}")

    segs_df = segmentations(crops)
    if segs_df is None:
        print(f"No segmentations found in the penalty interval [{min(crops.betas)}, {max(crops.betas)}]")
    else:
        print(f"Number of segmentations calculated: {len(segs_df)}")
        print(f"Least number of changepoints = {segs_df['m'].min()} : Maximum number of changepoints = {segs_df['m'].max()}")

def plot(crops: CropsClass, data: pd.DataFrame = None) -> None:
    df = segmentations(crops)
    if df is None:
        return

    fig, ax = plt.subplots()
    for idx, row in df.iterrows():
        # Vẽ các đường thẳng đứng từ trục x tại các giá trị trong row[4:]
        ax.vlines(x=row[4:], ymin=0, ymax=1, label=f"Beta = {row['beta']}")

    if data is not None:
        ax2 = ax.twinx()
        ax2.plot(data['x'], data['y'], 'r-')
        ax2.set_ylabel('Data Values')

    ax.set_xlabel('Location')
    ax.set_ylabel('Penalty')
    ax.legend()
    plt.show()

@lru_cache(None)
def crops(method: Callable[[float], Tuple[float, List[float]]],
          beta_min: float, beta_max: float,
          max_iterations: int = 20, algo: Pelt = None) -> CropsClass:
    
    def impl(f, beta_star, n):
        res = set()
        while beta_star and n > 0:
            sys.stdout.write("\rProcessing... Remaining iterations: {}".format(n))
            sys.stdout.flush()
            n -= 1
            try:
                beta = list(beta_star)[-1]
                beta_0, beta_1 = beta
                Qm_0, segs_0 = f(beta_0, algo)
                Qm_1, segs_1 = f(beta_1, algo)
                m_0 = len(segs_0)
                m_1 = len(segs_1)
                if m_0 > m_1 + 1:
                    beta_int = (Qm_1 - Qm_0) / (m_0 - m_1)
                    Qm_int, m_int = f(beta_int, algo)
                    if m_int != m_1:
                        beta_star = beta_star.union({(beta_int, beta_1), (beta_0, beta_int)})
                beta_star = beta_star.difference({beta})
                res.add(beta)
            except Exception as e:
                print("An error happen with {}".format(e))
        sys.stdout.write("\rDone! Processing completed.                    \n")
        sys.stdout.flush()
        return res

    def check_method_return_values(res):
        if not isinstance(res, tuple) or not isinstance(res[0], (int, float)) or not isinstance(res[1], list):
            raise ValueError("Method should return a tuple with a numeric value and a list.")

    check_crops_arguments(method, beta_min, beta_max, max_iterations)
    memoised_method = lru_cache(None)(method)
    CPT = memoised_method
    res = impl(CPT, {(beta_min, beta_max)}, max_iterations)
    return CropsClass(CPT, list(res), algo)

def check_crops_arguments(method, beta_min, beta_max, max_iterations):
    if not callable(method):
        raise ValueError("Method should be a callable function.")
    if not isinstance(beta_min, (int, float)) or beta_min < 0:
        raise ValueError("beta_min should be a positive numeric value.")
    if not isinstance(beta_max, (int, float)) or beta_max <= beta_min:
        raise ValueError("beta_max should be greater than beta_min.")
    if not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations should be a positive integer.")

def pelt_crops(data, beta_min, beta_max, max_iterations=20, jump= 5, min_size =2, method='l1'):
    algo = Pelt(model=method, min_size=min_size, jump=jump).fit(data)
    result = crops(pelt_for_crops, beta_min, beta_max, max_iterations, algo)
    # print_summary(result)
    data = np.array(data)
    # df_data = pd.DataFrame({'x': range(len(data)), 'y': data})
    # plot(result, data=df_data)
    df = segmentations(result)
    if str(type(df)).find("None")!=-1:
        pen = 15*np.log10(len(data))
        df = {pen:[0,len(data)]}
        return df
    df_dict = df.copy().transpose()
    df_dict.columns = df_dict.iloc[0]
    df_dict = df_dict.drop(index = df_dict.index[0:4])
    df_dict = df_dict.loc[:, ~df_dict.columns.duplicated()]
    df_dict = df_dict.to_dict(orient='list')
    cleaned_dict = {k: [v for v in lst if pd.notna(v)] for k, lst in df_dict.items()}
    return cleaned_dict
if __name__ == "__main__":
    N=100
    np.random.seed(1)
    data = np.concatenate([
        np.random.normal(loc=0, scale=1, size=N),   
        np.random.normal(loc=5, scale=0.5, size=N),  
        np.random.normal(loc=-3, scale=2, size=N),   
        np.random.normal(loc=10, scale=1, size=N),   
        np.random.normal(loc=2, scale=0.2, size=N),  
        np.random.normal(loc=-5, scale=3, size=N)   
    ])
    pelt_crops(data, 4, 15)
