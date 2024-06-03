import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import os, sys
from pathlib import Path
from pandarallel import pandarallel

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from ortools.sat.python import cp_model

########   Initialize and setup pandas methods   ########
pandarallel.initialize(nb_workers=os.cpu_count()-1, progress_bar=False, 
                       verbose=1, use_memory_fs=False) 
os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 

try:
    base_path = os.path.dirname(os.path.abspath(__file__))
except NameError:
    base_path = Path().resolve()
libs_path = os.path.join(base_path, '..', '..', 'libs')
libs_path = os.path.abspath(libs_path) 

if libs_path not in sys.path:
    sys.path.append(libs_path)

from scoring import Scoring
import pickle
from concurrent.futures import ProcessPoolExecutor
import os

import warnings
warnings.simplefilter('ignore')

def load_data():
    df_tp = pd.read_pickle('../TP/tp_prescore.pkl')
    df_tn = pd.read_pickle('../TN/tn_prescore.pkl')
    return df_tp, df_tn

def process_bootstrap(i):
    results = {}
    # For TP
    tp_val = df_tp_all.sample(frac=0.9, random_state=i)
    tp_test = df_tp_all.drop(tp_val.index)
    tp_val.to_pickle(f'pkls/tp_val_{i}.pkl')
    tp_test.to_pickle(f'pkls/tp_test_{i}.pkl')
    
    # For TN
    tn_val = df_tn_all.sample(frac=0.9, random_state=i)
    tn_test = df_tn_all.drop(tn_val.index)
    tn_val.to_pickle(f'pkls/tn_val_{i}.pkl')
    tn_test.to_pickle(f'pkls/tn_test_{i}.pkl')

    scoring_calibration = {}
    for index, solution in enumerate(all_solutions):
        if index % 200 == 0:
            print(f'set {i}: Solution {index + 1}')
        else:
            pass
        scoring_calibration[f"Solution {index + 1}"] = scoaring_calibraiton(
            solution, f'pkls/tp_val_{i}.pkl', f'pkls/tn_val_{i}.pkl')
    
    results[i] = scoring_calibration
    return results

class SolutionCollector(cp_model.CpSolverSolutionCallback):
    def __init__(self, variables):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__variables = variables
        self.__solutions = []

    def OnSolutionCallback(self):
        solution = {v.Name(): self.Value(v) for v in self.__variables}
        self.__solutions.append(solution)

    def GetAllSolutions(self):
        return self.__solutions

def find_all_solutions():
    # モデルを初期化
    model = cp_model.CpModel()

    # 変数の定義
    s = {i: model.NewIntVar(-10, 10, f's{i}') for i in range(1, 15)}

    # 制約の追加
    model.Add(s[3] == 0)                    #
    model.Add(s[1] + s[10] + s[14] == 9)    #
    model.Add(s[3] + s[11] + s[12] == 0)    #

    # Knowledge-based
    # Absolutes
    model.Add(s[2] >= 1)                #
    model.Add(s[1] - s[2] > 0)          #   
    
    # SpliceAI limitation
    # Absolutes
    model.Add(s[4] < 0)                 #
    model.Add(s[5] <= 0)                #
    model.Add(s[6] >= 0)                #
    # model.Add(s[7] >= 1)
    model.Add(s[12] < 0)                #
    model.Add(s[13] <= 0)               #
    model.Add(s[14] >= 0)               #
    # Order
    model.Add(s[5] - s[4] > 0)          #
    model.Add(s[6] - s[5] > 0)          #
    model.Add(s[7] - s[6] > 0)          #
    model.Add(s[13] - s[12] > 0)        #
    model.Add(s[14] - s[13] > 0)        #
    model.Add(s[4] >= s[12])            #

    # Absolutes
    model.Add(s[11] >= 0)               #
    # Score order
    model.Add(s[9] >= s[7])             #
    model.Add(s[8] - s[9] > 0)          #
    model.Add(s[10] - s[8] > 0)         #
    # model.Add(s[9] - s[11] > 0)         #
    model.Add(s[9] - s[11] >= 0)         #

    solver = cp_model.CpSolver()
    solution_collector = SolutionCollector([s[i] for i in range(1, 15)])
    solver.SearchForAllSolutions(model, solution_collector)
    
    return solution_collector.GetAllSolutions()

# Code below is adapted from Netflix's VMAF and BesenbacherLab's ROC-utils
# https://github.com/Netflix/vmaf/
# https://github.com/BesenbacherLab/ROC-utils
# Modifications: np.float -> np.float64

def compute_midrank(x):
    """Computes midranks.
    Args:
       x - a 1D numpy array
    Returns:
       array of midranks
    """
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=np.float64)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5*(i + j - 1)
        i = j
    T2 = np.empty(N, dtype=np.float64)
    # Note(kazeevn) +1 is due to Python using 0-based indexing
    # instead of 1-based in the AUC formula in the paper
    T2[J] = T + 1
    return T2


def fastDeLong(predictions_sorted_transposed, label_1_count):
    """
    The fast version of DeLong's method for computing the covariance of
    unadjusted AUC.
    Args:
       predictions_sorted_transposed: a 2D numpy.array[n_classifiers, n_examples]
          sorted such as the examples with label "1" are first
    Returns:
       (AUC value, DeLong covariance)
    Reference:
     @article{sun2014fast,
       title={Fast Implementation of DeLong's Algorithm for
              Comparing the Areas Under Correlated Receiver Operating Characteristic Curves},
       author={Xu Sun and Weichao Xu},
       journal={IEEE Signal Processing Letters},
       volume={21},
       number={11},
       pages={1389--1393},
       year={2014},
       publisher={IEEE}
     }
    """
    # Short variables are named as they are in the paper
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive_examples = predictions_sorted_transposed[:, :m]
    negative_examples = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty([k, m], dtype=np.float64)
    ty = np.empty([k, n], dtype=np.float64)
    tz = np.empty([k, m + n], dtype=np.float64)
    for r in range(k):
        tx[r, :] = compute_midrank(positive_examples[r, :])
        ty[r, :] = compute_midrank(negative_examples[r, :])
        tz[r, :] = compute_midrank(predictions_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def calc_pvalue(aucs, sigma):
    """Computes log(10) of p-values.
    Args:
       aucs: 1D array of AUCs
       sigma: AUC DeLong covariances
    Returns:
       log10(pvalue)
    """
    l = np.array([[1, -1]])
    z = np.abs(np.diff(aucs)) / np.sqrt(np.dot(np.dot(l, sigma), l.T))
    return np.log10(2) + scipy.stats.norm.logsf(z, loc=0, scale=1) / np.log(10)


def compute_ground_truth_statistics(ground_truth):
    assert np.array_equal(np.unique(ground_truth), [0, 1])
    order = (-ground_truth).argsort()
    label_1_count = int(ground_truth.sum())
    return order, label_1_count


def delong_roc_variance(ground_truth, predictions):
    """
    Computes ROC AUC variance for a single set of predictions
    Args:
       ground_truth: np.array of 0 and 1
       predictions: np.array of floats of the probability of being class 1
    """
    order, label_1_count = compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = predictions[np.newaxis, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count)
    assert len(aucs) == 1, "There is a bug in the code, please forward this to the developers"
    return aucs[0], delongcov


def delong_roc_test(ground_truth, predictions_one, predictions_two):
    """
    Computes log(p-value) for hypothesis that two ROC AUCs are different
    Args:
       ground_truth: np.array of 0 and 1
       predictions_one: predictions of the first model,
          np.array of floats of the probability of being class 1
       predictions_two: predictions of the second model,
          np.array of floats of the probability of being class 1
    """
    order, label_1_count = compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = np.vstack((predictions_one, predictions_two))[:, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count)
    return calc_pvalue(aucs, delongcov)

# Calculate AUC confidence interval (95%)
def compute_auc_confidence_interval(auc, var, confidence_level=0.95):
    alpha = 1 - confidence_level
    z_score = scipy.stats.norm.ppf(1 - alpha/2)  # 2-tailed z score
    se = np.sqrt(var)  # Calculate SE from variance
    lower_bound = auc - z_score * se
    upper_bound = auc + z_score * se
    return lower_bound, upper_bound


def scoaring_calibraiton(ths: dict, tp_pkl: str, tn_pkl: str):

    ths_scores = {'clinvar_same_pos': ths['s1'],
             'clinvar_same_motif': ths['s2'],
             'clinvar_else': ths['s3'],
             'non_canon_splai_lte_0.1_outside': ths['s4'],    
             'non_canon_splai_lte_0.1_other': ths['s5'],
             'non_canon_splai_bet_0.1_0.2': ths['s6'],
             'non_canon_splai_gte_0.2': ths['s7'],
             'canon_strong': ths['s8'], 
             'canon_moderate': ths['s9'], 
             'frameshift_nmd_eloF': ths['s10'], 
             'frameshift_nmd_not_eloF': ths['s11'],
             'canon_splai_lte_0.1': ths['s12'],
             'canon_splai_bet_0.1_0.2': ths['s13'],
             'canon_splai_gte_0.2': ths['s14']}

    scoring = Scoring(ths=ths_scores)
    # df = pd.read_pickle('../TP/tp_prescore_train.pkl')
    df = pd.read_pickle(tp_pkl)
    df['insilico_screening'] = df.parallel_apply(scoring.insilico_screening, axis=1)
    df['clinvar_screening'] = df.parallel_apply(scoring.clinvar_screening, axis=1)
    tp = scoring.calc_priority_score(df)

    # df = pd.read_pickle('../TN/tn_prescore_train.pkl')
    df = pd.read_pickle(tn_pkl)
    df['insilico_screening'] = df.parallel_apply(scoring.insilico_screening, axis=1)
    df['clinvar_screening'] = df.parallel_apply(scoring.clinvar_screening, axis=1)
    tn = scoring.calc_priority_score(df)
    # Exclude Y chromosome
    tn = tn[tn['CHROM'] != 'Y']

    tp.loc[:,'LABEL'] = 1
    tn.loc[:,'LABEL'] = 0
    tp = tp[['LABEL', 'PriorityScore', 'maxsplai']]
    tn = tn[['LABEL', 'PriorityScore', 'maxsplai']]

    # Combine tp and tn
    data = pd.concat([tp, tn], ignore_index=True)

    ground_truth = np.array(data['LABEL'])
    predictions_fw = np.array(data['PriorityScore'])

    auc1, var1 = delong_roc_variance(ground_truth, predictions_fw)
    cilower1, ciupper1 = compute_auc_confidence_interval(auc1, var1)

    return f"{auc1:.10f} [{cilower1:.12f}-{ciupper1:.12f}]"

################################################################################
# Set number of bootstrap and parallel processing 
start = sys.argv[1]
end = sys.argv[2]
bootstrap = int(end) - int(start) + 1
print(f"Bootstrap: {bootstrap} (start: {start}, end: {end})")
num_workers = os.cpu_count()-4
df_tp_all, df_tn_all = load_data()
all_solutions = find_all_solutions()
print(f'Total solutions found: {len(all_solutions)}')

# Output all_solutions as tsv
with open(f"all_solutions_{start}_{end}.tsv", 'w') as f:
    for solution in all_solutions:
        f.write(f"{solution}\n")

# for index, solution in enumerate(all_solutions):
#     print(f'Solution {index + 1}: {solution}')

if __name__ == '__main__':

    all_results: dict = {}
    start = int(start)
    end = int(end)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_results = list(executor.map(process_bootstrap, range(start, end + 1)))
    for result in future_results:
        all_results.update(result)

    # Save all_results object as pickle
    import pickle
    with open(f"all_results_{start}_{end}.pkl", 'wb') as f:
        pickle.dump(all_results, f)

    # load all_results object from pickle
    with open(f"all_results_{start}_{end}.pkl", 'rb') as f:
        all_results_restored = pickle.load(f)

    df = pd.DataFrame(columns=['Dataset', 'Candidate', 'Score'])
    for set in all_results_restored.keys():
        buf_df = pd.DataFrame(
            list(all_results_restored[set].items()), columns=['Candidate', 'Score'])
        buf_df['Dataset'] = set
        df = pd.concat([df, buf_df], ignore_index=True)

    df['auROC'] = df['Score'].str.split(' ', expand=True)[0].astype(float)
    df['95%CI'] = df['Score'].str.split(' ', expand=True)[1]
    df['CI_lower'] = df['95%CI'].str.extract(r'(\d\.\d*)').astype('float')
    df['CI_upper'] = df['95%CI'].str.extract(r'((?<=-)\d\.\d*(?=]))').astype(float)
    df.drop(['Score', '95%CI'], axis=1, inplace=True)

    solutions = pd.DataFrame(all_solutions)
    solutions_list = [f"Solution {i}" for i in range(1, 4851)]
    patterns_sr = pd.Series(solutions_list, name='Solution')
    solutions = pd.concat([patterns_sr, solutions], axis=1)

    # Merge the two dataframes
    df = pd.merge(df, solutions, left_on='Candidate', right_on='Solution')
    df.drop('Solution', axis=1, inplace=True)

    # Calculate the sample variance of the scores for each solution
    df['SampleVariance'] = df.iloc[:, 5:20].var(axis=1, ddof=0)

    # Group by dataset and extract the best auROC
    ddf = df.groupby('Dataset')
    maxdf = df.loc[ddf['auROC'].idxmax(),:]

    # Extract best auROC for each dataset
    # Extract solutions for each dataset with the best auROC value defined above
    # Store a dictionary
    set_max = []
    for i in range(start, end + 1):
        set_max.append(df.loc[(df['Dataset'] == i) & (df['auROC'] == maxdf.loc[maxdf['Dataset'] == i, 'auROC'].values[0]), :])

    # Save set_max as pickle
    with open(f"pkls/set_max_{start}.pkl", 'wb') as f:
        pickle.dump(set_max, f)

    # Extract the solution No. with most highest sample variance from each set_max[i]. 
    # The results put a dictionary with set number.
    best = {}
    for i, set_num in enumerate(range(start, end + 1)):
        highest_variance = set_max[i].loc[set_max[i]['SampleVariance'].idxmax(), 'SampleVariance']
        # If other candidate has the same variance, add it to the 
        best[f'set {set_num}'] = set_max[i].loc[set_max[i]['SampleVariance'] == highest_variance, 'Candidate'].values

    # Save best as pickle
    with open(f"pkls/best_{start}.pkl", 'wb') as f:
        pickle.dump(best, f)

