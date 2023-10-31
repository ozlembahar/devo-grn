from glob import glob
import pandas as pd
import pickle
from collections import Counter 
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
from os.path import basename, normpath
from collections import Counter


def regulons_to_df(dir_glob) -> list:
    regulons_dfs = {}
    for folder in glob(dir_glob):
        run = basename(normpath(folder))
        with open(folder + 'regulons.pkl', 'rb') as f:
            regulon = pickle.load(f)
            reg_df = pd.DataFrame({
                'TF': [reg.transcription_factor for reg in regulon],
                'genes': [list(reg.genes) for reg in regulon],
                'score': [reg.score for reg in regulon]
            })
        reg_df['size'] = reg_df['genes'].apply(lambda x: len(x))
        reg_df = reg_df.set_index('TF')
        regulons_dfs[run] = reg_df
        
    return regulons_dfs


def dict_tocsv(df_dict, path) -> None:
    if not os.path.exists(path):
        os.makedirs(path)
    
    for run, df in df_dict.items():
        df.to_csv(path + run + ".csv")


def count_TFs(regulons_dict) -> pd.Series:
    tfs_counter = Counter
    for tf, df in regulons_dict.items():
        tfs_counter.update(df.index)
    
    return pd.Series(tfs_counter).sort_values(ascending=False)


def main():
    # results_dir_glob = "results/run_*/"
    # results_dir = "results/"

    results_dir_glob = './run_*/'
    results_dir = "./"


    # regulons to csv, including size of regulon
    regulons_dict = regulons_to_df(results_dir_glob)
    dict_tocsv(regulons_dict, path=results_dir + 'regulons/')
    
    # count TFs
    tfs_counter =count_TFs(regulons_dict=regulons_dict)
    # display count distributions
    f, (ax_box, ax_hist) = plt.subplots(2, sharex=True, gridspec_kw={"height_ratios": (.15, .85)})
    sns.boxplot(tfs_counter, ax=ax_box)
    sns.histplot(tfs_counter, ax=ax_hist, bins=50) # set bin number differently
    sns.despine(ax=ax_box, left=True)
    plt.axvline(40, color='red')
    plt.savefig("fig.png")




if __name__ == "__main__":
    try:
        working_dir = sys.argv[1]  # get working dir from user # Assuming project directory
        os.chdir(working_dir)
    except Exception:
        pass
    main()