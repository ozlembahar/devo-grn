
import pickle

from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell, derive_auc_threshold

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

import pandas as pd
import os
import seaborn as sns
import matplotlib as plt

from datetime import datetime
startTime = datetime.now()


ex_matrix = pd.read_csv("./data/expression_mat.csv", index_col=0)  # load count matrix 
tf_names = load_tf_names("./data/allTFs_dmel.txt") # Derive list of Transcription Factors(TF) for Drosophila

db_fnames = "./data/dm6-5kb-upstream-full-tx-11species.mc8nr.genes_vs_motifs.rankings.feather"
dbs = [RankingDatabase(fname=db_fnames, name=os.path.splitext(os.path.basename(db_fnames))[0])]

motif_annotation_file = "./data/motifs-v8-nr.flybase-m0.001-o0.0.tbl"



if not os.path.exists("results/"):
    os.mkdir("results")



n = 2  # TODO: decide on n based on computational resources 
all_results = [None] * n
for i in range(0, n): 
    run_num = i+1
    if not os.path.exists(f"results/run_{run_num}"): os.mkdir(f"results/run_{i+1}")

    """ phase 1 - GRN inference, generation of co-expression modules """
    adjacencies = grnboost2(ex_matrix, tf_names, verbose=True) # adjacencies table of tf, target and importance weight
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix)) # module generation - candidate regulons from TF-target gene interactions 
    # save to files:
    adjacencies.to_csv(f"results/run_{run_num}/adjacencies.csv", index=False, sep='\t')
    with open(f"results/run_{run_num}/modules.pkl", 'wb') as f:
        pickle.dump(modules, f)
    print("phase 1")
    """ phase 2+3 - Regulon prediction """
    df = prune2df(dbs, modules, motif_annotations_fname=motif_annotation_file) # Prune modules for targets with cis regulatory footprints (RcisTarget)
    regulons = df2regulons(df) # convert data frame to rergulons
    # save to files:
    df.to_csv(f"results/run_{run_num}/motifs.csv")
    with open(f"results/run_{run_num}/regulons.pkl", 'wb') as f:
        pickle.dump(regulons, f)
    print("phase 2+3")
    """ phase 4 - cellular enrichment """
    auc_mtx = aucell(ex_matrix, regulons, num_workers=1)  # Calculate enrichment of gene signatures for single cells. # TODO: change num_workers
    auc_mtx.to_csv(f"results/run_{run_num}/AUCell_mat.csv")
    # auc_mtx.to_pickle(f"results/results_{i+1}.pkl")# pickle results to results/ folder for later analysis
    # AUCell returns A dataframe with the AUCs (n_cells x n_modules).

    types_df = pd.read_csv("../data/cell_type.csv", index_col=0) 

    # heatmap
    lut = dict(zip(types_df.type.unique(), sns.color_palette("hls", len(types_df.type.unique()))))
    cell_colors = types_df.type.map(lut)
    row_colors = auc_mtx.merge(cell_colors, how='left', left_index=True, right_index=True).type
    ax= sns.clustermap(auc_mtx, figsize=(12,12,),yticklabels=True, xticklabels=True, row_colors=row_colors)
    ax.savefig(f"results/run_{run_num}/AUCell_heatmap.png")

print(datetime.now() - startTime)