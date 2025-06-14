import loompy as lp
import numpy as np
import scanpy as sc

x = sc.read_csv("count.csv")
row_attrs = {"Gene": np.array(x.var_names)}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("sample.loom", x.X.transpose(), row_attrs, col_attrs)