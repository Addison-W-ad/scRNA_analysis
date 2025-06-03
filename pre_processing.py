import scanpy as sc
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
import loompy as lp
import skmisc

def setup_settings():
    """Set up scanpy and matplotlib settings"""
    sc.settings.set_figure_params(
        dpi=300,
        dpi_save=300,
        frameon=False,
        color_map='viridis',
        format='png',
        facecolor='white',
        figsize=(5, 5)
    )
    sc.settings.verbosity = 1


def concat_data(data_dir:str,results_dir:str,samples:list,label:str,min_genes:int,min_cells:int):
    # Setup
    setup_settings()
    adatas = {}
    for sample in samples:
        print(f"Processing sample {sample}...")
        sample_path = os.path.join(data_dir, sample, "filtered_feature_bc_matrix")
        # Read data
        adata = sc.read_10x_mtx(sample_path, var_names='gene_symbols', cache=False)
        # Add sample ID
        adata.obs['sample'] = sample
        # Perform QC and filtering
        adata = perform_qc_and_filtering(adata,min_genes,min_cells)
        adatas[sample] = adata
    # Concatenate all samples
    combined = sc.concat(
        adatas.values(),
        keys=adatas.keys(),
        join='outer',
        label=label,
        index_unique='-')
    return combined

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
    np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def perform_qc_and_filtering(adata,min_genes:int,min_cells:int,log_path:str):
    """Perform QC metrics calculation and basic filtering"""
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl'))
    adata.var['hb'] = adata.var_names.str.contains("^Hb[^(p)]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], 
                             percent_top=[20,20,20], inplace=True)
    
    stats = [
        "\nQC metrics before filtering:",
        f"Number of cells: {adata.n_obs}",
        f"Number of genes: {adata.n_vars}",
        f"Median UMI counts per cell: {np.median(adata.obs.total_counts):.1f}",
        f"Median genes per cell: {np.median(adata.obs.n_genes_by_counts):.1f}", 
        f"Median MT content: {np.median(adata.obs.pct_counts_mt):.1f}%"
    ]
    
    with open(log_path, 'a') as f:
        f.write("\n".join(stats) + "\n")
    
    # Filter cells
    # sc.pp.filter_cells(adata, min_genes=20) 
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells) # 
    adata = adata[adata.obs.pct_counts_mt < 5].copy()  # Make a copy to avoid view
    adata = adata[adata.obs.pct_counts_ribo < 50].copy() # optional 
    # adata = adata[adata.obs.pct_counts_ribo < 50] # optional 
    # adata = adata[adata.obs.pct_counts_hb < 50]
    # remove outliers using 5MADs and 3 MADs for mitochondrial and total counts     

    # identify cells  more than 5 MADs (Median Absolute Deviations) away from the median.  
    adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5) 
    | is_outlier(adata, "log1p_n_genes_by_counts", 5)
    | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )

    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
    adata.obs["pct_counts_mt"] > 5
    )   # pct_counts_Mt is filtered with 3 MADs. Additionally, cells with a percentage of mitochondrial counts exceeding 5 % are filtered out.

    stats = [
        "\nAfter filtering:",
        f"Number of cells: {adata.n_obs}",
        f"Number of genes: {adata.n_vars}",
        f"Number of outliers: {adata.obs.outlier.sum() if 'outlier' in adata.obs.columns else 0}",
        f"Number of mt_outliers: {adata.obs.mt_outlier.sum() if 'mt_outlier' in adata.obs.columns else 0}",
        f"Total number of cells: {adata.n_obs}",
    ]
    
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    stats.append(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    with open(log_path, 'a') as f:
        f.write("\n".join(stats) + "\n")    
    return adata


def create_loom_file(adata, output_filename):
    """
    Create a loom file from the filtered AnnData object for SCENIC analysis
    
    Parameters:
    -----------
    adata : AnnData
        Filtered and processed AnnData object
    output_filename : str 
        Name of output loom file
    """
    import loompy as lp
    import os
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    
    # Create row attributes (genes)
    row_attrs = {
        "Gene": np.array(adata.var_names),
    }
    
    # Create column attributes (cells)
    col_attrs = {
        "CellID": np.array(adata.obs_names),
        "nGene": np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(),
        "nUMI": np.array(np.sum(adata.X.transpose(),axis=0)).flatten(),
    }
    
    # Create loom file
    lp.create(
        output_filename,
        adata.X.transpose(),
        row_attrs,
        col_attrs
    )


def preprocess_data(adata):
    adata.layers["counts"] = adata.X.copy()
    """Normalize and process the data"""
    # Normalize to depth
    sc.pp.normalize_total(adata, target_sum=1e6) # CPM normalization 1e6
    # Log transform
    sc.pp.log1p(adata)
    adata.layers["logCounts"] = adata.X.copy()
    # Calculate highly variable genes
    # sc.pp.highly_variable_genes(adata, n_top_genes=2000,flavor='seurat_v3',layer='counts') 
    #seurat_v3 is used for batch data variable and used the raw counts in the count layer
    # Scale data
    sc.pp.scale(adata, max_value=10) # above 10 or below -10 are set to 10 or -10
    return adata

def run_clustering_and_umap(adata):
    """Perform PCA, clustering and UMAP"""
    # Run PCA Singular Value Decomposition
    sc.pp.highly_variable_genes(adata, n_top_genes=2000,flavor='seurat_v3',layer='counts') 
    sc.tl.pca(adata, svd_solver='arpack',use_highly_variable=True)
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40) # preserve the local data i.e find out rare effect 
    # Run UMAP
    sc.tl.umap(adata)
    # Run Leiden clustering
    sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
    return adata

def find_markers(adata,index:str,resolution="leiden_res0_25",top_n=100):
    """Find marker genes for each cluster"""
    # Find marker genes
    sc.tl.rank_genes_groups(adata,resolution, method='wilcoxon')
    sc.tl.filter_rank_genes_groups(adata)
    
    # Create dataframe of top markers
    cluster_genes_filtered = {}
    for cluster in adata.obs[resolution].unique():
        genes_df = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        top_n_genes = genes_df.sort_values('scores', ascending=False).head(top_n)['names'].tolist()
        cluster_genes_filtered[f'Cluster_{cluster}'] = top_n_genes
    
    markers_df = pd.DataFrame(cluster_genes_filtered)
    markers_df.to_csv(f'{index}_cluster_markers_{resolution}_{top_n}.csv')
    return markers_df

def pcs_find(adata):
    """Find the optimal number of PCs"""
    # Run PCA Singular Value Decomposition
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    
def n_neighbors_find(adata,n_neighbors_range:list,n_pcs:int):

    fig, axes = plt.subplots(1, len(n_neighbors_range), figsize=(20, 4))

    for idx, n_neigh in enumerate(n_neighbors_range):
        # Copy the AnnData object to avoid modifying the original
        adata_temp = adata.copy()
        
    # Compute neighbors and UMAP
        sc.pp.neighbors(adata_temp, n_neighbors=n_neigh, n_pcs=40)
        sc.tl.umap(adata_temp)
    
    # Plot
        sc.pl.umap(adata_temp, 
               color='sample',  # or any other relevant column
               title=f'n_neighbors={n_neigh}',
               show=False,
               ax=axes[idx])

    plt.tight_layout()
    plt.show()

def bbknn_integration(adata):
    adata.X = adata.layers["logCounts"].copy()
    adata.raw = adata.copy()
    sc.pp.pca(adata)
    import bbknn
    bbknn.bbknn(
        adata, batch_key='Sample', neighbors_within_batch=3
    )
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.25, key_added="leiden_resolution_0.25")
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_resolution_0.5")
    sc.tl.leiden(adata, resolution=1, key_added="leiden_resolution_1")
    sc.tl.leiden(adata, resolution=1.5, key_added="leiden_resolution_1.5")
    sc.tl.leiden(adata, resolution=2, key_added="leiden_resolution_2")
    return adata

def run_leiden_clustering(adata, resolutions=[0.25, 0.5,0.75, 1,1.25 ,1.5,1.75, 2]):
    """
    Run Leiden clustering at multiple resolutions.
    
    Parameters:
    -----------
    adata : AnnData object
        Annotated data matrix
    resolutions : list
        List of resolution values to use for clustering
    """
    for res in resolutions:
        sc.tl.leiden(adata, 
                    resolution=res,
                    key_added=f'leiden_resolution_{res}')

# Run the clustering
