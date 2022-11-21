"""Base functions for XClone clustering analysis.
"""

def XClustering(Xdata, Xlayer, PCA_comp = 30, n_clusters = 5, method = "k-means"):
    """
    """
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans
    from sklearn.cluster import AgglomerativeClustering
    import umap
    
    cell_num = Xdata.shape[0]
    X_use = Xdata.layers[Xlayer].reshape(cell_num, -1)
    
    ## PCA
    reduced_data = PCA(n_components=PCA_comp).fit_transform(X_use)
    
    ## Clustering
    if method == "k-means":
        kmeans = KMeans(init="k-means++", n_clusters=n_clusters, n_init=4)
        kmeans.fit(reduced_data)
        Z = kmeans.predict(reduced_data)
        
    if method == "Hierarchical":
        clustering = AgglomerativeClustering(n_clusters = n_clusters).fit(reduced_data)
        Z = clustering.labels_
    
    ## UMAP
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(reduced_data)
    
    
    return Z, embedding