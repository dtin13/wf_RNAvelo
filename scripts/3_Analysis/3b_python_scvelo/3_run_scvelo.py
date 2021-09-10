############ Basic Model #################################################

sdata.var_names = genes
sdata.var_names

#Display spliced/unspliced proportions (10%-25% of unspliced recommended)
scv.utils.show_proportions(sdata)

#Preprocess the data
scv.pp.filter_genes(sdata, min_shared_counts=20)
scv.pp.normalize_per_cell(sdata)
scv.pp.filter_genes_dispersion(sdata, n_top_genes=2000)
scv.pp.log1p(sdata)

scv.pp.filter_and_normalize(sdata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(sdata, n_pcs=30, n_neighbors=30)

#Estimate RNA Velocity
scv.tl.velocity(sdata)
scv.tl.velocity_graph(sdata)

#Project the velocities
scv.pl.velocity_embedding_stream(sdata, basis='umap')
scv.pl.velocity_embedding(sdata, arrow_length=3, arrow_size=2, dpi=120)

#Identify important genes
scv.tl.rank_velocity_genes(sdata, groupby='clusters', min_corr=.3)

df = scv.DataFrame(sdata.uns['rank_velocity_genes']['names'])
df.head()

#Speed and coherence
scv.tl.velocity_confidence(sdata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(sdata, c=keys, cmap='coolwarm', perc=[5, 95])

df = sdata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

# Pseudotime (PAGA)
!pip install python-igraph --upgrade --quiet
# this is needed due to a current bug - bugfix is coming soon.
sdata.uns['neighbors']['distances'] = sdata.obsp['distances']
sdata.uns['neighbors']['connectivities'] = sdata.obsp['connectivities']

scv.tl.paga(sdata, groups='clusters')
df = scv.get_df(sdata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(sdata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)


############ Dynamical Model #################################################

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap')

#Kinetic rate parameters
df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata, 'fit*', dropna=True).head()

#Latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80)

#Top-likelihood genes
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False)
var_names = ['Actn4', 'Ppp3ca', 'Cpe', 'Nnat']
scv.pl.scatter(adata, var_names, frameon=False)
scv.pl.scatter(adata, x='latent_time', y=var_names, frameon=False)

#Cluster-specific top-likelihood genes
scv.tl.rank_dynamical_genes(adata, groupby='clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(5)


for cluster in ['Ductal', 'Ngn3 high EP', 'Pre-endocrine', 'Beta']:
    scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, frameon=False)





