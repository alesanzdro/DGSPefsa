import pandas as pd
import plotly.figure_factory as ff

def create_heatmap_matrix(df):
    unique_genes = set()
    max_alleles = {}
    for mlst_str in df['MLST']:
        mlst_pairs = mlst_str.split(';')
        for pair in mlst_pairs:
            gen, alelo = pair.split('(')
            alelo = int(alelo.replace(')', ''))
            unique_genes.add(gen)
            if gen not in max_alleles or alelo > max_alleles[gen]:
                max_alleles[gen] = alelo
    unique_genes = sorted(unique_genes)
    columns = [('muestra', '')] + [(gen, i+1) for gen in unique_genes for i in range(max_alleles[gen])]
    heatmap_data = []
    for _, row in df.iterrows():
        sample = row['Sample']
        mlst_str = row['MLST']
        data_row = {('muestra', ''): sample}
        mlst_pairs = mlst_str.split(';')
        for pair in mlst_pairs:
            gen, alelo = pair.split('(')
            alelo = int(alelo.replace(')', ''))
            data_row[(gen, alelo)] = 1.0
        heatmap_data.append(data_row)
    heatmap_df = pd.DataFrame(heatmap_data, columns=pd.MultiIndex.from_tuples(columns))
    # Convertir NaN a 0
    heatmap_df = heatmap_df.fillna(0.0)
    return heatmap_df

data = {
    'Sample': ['17_STEC_11111', '17_STEC_10960'],
    'MLST': ['adk(3);fumC(2);gyrB(4)', 'adk(2);fumC(5)']
}

df = pd.DataFrame(data)
heatmap_df = create_heatmap_matrix(df)

for col in heatmap_df.columns:
    if col[0] != 'muestra':
        heatmap_df[col] = heatmap_df[col].astype(float)

colorscale = [
    [0, 'lightgrey'], # Color para el valor 0
    [1, 'purple']     # Color para el valor 1
]

fig = ff.create_annotated_heatmap(
    z=heatmap_df.iloc[:, 1:].values,
    x=[f"{col[0]}_{col[1]}" for col in heatmap_df.columns[1:].tolist()],
    y=heatmap_df[('muestra', '')].tolist(),
    colorscale=colorscale,
    xgap=3,
)

fig.update_layout(
    xaxis=dict(
        tickvals=list(range(len(heatmap_df.columns[1:]))),
        ticktext=[f"{col[0]}_{col[1]}" for col in heatmap_df.columns[1:].tolist()],
        tickangle=-45,
        tickmode='array'
    ),
    yaxis=dict(
        tickvals=list(range(len(heatmap_df[('muestra', '')]))), 
        ticktext=heatmap_df[('muestra', '')].tolist()
    )
)

fig.show()
