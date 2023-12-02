from datetime import datetime
import pandas as pd
import dash
from dash import html, dcc, dash_table
from dash.dependencies import Input, Output
import plotly.figure_factory as ff
import os

app = dash.Dash(__name__)
app.title = 'DGSP'
app._favicon = ("assets/logo.ico")

# Leer el archivo CSV
df = pd.read_csv('csvfile.csv', sep=',', decimal='.', dtype={"Completeness": float})

# Extraemos los años de la columna 'Sample'
df['Year'] = ('20' + df['Sample'].str.split('_').str[0]).astype(int)

app.layout = html.Div([
    html.H1("RESULTADOS SECUENCIACIÓN DGSP", style={'text-align': 'center'}),
    html.H2(datetime.now().strftime('%Y-%m-%d'), style={'text-align': 'center'}),
    
    # Añadimos una entrada de texto para la búsqueda
    html.Div([
        dcc.Input(id='search-box', type='text', placeholder='Introduce términos de búsqueda...'),
    ], style={'text-align': 'center'}),
    
    dcc.RangeSlider(
        id='year-slider',
        min=df['Year'].min(),
        max=df['Year'].max(),
        value=[df['Year'].min(), df['Year'].max()],
        marks={i: str(i) for i in range(df['Year'].min(), df['Year'].max()+1)}
    ),
    html.Div(id='slider-output-container'), 
    dash_table.DataTable(
        id='table',
        #columns=[{'name': i, 'id': i} for i in df.columns],
        columns=[
            {'name': 'Run', 'id': 'Run', 'type': 'text'},
            {'name': 'Sample', 'id': 'Sample', 'type': 'text'},
            {'name': 'Sp_detected', 'id': 'Sp_detected', 'type': 'text'},
            {'name': 'Cov_Q30', 'id': 'Cov_Q30', 'type': 'numeric'},
            {'name': 'Completeness', 'id': 'Completeness', 'type': 'numeric'},
            {'name': 'Contamination', 'id': 'Contamination', 'type': 'numeric'},
            {'name': 'contigs', 'id': 'contigs', 'type': 'numeric'},
            {'name': 'N50_(contigs)', 'id': 'N50_(contigs)', 'type': 'numeric'},
            {'name': 'ST', 'id': 'ST', 'type': 'any'},
            {'name': 'MLST', 'id': 'MLST', 'type': 'text'},
            {'name': 'Antigenic_profile', 'id': 'Antigenic_profile', 'type': 'text'},
            {'name': 'ST_patho', 'id': 'ST_patho', 'type': 'text'},
            {'name': 'cgmlst_ST', 'id': 'cgmlst_ST', 'type': 'text'},
            {'name': 'Resistance_genes', 'id': 'Resistance_genes', 'type': 'text'},
            {'name': 'Antimicrobial(class)', 'id': 'Antimicrobial(class)', 'type': 'text'},
            {'name': 'Gene_mut(resistance)', 'id': 'Gene_mut(resistance)', 'type': 'text'},
            {'name': 'Status_Q30', 'id': 'Status_Q30', 'type': 'text'},
            {'name': 'Status_SNV', 'id': 'Status_SNV', 'type': 'text'},
            {'name': 'Status_Contamination', 'id': 'Status_Contamination', 'type': 'text'},
            {'name': 'Assembly_quality', 'id': 'Assembly_quality', 'type': 'text'},
            {'name': 'Date_results', 'id': 'Date_results', 'type': 'datetime', 'editable': False},
    
        ],
        data=df.to_dict('records'),
        page_size=25,  # Número de filas por página
        style_table={'overflowX': 'scroll', 'overflowY': 'scroll'}, # Barras de desplazamiento
        filter_action='native',  # Permitir filtrado por columnas
        sort_action='native',  # Permitir ordenar por columnas
        style_header={
            'backgroundColor': 'lightblue',
            'fontWeight': 'bold',
            'fontFamily': 'Calibri, Times New Roman'
        },
        style_data_conditional=[
            {
                'if': {'column_id': 'Status_Q30', 'filter_query': '{Status_Q30} eq "PASS"'},
                'backgroundColor': 'green',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Status_Q30', 'filter_query': '{Status_Q30} eq "FAIL"'},
                'backgroundColor': 'red',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Status_SP', 'filter_query': '{Status_SP} eq "PASS"'},
                'backgroundColor': 'green',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Status_SP', 'filter_query': '{Status_SP} eq "FAIL"'},
                'backgroundColor': 'red',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Status_SNV', 'filter_query': '{Status_SNV} eq "PASS"'},
                'backgroundColor': 'green',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Status_SNV', 'filter_query': '{Status_SNV} eq "FAIL"'},
                'backgroundColor': 'red',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Status_Contamination', 'filter_query': '{Status_Contamination} eq "PASS"'},
                'backgroundColor': 'green',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Status_Contamination', 'filter_query': '{Status_Contamination} eq "FAIL"'},
                'backgroundColor': 'red',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Assembly_quality', 'filter_query': '{Assembly_quality} eq "GOOD"'},
                'backgroundColor': 'green',
                'color': 'white'
            },
            {
                'if': {'column_id': 'Assembly_quality', 'filter_query': '{Assembly_quality} eq "BAD"'},
                'backgroundColor': 'red',
                'color': 'white'
            },
        ],
        style_data={
            'backgroundColor': 'rgba(173, 216, 230, {Completeness})',
            'border': '1px solid white'
        },
        style_as_list_view=True,
        export_format='csv',  # Permitir exportar como CSV
    ),
    html.Div(id='output_div'),
    html.Div([
        dcc.Dropdown(
            id='sample-selection',
            options=[
                {'label': 'LMON', 'value': 'LMON'},
                {'label': 'SALM', 'value': 'SALM'},
                {'label': 'STEC', 'value': 'STEC'},
                {'label': 'CAMP', 'value': 'CAMP'}
            ],
            value='LMON',  # valor inicial
            multi=False  # permite seleccionar múltiples valores
        ),
        dcc.Graph(id='heatmap')
    ]),


], style={'backgroundColor': '#f0f0f0'})


@app.callback(
    Output('slider-output-container', 'children'),
    Input('year-slider', 'value')
)
def update_slider_output(value):
    return f'Muestras según años: {value[0]} - {value[1]}'

@app.callback(
    Output('output_div', 'children'),
    Input('table', 'derived_virtual_data')
)
def update_output(rows):
    if rows is None:
        raise dash.exceptions.PreventUpdate
    else:
        return f"Número de muestras mostradas: {len(rows)}"

@app.callback(
    Output('table', 'data'),
    [Input('year-slider', 'value'), Input('search-box', 'value')]
)
def update_table(selected_range, search_value):
    filtered_df = df[(df['Year'] >= selected_range[0]) & (df['Year'] <= selected_range[1])]
    
    if search_value:
        search_terms = search_value.split(' ')
        for term in search_terms:
            filtered_df = filtered_df[filtered_df.apply(lambda row: row.astype(str).str.contains(term, case=False).any(), axis=1)]
    
    return filtered_df.to_dict('records')

# Callback para actualizar el heatmap
@app.callback(
    Output('heatmap', 'figure'),
    [Input('sample-selection', 'value'), 
     Input('table', 'derived_virtual_data')]
)

conda activate dgsp_efsa_dash
cd /software/DGSPefsa/scripts
python report_dash.py


def update_heatmap(selected_samples, table_data):
    if not selected_samples or not table_data:
        return {
            'data': [],
            'layout': {
                'title': 'No hay datos disponibles para el heatmap',
                'xaxis': {'title': 'Genes'},
                'yaxis': {'title': 'Muestras'}
            }
        }

    # Convertir los datos de la tabla a DataFrame
    table_df = pd.DataFrame(table_data)

    # Filtramos el DataFrame con las muestras seleccionadas (sin filtrar por "Completeness")
    filtered_df = table_df[table_df['Sample'].str.contains(selected_samples, case=False, na=False)]
    unique_genes = set()
    max_alleles = {}
    for mlst_str in filtered_df['MLST']:
        mlst_pairs = mlst_str.split(';')
        for pair in mlst_pairs:
            gen, alelo = pair.split('(')
            alelo = int(alelo.replace(')', ''))
            unique_genes.add(gen)
            if gen not in max_alleles or alelo > max_alleles[gen]:
                max_alleles[gen] = alelo

    unique_genes = sorted(unique_genes)
    # Aquí, he modificado el loop para que cree columnas desde 1 hasta el alelo máximo para cada gen
    columns = [('muestra', '')] + [(f"<b>{gen}</b>", i+1) for gen in unique_genes for i in range(max_alleles[gen])]

    heatmap_data = []
    for _, row in filtered_df.iterrows():
        sample = row['Sample']
        mlst_str = row['MLST']
        if pd.isnull(mlst_str):
            continue
        data_row = {('muestra', ''): sample}
        for col in columns[1:]:
            data_row[col] = 0.0

        mlst_pairs = mlst_str.split(';')
        for pair in mlst_pairs:
            gen, alelo = pair.split('(')
            alelo = int(alelo.replace(')', ''))
            data_row[(f"<b>{gen}</b>", alelo)] = 1.0
        heatmap_data.append(data_row)

    dict_keys = list(heatmap_data[0].keys())
    heatmap_df = pd.DataFrame(heatmap_data, columns=pd.MultiIndex.from_tuples(dict_keys))
    
    heatmap_df = heatmap_df.fillna(0.0)
    for col in heatmap_df.columns:
        if col[0] != 'muestra':
            heatmap_df[col] = heatmap_df[col].astype(float)

    colorscale = [
        [0, 'lightgrey'], # Color para el valor 0
        [1, 'purple']     # Color para el valor 1
    ]

    if heatmap_df.empty:
        raise ValueError("heatmap_df está vacío, no se puede crear un heatmap")

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
            ticktext=[f"{html.B(col[0])}<br>{col[1]}" for col in heatmap_df.columns[1:].tolist()],
            tickangle=-45,
            tickmode='array'
        ),
        yaxis=dict(
            tickvals=list(range(len(heatmap_df[('muestra', '')]))),
            ticktext=heatmap_df[('muestra', '')].tolist(),
            tickfont=dict(
                size=14,
                color='black',
                family='Arial Black, sans-serif'
            )
        )
    )
    return fig



#figura=update_heatmap(selected_samples)


if __name__ == '__main__':
    app.run_server(debug=True)
