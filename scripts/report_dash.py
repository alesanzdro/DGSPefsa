import os
import io
import dash
import numpy as np
import pandas as pd
import plotly.figure_factory as ff
from dash import html, dcc, dash_table, Output, Input, State
from datetime import datetime
from flask import send_file
from pandas.errors import ParserError


# conda activate dgsp_efsa_dash
# cd /software/DGSPefsa/scripts
# python report_dash.py
# Crear una función para procesar la columna MLST
def procesar_mlst(mlst):
    # Si mlst es igual a "-", se deja sin cambios
    if mlst == "-":
        return mlst

    # Dividir la cadena por ";"
    genes = mlst.split(';')
    nuevos_genes = []

    # Iterar a través de los genes
    for gen in genes:
        # Dividir cada gen por "("
        partes = gen.split('(')
        gen_nombre = partes[0]
        alelos = partes[1].rstrip(')').split('-')

        # Filtrar alelos mayores de 25
        alelos = [alelo for alelo in alelos if int(alelo) <= 25]

        # Reconstruir el gen con alelos separados por ","
        nuevo_gen = f"{gen_nombre}({','.join(alelos)})"
        nuevos_genes.append(nuevo_gen)

    # Unir los genes nuevamente por ";"
    mlst_modificado = ';'.join(nuevos_genes)

    return mlst_modificado



# Primero, definimos la función que procesará cada valor en la columna 'Sample'
def process_sample(sample):
    if 'CNP' in sample:
        parts = sample.split('_')
        if len(parts) in [2, 3]:
            return int('20' + parts[1][-2:])
        else:
            return None  # o cualquier valor predeterminado que desees usar
    else:
        return int('20' + sample.split('_')[0])



app = dash.Dash(__name__, suppress_callback_exceptions=True)
app.title = 'DGSP'
app._favicon = ("assets/logo.ico")

# Leer el archivo CSV
#df = pd.read_csv('csvfile.csv', sep=',', decimal='.', dtype={"Completeness": float})

# Definir la carpeta donde se encuentran los archivos CSV
folder_path = 'RECO'

# Crear una lista para almacenar los dataframes individuales
dataframes = []

# Iterar a través de los archivos en la carpeta
for filename in os.listdir(folder_path):
    if filename.endswith('.csv') and '_summary_' in filename:
        try:
            # Leer el archivo CSV
            predf = pd.read_csv(os.path.join(folder_path, filename), sep=',', decimal='.')
            # Aplicar la función procesar_mlst a la columna MLST

            predf['MLST'] = predf['MLST'].apply(procesar_mlst)

            # Dividir el nombre del archivo para obtener 'Run' y 'Date_results'
            parts = filename.replace('.csv', '').split('_summary_')
            if len(parts) == 2:
                Date_results, Run = parts
                # Crear un nuevo dataframe con las columnas 'Run' y 'Date_results'
                predf['Run'] = Run
                predf['Date_results'] = pd.to_datetime(Date_results, format='%y%m%d')
                
                # Agregar el dataframe a la lista
                dataframes.append(predf)
        except ParserError as e:
            print(f"Error al procesar el archivo '{filename}': {str(e)}")
        except Exception as e:
            print(f"Error general al procesar el archivo '{filename}': {str(e)}")

# Combinar todos los dataframes en uno solo
if dataframes:
    df = pd.concat(dataframes, ignore_index=True)
else:
    df = pd.DataFrame()  # En caso de que no se encuentren archivos

# Extraemos los años de la columna 'Sample'
#df['Year'] = ('20' + df['Sample'].str.split('_').str[0]).astype(int)
df['Year'] = df['Sample'].apply(process_sample)


app.layout = html.Div([
    html.H1("RESULTADOS SECUENCIACIÓN DGSP", style={'text-align': 'center'}),
    html.H2(datetime.now().strftime('%Y-%m-%d'), style={'text-align': 'center'}),
    html.Div([
        html.Button("Descargar EXCEL selección", id="btn_xlsx"),
        dcc.Download(id="download-dataframe-xlsx"),
    ], style={'float': 'right'}),
    # Añadimos una entrada de texto para la búsqueda
    html.Div([
        dcc.Input(
            id='search-box',
            type='text',
            placeholder='Introduce términos de búsqueda...',
            style={
                'backgroundColor': '#add8e6',  # Cambiar el color de fondo a verde claro
                'width': '40%',  # Ajustar el ancho del cuadro de entrada
                'padding': '10px'  # Agregar un espacio interno para hacerlo un poco más grande
            }
        ),
    ]),
    html.Br(),
    dcc.RangeSlider(
        id='year-slider',
        min=df['Year'].min(),
        max=df['Year'].max(),
        value=[df['Year'].min(), df['Year'].max()],
        marks={i: str(i) for i in range(df['Year'].min(), df['Year'].max()+1)}
    ),
    html.Div(id='slider-output-container', style={'float': 'right'}),
    html.Br(),
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
        style_as_list_view=True
        #export_format='xlsx',  # Permitir exportar como Excel
    ),
    html.Div(id='output_div'),

    html.Br(),
    html.Br(),
    html.H1('Heatmap MLST con muestras seleccionadas', style={'text-align': 'center'}),
    html.Br(),
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
        multi=False,  # permite seleccionar múltiples valores
        style={'width': '30%', 'float': 'left', 'margin-right': '10px'}
        ),
        html.Br(),
        html.Br(),
        dcc.Graph(id='heatmap')
    ]),
], style={'backgroundColor': '#f0f0f0'})


@app.callback(
    Output('slider-output-container', 'children'),
    Input('year-slider', 'value')
)
def update_slider_output(value):
    return f'Seleción: {value[0]} - {value[1]}'

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


@app.callback(
    Output("download-dataframe-xlsx", "data"),
    Input("btn_xlsx", "n_clicks"),  # Esto será el único que desencadene la callback
    State("year-slider", "value"),  # Esto sólo proporcionará su valor actual cuando se desencadene la callback
    State("search-box", "value"),  # Esto sólo proporcionará su valor actual cuando se desencadene la callback
    prevent_initial_call=True,
)
def func(n_clicks, selected_range, search_value):
    # Asegúrate de que selected_range y search_value están recibiendo valores desde las entradas correctas
    filtered_df = df[(df['Year'] >= selected_range[0]) & (df['Year'] <= selected_range[1])]
    
    if search_value:
        search_terms = search_value.split(' ')
        for term in search_terms:
            filtered_df = filtered_df[filtered_df.apply(lambda row: row.astype(str).str.contains(term, case=False).any(), axis=1)]

    filename = datetime.now().strftime("%y%m%d") + "_dgsp_export.xlsx"
    return dcc.send_data_frame(filtered_df.to_excel, filename, sheet_name=datetime.now().strftime("%y%m%d"))


# Callback para actualizar el heatmap
@app.callback(
    Output('heatmap', 'figure'),
    [Input('sample-selection', 'value'), 
     Input('table', 'derived_virtual_data')]
)
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

    # Filtrar el DataFrame para seleccionar las muestras deseadas y eliminar aquellas con MLST="-"
    filtered_df = table_df[table_df['Sample'].str.contains(selected_samples, case=False, na=False) & (table_df['MLST'] != "-")]

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
        annotation_text=[['' for _ in range(len(heatmap_df.columns) - 1)] for _ in range(len(heatmap_df))],
    )
    fig.update_layout(
        xaxis=dict(
            tickvals=list(range(len(heatmap_df.columns[1:]))),
            ticktext=[f"{col[0]}<br>{col[1]}" for col in heatmap_df.columns[1:].tolist()],
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
        ),
        shapes=[
            # Líneas de la cuadrícula horizontal
            dict(
                type='line',
                xref='paper',
                x0=0,
                x1=1,
                yref='y',
                y0=i-0.5,
                y1=i-0.5,
                line=dict(
                    color='darkgrey',
                    width=0.5
                )
            ) for i in range(1, heatmap_df.shape[0])
        ]
    )



    return fig



#figura=update_heatmap(selected_samples)


if __name__ == '__main__':
    app.run_server(debug=True)
