# Importo los módulos necesarios
import pandas as pd
import glob
import warnings
from utils import parse_file, carga_txt_files, clases_orf
from utils import primera_grafica, result_respir, num_clases_protein
from utils import num_clases_hydro, promedio_patron_protein
from utils import clases_condicional, gen_dataframe, promedio_patron_hydro
warnings.simplefilter(action='ignore', category=FutureWarning)


# Pandas me da el siguiente error:
# A value is trying to be set on
# a copy of a slice from a DataFrame.
# Try using .loc[row_indexer,col_indexer] = value instead
# disable chained assignments
pd.options.mode.chained_assignment = None

if __name__ == "__main__":
    # Load dataset
    # Genero el dataframe correspondiente a los archivos txt
    txt_files = glob.glob("data/orfs/*.txt")
    # aplico la función carga_txt_files
    df_total = carga_txt_files(txt_files)

    # Genero el dataframe con el fichero 'tb_fuctions.plt'
    df = parse_file('data/tb_functions.pl')

    # Imprimimos el listado de nº de ORFs que tiene cada clase
    clases_orf(df)

    # Graficamos el nº de ORFs que tiene cada clase
    primera_grafica(df)

    # Mostrar por pantalla cuántos ORFs pertenecen a la clase
    # que tiene Respiration como descripción
    print("El Nº de ORF incluidos en la clase con descripción "
          "'Respiration' es: {} \n\n".
          format(result_respir(df)))

    # Mostrar el número de clases que tiene 'protein' en
    # la descripción
    num_clases_protein(df)

    # Mostrar el número de clases que tiene palabras de
    # 13 caracteres y contiene 'hydro' en la descripción
    num_clases_hydro(df)

    # Obtener promedio ORFs con patron 'protein'
    promedio_patron_protein(df_total)

    # Obtener promedio ORFs con patron 'hydro'
    promedio_patron_hydro(df_total)

    # Numero de clases según condiciones indicadas
    clases_condicional(df)

    # Grafico el número de clases segun las condiciones
    gen_dataframe()

