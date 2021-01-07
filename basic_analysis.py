# Importo los módulos necesarios
import pandas as pd
import glob
import re
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
from utils import parse_file, primera_grafica, result_respir, num_clases_protein
from utils import num_clases_hydro, promedio_patron_protein, promedio_patron_hydro
from utils import clases_condicional

# Pandas me da el siguiente error:
# A value is trying to be set on
# a copy of a slice from a DataFrame.
# Try using .loc[row_indexer,col_indexer] = value instead
# disable chained assignments
pd.options.mode.chained_assignment = None

if __name__ == "__main__":

    # Load dataset
    # Carga ficheros "txt"
    # obtengo una lista con todos los txt contenidos en un directorio
    txt_files = glob.glob("data/orfs/*.txt")


    def read_lists():
        """
      Parsear un texto dado por un archivo
      Parameters
      ----------
      Returns
      -------
      data : lista de distas
          Parsed data
      """
        for i in txt_files:
            with open(i) as file:
                sublist = []
                previous_line = ''
                # Genero un patron donde solo me busque la lineas que contengan
                # las palabras "begin, funciton, tb_protein, tb_to_tb_evalue, end"
                pattern = re.compile(
                    r'\bbegin\b | \bfunction\b | \btb_protein\b | \btb_to_tb_evalue\b | \bend\b',
                    flags=re.I | re.X)
                for line in file:
                    # realizo transformaciones
                    line = line.replace('.', '')
                    line = line.replace('(', ";(")
                    line = line.replace('begin;(model', "begin_model")
                    line = line.replace('end;(model', "end_model")
                    line = line.replace('))', ")")
                    line = line.replace(')', "")
                    line = line.replace('(', "")
                    line = line.strip()
                    # Le indico donde comienza y termina cada sublista
                    if line.startswith(
                            'begin_model'
                    ) and previous_line.startswith(
                        'end_model'
                    ):
                        yield sublist
                        sublist = []
                    # Le informo que debe buscar las palabras del patros para
                    # añadir ilas a la sublista
                    if pattern.search(line) != None:
                        # voy añadiendo los datos a la sublista
                        sublist.append(line.split(';'))
                    # sublist.append(line)
                    previous_line = line
                # obtengo la sublista
                return sublist


    # creo una lista vacia para recolectar los datos
    data = []

    # Realizo una iteración para generar una lista de listas donde en
    # cada lista contiene la sublista generadas anteriormente
    for sublist in read_lists():
        data.append(sublist)


    # convierto la lista en un dataframe
    def appen_data(data):
        """
      Con esta función lo que hago es convertir las listas de la función
      anterior en un dataframe.
      :param data:
      :return:
      """

        def countList(data):
            """
          Función destinada a obtener el número de listas generadas
          :param data: la lista generada anteriormente
          :return: el número de listas
          """
            return len(data)

        # Genero una función para unir las columnas con
        # mismo nombre se unan y se separen por ';'
        def sjoin(x):
            """
          Función generada para fusionar las columnas con mismo nombre
          se unan incorporando el detalle de cada una de ellas
          y se separen por ";"
          :param x: columna a unir
          :return: los valores de las columnas unidas
          """
            return ';'.join(x[x.notnull()].astype(str))

        # Creo un dataframe vacio para agrupar los dataframes
        # que generamos con las listas
        df_a = pd.DataFrame()
        # genero una iteración donde convierto las listas en un dataframe
        for i in range(countList(data)):
            df_n = pd.DataFrame(data[i]).set_index(0)
            # pivoto el dataframe
            df_n = df_n.transpose()
            # uno las columnas con el mismo nombre
            df_n = df_n.groupby(level=0, axis=1).apply(
                lambda x: x.apply(sjoin, axis=1))
            df_a = df_a.append(df_n)
        return df_a


    txt_data = appen_data(data)

    # selecciono las columnas que me interesa
    data_relac_txt = txt_data.iloc[:, [1, 2, 0]]

    # las renombro
    data_relac_txt.columns = ['ORF', 'relacion_ORF', 'function']

    # Separo la columna función en otras
    data_relac_txt[['clase',
                    'name_gem',
                    'descripcion']] = data_relac_txt['function'].str.split(
        ",'", expand=True)

    # Genero nuevo dataframe seleccionando las columnas que me interesan
    df_total = data_relac_txt[['clase',
                               'name_gem',
                               'descripcion',
                               'ORF',
                               'relacion_ORF']]
    # Elimino valores na
    df_total = df_total.dropna()
    # genero nueva columna donde indico el número de relaciones
    # sumando el número de "tb" que hay
    df_total['numero de relaciones'] = df_total['relacion_ORF'].str.findall(
        r'tb').str.len()

    # Carga fichero "tb_functions.pl"
    file = 'data/tb_functions.pl'

    # Genero el dataframe con el fichero 'tb_fuctions.plt'
    df = parse_file('data/tb_functions.pl')

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
