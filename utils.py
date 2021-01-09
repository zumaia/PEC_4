import matplotlib.pyplot as plt
import re
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# Carga ficheros txt
def carga_txt_files(txt_files):
    """
    Función cargar, tranformar y convertir en dataframe
    los archivos txt
    Parameters
    ----------
    txt_files
    Returns
    -------
    el dataframe con los archivos txt
    """
    # Carga ficheros "txt"
    # obtengo una lista con todos los txt contenidos en un directorio

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
                # las palabras "begin, funciton, tb_protein, tb_to_tb_evalue,
                # end"
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
    return df_total


# Carga fichero "tb_functions.pl"
# file = 'data/tb_functions.pl'

def parse_file(file):
    """
    Parse text at given filepath
    Parameters
    ----------
    Filepath for file_object to be parsed
    Returns
    -------
    data : pd.DataFrame
        Parsed data
    """
    # Cargo el archivo de texto dentro de "content"
    with open(file) as f:
        content = f.readlines()

        content = [x.strip() for x in content]

        # creo la lista vaciacreate an empty list to collect the data
        events_cl = []
        # create an empty list to collect the data
        events_fc = []

        # this is to filter the lines in the data
        for line in content:
            # extracto la clase
            if 'class([' in line:
                line = line.replace('.', '')
                line = line.replace('"', '')
                line = line.replace(')', '')
                line = line.replace('],', '];')
                events_cl.append(line.split('('))
            # extracto la función
            if 'function(' in line:
                line = line.replace('").', '"')
                line = line.replace(',[', ';[')
                line = line.replace("],'", '];"')
                line = line.replace("',", '";')
                events_fc.append(line.split('('))

        # para el caso de las clases
        # convierto la lista en dataframe
        data_class = pd.DataFrame(events_cl)
        # divido la columna en 2 columnas nuevas
        data_class[['clase',
                    'desc']] = data_class[1].str.split(';',
                                                       1,
                                                       expand=True)
        # Selecciono las columnas necesarias
        data_class = data_class[['clase', 'desc']]
        # Para el caso de las funciones
        # convierto la lista en dataframe
        data_fc = pd.DataFrame(events_fc)
        # divido la columna en 2 columnas nuevas
        data_fc[['ORF',
                 'clase',
                 'name_gem',
                 'desc_ORF']] = data_fc[1].str.split(';',
                                                     3,
                                                     expand=True)
        # Selecciono las columnas necesarias
        data_fc = data_fc[['ORF', 'clase', 'name_gem', 'desc_ORF']]
        # los 2 dataframes
        df_merged = pd.merge(data_class, data_fc, on='clase', how='outer')

    return df_merged


def get_unique_values(df, col):
    """
    Return the number of unique values in the column `col`
        for the dataframe `df`.
    """
    return len(df[col].unique())


def grafico(x):
    """
    Función destinada a generar la gráfica para reproducir los
    ORFs que pertenecen a cada clase
    :param x: introducimos el valor
    :returns obtenemos la gráfica y un fichero
    """
    ax = x.plot(kind='barh', figsize=(8, 20), color='#86bf91',
                zorder=2, width=0.75)

    # Despine
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # Switch off ticks
    ax.tick_params(axis="both", which="both",
                   bottom="off", top="off", labelbottom="on",
                   left="off", right="off", labelleft="on")
    # Draw vertical axis lines
    vals = ax.get_xticks()
    for tick in vals:
        ax.axvline(x=tick, linestyle='dashed', alpha=0.3,
                   color='#eeeeee', zorder=1)

    [ax.text(v, i, '{}'.format(v)) for i, v in enumerate(x)]
    # Set titulo
    ax.set_title("Numero de ORFs por clase", fontsize=14, weight="bold")
    # Set x-axis label
    ax.set_xlabel("Número de ORF", labelpad=20, weight='bold', size=12)
    # Set y-axis label
    ax.set_ylabel("Clase",  labelpad=20, weight='bold', size=12)
    plt.savefig(f"./graficos/numero_ORF.png")
    plt.show()
    plt.close()


def df_dict(df):
    """
    Función destinada a convertir en un diccionario
    la suma de OFS agrupados or clase
    Parameters
    ----------
    df
    Returns
    -------
    diccionario
    """
    # Convierto en un diccionario
    g_c_count = df.groupby(['clase'])['ORF'].count().to_dict()
    # Imprimo para cada fila
    return g_c_count


def clases_orf(df):
    """
    Función destinada a listar el nº de ORFs que tiene
    cada clase
    Parameters
    ----------
    df
    Returns
    -------
    Listado de clases con el nº de ORFs
    """
    # Convierto en un diccionario
    g_c_count = df.groupby(['clase'])['ORF'].count().to_dict()
    # Imprimo para cada fila
    for c in g_c_count:
        print("La clase {clase} tiene {ORF} ORFs".format(clase=c,
                                                         ORF=g_c_count[c]))


def primera_grafica(df):
    """
    Función destinada a realizar un recuento de las ORFs
    que tiene cada clase y su resultado mostrarlo en una
    gráfica utilizando la función grafico
    :param df: dataframe
    :returns obtenemos la gráfica
    """
    # Agrupo las clases y cuento el número de ORF
    group_clase_count = df.groupby(['clase'])['ORF'].count()
    # Aplico la función al resultado.
    grafico(group_clase_count)


def result_respir(df):
    """
    Función para mostrar en pantalla el nº de ORF que
    pertenecen a la clase que tiene Respiration como
    descripcion
    :param df
    :returns resultado consulta
    """
    result_respiration = df[df.desc.str.contains(
        'Respiration', case=False)].sum()
    return result_respiration['ORF']


# Genero una función para mostrar en pantalla el
# resultado de este ejercicio
def plot_pie(x, names):
    """
    Genero una función donde obtengo en formato tarta el
    número de clases
    :param x
    :param names: el patrón correspondiente
    :return:
    """
    # eliminio el index
    x = x.reset_index()
    # convierto un dataframe solo con el campo seleccionado
    x = pd.DataFrame(x)['desc_ORF']
    x.index += 1
    x.plot.pie(figsize=(8, 8))
    plt.title(f"El Nº de clases con mínimo de 1 ORF con patron "
              f"${names}$", fontsize=14, weight="bold")
    plt.ylabel("Porcion clase sobre total")
    plt.xlabel("Total clases: {}".format(len(x)))
    # guardo la gráfoca en un lugar adecuado
    plt.savefig(f"./graficos/grafica_patron_{names}.png")
    # muestro en pantalla
    plt.show()
    plt.close()
    return plot_pie


def num_clases_protein(df):
    """
    Función para mostrar en pantalla el nº de clases
    con un mínimo de un ORF con el patron 'termino'
    proteina
    descripcion
    :param df: dataframe
    :returns resultado consulta
    """
    # Realizo la busqueda
    des_ORF_protein = df[df['desc_ORF'].str.contains(r'protein',
                                                     re.IGNORECASE,
                                                     regex=False, na=False)]
    # agrupo y cuento las clases
    protein = des_ORF_protein.groupby(['clase'])['desc_ORF'].count()
    # obtengo el número
    numero = len(protein)
    # imprimo el resultado
    print(
        "El número de clases que contienen como mínimo un ORF con el patron "
        "'que contiene el término protein'.  {} \n". format(numero))

    # aplico la función al resultado anterior
    plot_pie(protein, 'proteina')


def num_clases_hydro(df):
    """
    Función para mostrar en pantalla el nº de clases
    con un mínimo de un ORF con el patron 'termino'
    hydro y contiene una palabra de 13 caracteres en
    descripcion
    :param df: dataframe
    :returns resultado consulta
    """
    # Busco las palabras que contengan 13 caracteres
    ORF_hydro_13 = df[df['desc_ORF'].str.contains(r'\b\w{13}\b') == True]
    # dentro de ese resultado, que me ofrezca los que contienen
    # la palabra 'hydro'
    ORF_hydro_13 = ORF_hydro_13[ORF_hydro_13['desc_ORF'].str.contains(
                                                        r'hydro') == True]
    # agrupo y cuento las clases
    hydro = ORF_hydro_13.groupby(['clase'])['desc_ORF'].count()
    # obtengo el número
    numero = len(hydro)
    # imprimo el resultado
    print("El número de clases según el segundo patron 'La descripción"
          "contiene una palabra de 13 caracteres y esta contiene el "
          "término *hydro*' es :{} \n".format(numero))
    # aplico la función al resultado anterior
    plot_pie(hydro, 'hydro')


def numero_protein(df):
    """
    Función para obtener el nº de clases
    con un mínimo de un ORF con el patron 'termino'
    proteina
    descripcion
    :param df: dataframe
    :returns resultado consulta
    """
    # Realizo la busqueda
    des_ORF_protein = df[df['desc_ORF'].str.contains(r'protein',
                                                     re.IGNORECASE,
                                                     regex=False,
                                                     na=False)]
    # agrupo y cuento las clases
    protein = des_ORF_protein.groupby(['clase'])['desc_ORF'].count()
    # obtengo el número
    numeros = len(protein)
    return numeros


def numero_hydro(df):
    """
    Función para obtener el  nº de clases
    con un mínimo de un ORF con el patron 'termino'
    hydro y contiene una palabra de 13 caracteres en
    descripcion
    :param df
    :returns resultado consulta
    """
    # Busco las palabras que contengan 13 caracteres
    ORF_hydro_13 = df[df['desc_ORF'].str.contains(r'\b\w{13}\b') == True]
    # dentro de ese resultado, que me ofrezca los que contienen
    # la palabra 'hydro'
    ORF_hydro_13 = ORF_hydro_13[ORF_hydro_13['desc_ORF'].str.contains(
                                                        r'hydro') == True]
    # agrupo y cuento las clases
    hydro = ORF_hydro_13.groupby(['clase'])['desc_ORF'].count()
    # obtengo el número
    numero = len(hydro)
    return numero


# Genero la función para las gráficas de los promedios
def graf_relac(x, names):
    """
    Función destinada a generar las gráficas de los
    promedios de relacionessegún el patrón correspondiente
    :param x: el datafame de los ORF con
    el número de relaciones
    :param names: Indicamos el patron
    :return: la gráfica
    """
    x.plot.barh('ORF', 'numero de relaciones', figsize=(8, 10))
    # muestro en patalla la media
    plt.axvline(x['numero de relaciones'].mean(),
                color='r', linestyle='--', label='Promedio')
    # la leyenda en la parte superior derecha
    plt.legend(loc='upper right')
    # titulo
    plt.title(f"Nº de ORF relacionados con otro ORF y su promedio ${names}$",
              fontsize=14,
              weight="bold")
    # etiquetas
    plt.xlabel("Nº de ORF")
    plt.ylabel("ORF")
    # guardo en fichero png la gŕafica
    plt.savefig(f"./graficos/grafica_relaciones_patron_{names}.png")
    # muestro en pantalla
    plt.show()
    # cierro
    plt.close()
    return plot_pie


def promedio_patron_protein(df_total):
    """
    Función que generamos una grafica a partir del
    filtro según patrón "protein" en el apartado descripción
    del dataset
    :param: dataframe
    :return: Resultado impreso y la gráfica
    """
    # filtro los valores que contengan la palabra 'protein'
    df_total_protein = df_total[df_total['descripcion'].str.contains(
        r'protein', regex=False, na=False)]
    # selecciono 2 columnas
    prome_relac_protein = df_total_protein[['ORF', 'numero de relaciones']]
    # obtengo la media
    n_prome_relac_protein = float(
                    df_total_protein[['numero de relaciones']].mean())
    # muestro el resultado
    print("El número promedio de ORFs con los cuales se relacionan los"
          "ORFs con el patron 'que contiene el término protein'. {} \n".format(
            n_prome_relac_protein))
    # aplico la función a los datos obtenidos
    graf_relac(prome_relac_protein, 'proteina')


def promedio_protein(df_total):
    """
    Función que generamos una grafica a partir del
    filtro según patrón "protein" en el apartado descripción
    del dataset
    :param: dataframe
    :return: Resultado impreso y la gráfica
    """
    # filtro los valores que contengan la palabra 'protein'
    df_total_protein = df_total[df_total['descripcion'].str.contains(
        r'protein', regex=False, na=False)]
    # selecciono 2 columnas
    prome_relac_protein = df_total_protein[['ORF', 'numero de relaciones']]
    # obtengo la media
    n_prome_relac_protein = float(
                    df_total_protein[['numero de relaciones']].mean())
    # muestro el resultado
    return n_prome_relac_protein


def promedio_patron_hydro(df_total):
    """
    Función que generamos una grafica a partir del
    filtro según patrón "hydro" en el apartado descripción
    del dataset
    Parameters
    ----------
    df_total

    Returns
    -------
    Se obtiene el resultado del promedio y la gráfica correspondiente
    """
    df_total_hydro = df_total[df_total['descripcion'].str.contains(
                                            r'\b\w{13}\b') == True]
    df_total_hydro = df_total_hydro[df_total_hydro['descripcion'].str.contains(
                                                            r'hydro') == True]
    n_prome_relac_hydro = float(
                            df_total_hydro[['numero de relaciones']].mean())
    print("El número promedio de ORFs con los cuales se relacionan"
          "los ORFs según el segundo patron 'La descripción"
          "contiene una palabra de 13 caracteres y esta contiene el "
          "término *hydro*' es: {} \n".format(n_prome_relac_hydro))
    graf_relac(df_total_hydro, 'hydro')


def promedio_hydro(df_total):
    """
    Función que generamos un resultado mediio a partir del
    filtro según patrón "hydro" en el apartado descripción
    del dataset
    Parameters
    ----------
    df_total

    Returns
    -------
    Se obtiene el resultado del promedio y la gráfica correspondiente
    """
    df_total_hydro = df_total[df_total['descripcion'].str.contains(
                                            r'\b\w{13}\b') == True]
    df_total_hydro = df_total_hydro[df_total_hydro['descripcion'].str.contains(
                                                            r'hydro') == True]
    n_prome_relac_hydro = float(
                            df_total_hydro[['numero de relaciones']].mean())
    return n_prome_relac_hydro


# Genero una lista vacia
result = []


def clases_condicional(df):
    """
    Función destinada a obtener el número de clases que tienen
    como mínimo una dimensión mayor estricta (>) que 0 y a la
    vez múltiple de M entre 2 y 9.
    :parameter: dataframe
    :return:
    número de clases que tienen como mínimo una dimensión
    mayor estricta (>) que 0 y a la vez múltiple de M entre 2 y 9.
    """
    # Selecciono la columna clase del dataframe y realizo unas modificaciones
    df['clase'] = df['clase'].str.replace(
        '[', '').str.replace(']', '').str.replace("''", "")
    # obtengo una lista de las cláses únicas
    clase_uniq = df['clase'].unique().tolist()
    # Genero una lista vacia
    lista = []
    print("Para cada entero M entre 2 y 9 (ambos incluidos),"
          "calcula el número de clases que tienen como mínimo"
          "una *dimensión* mayor estricta (>) que 0 y a la vez"
          "múltiple de M.")
    # Realizo una iteración para generar una lista de listas
    # selecciono cada uno de los valores de las lista
    for j in range(len(clase_uniq)):
        # divido el elemento correspondiente
        x = clase_uniq[j].split(',')
        # convierto el str en int
        x = [int(i) for i in x]
        # añado a la nueva lista
        lista.append(x)
    # genero una iteración para obtener el resultado solicitado
    for i in range(2, 10):
        # any() usado para filtrar elementos
        # selecciono las que tienen al menos un 0
        lista_zero = [sub for sub in lista if(any(ele > 0 for ele in sub))]
        # para evitar que a la hora de ejecutar el siguiente condicionante
        # me selecciono el valor 0 ya que '0%i==0' y me lo contaría
        # sustituyo su valor por -1
        lista_zeros = [[-1 if x == 0 else x for x in y] for y in lista_zero]
        # seleccion las listas que es multiple por M
        res = [subs for subs in lista_zeros if(any(
                                ele % i == 0 for ele in subs))]
        # genero una lista con los datos
        a = len(res)
        result.append(a)
        # imprimo el resultado según el formato solicitado
        print("M = {}:  ".format(i) + str(len(res)) + " clases")


def numero_lista(df):
    """
    Función destinada a obtener el número de lista de clases
    :parameter: dataframe
    :return:
    numero de lista de clases
    """
    # Selecciono la columna clase del dataframe y realizo unas modificaciones
    df['clase'] = df['clase'].str.replace(
        '[', '').str.replace(']', '').str.replace("''", "")
    # obtengo una lista de las cláses únicas
    clase_uniq = df['clase'].unique().tolist()
    # Genero una lista vacia
    lista = []
    # Realizo una iteración para generar una lista de listas
    # selecciono cada uno de los valores de las lista
    for j in range(len(clase_uniq)):
        # divido el elemento correspondiente
        x = clase_uniq[j].split(',')
        # convierto el str en int
        x = [int(i) for i in x]
        # añado a la nueva lista
        lista.append(x)
    return len(lista)


def grafico_n_clases(d_result):
    """
    Función destinada a generar la gráfica para reproducir las
    clases según su multiplo
    :param d_result: introducimos el valor
    :returns obtenemos la gráfica y un fichero
    """
    ax = d_result.plot(kind='barh', figsize=(8, 8), color='#86bf91',
                       zorder=2, width=0.85)
    # Despine
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    # Switch off ticks
    ax.tick_params(axis="both", which="both",
                   bottom="off", top="off", labelbottom="on",
                   left="off", right="off", labelleft="on")
    # Draw vertical axis lines
    vals = ax.get_xticks()
    for tick in vals:
        ax.axvline(x=tick, linestyle='dashed', alpha=0.3,
                   color='#eeeeee', zorder=2)
    [ax.text(v, i+2, str(v)) for i, v in enumerate(result)]
    # Set titulo
    ax.set_title("Numero de clases por múltiplo", fontsize=14, weight="bold")
    # Set x-axis label
    ax.set_xlabel("nº de clases", labelpad=20, weight='bold', size=12)
    # Set y-axis label
    ax.set_ylabel("Multiplo", labelpad=20, weight='bold', size=12)
    plt.savefig(f"./graficos/clases_multiplo.png")
    plt.show()
    plt.close()


def gen_dataframe():
    """
    Función destinada a generar un dataframe para utilizar en la
    función posterior y poder obtener una gráfica donde se representa
    los multiplos y el nº de clases
    Parameters
    Returns
    dataframe
    """
    # genero una lista
    l_lista = [i for i in range(2, 10)]
    # creo un diccionario con dos listas tales que los multiplos
    # como claves y las listas como valores
    d = {'M': l_lista, 'clases': result}
    # convierto en dataframe
    d_result = pd.DataFrame(d).set_index('M')
    # Genero la gráfica
    grafico_n_clases(d_result)
