
"""E1. Actividad Integradora 1 - Puntos 1, 2 y 3.

KMP para busqueda de genes y proteinas.
Manacher adaptado para palindromos de ADN.
Traduccion de codones y seis marcos de lectura.
Cambio de reading frame en ORF1ab.
No se usan funciones de busqueda de la libreria estandar.
"""

import os

DIR = os.path.dirname(os.path.abspath(__file__))
GENOMA_WUHAN = "SARS-COV-2-MN908947.3.txt"

GENES = {
    "M": "gen-M.txt",
    "S": "gen-S.txt",
    "ORF1AB": "gen-ORF1AB.txt"
}


def leer_secuencia(nombre):
    """Lee un archivo FASTA o texto plano."""
    with open(os.path.join(DIR, nombre), "r") as f:
        return "".join(
            l.strip() for l in f if not l.startswith(">")
        ).upper()


def calcular_lps(patron):
    """Calcula el arreglo LPS para KMP. O(m)."""
    m = len(patron)
    lps = [0] * m
    largo, i = 0, 1

    while i < m:
        if patron[i] == patron[largo]:
            largo += 1
            lps[i] = largo
            i += 1
        elif largo > 0:
            largo = lps[largo - 1]
        else:
            lps[i] = 0
            i += 1

    return lps


def kmp(texto, patron):
    """Busca todas las apariciones de un patron. O(n+m)."""
    n, m = len(texto), len(patron)

    if m == 0:
        return []

    lps = calcular_lps(patron)
    resultado = []
    i = j = 0

    while i < n:
        if texto[i] == patron[j]:
            i += 1
            j += 1
        elif j > 0:
            j = lps[j - 1]
        else:
            i += 1

        if j == m:
            resultado.append(i - m)
            j = lps[j - 1]

    return resultado


def palindromo_mas_largo(s):
    """Manacher adaptado a palindromos de ADN. O(n)."""
    complemento = {
        "A": "T", "T": "A",
        "C": "G", "G": "C"
    }

    n = len(s)
    P = [0] * n
    izquierda = 0
    derecha = -1

    for i in range(n):
        if i > derecha:
            k = 0
        else:
            simetrica = izquierda + derecha - i + 1
            k = min(P[simetrica], derecha - i + 1)

        while i - k - 1 >= 0 and i + k < n:
            if complemento[s[i - k - 1]] != s[i + k]:
                break
            k += 1

        P[i] = k

        if i + k - 1 > derecha:
            izquierda = i - k
            derecha = i + k - 1

    max_indice = 0

    for i in range(n):
        if P[i] > P[max_indice]:
            max_indice = i

    inicio = max_indice - P[max_indice]
    longitud = 2 * P[max_indice]

    return inicio, longitud


# Tabla de codones del punto 3

CODONES = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


def traducir_adn(secuencia):
    """Convierte los codones de ADN en aminoacidos."""
    proteina = ""

    for i in range(0, len(secuencia) - 2, 3):
        codon = secuencia[i:i + 3]
        aminoacido = CODONES[codon]
        proteina += aminoacido

    return proteina


def leer_proteinas(nombre):
    """Lee las proteinas del archivo FASTA."""
    proteinas = {}
    nombre_proteina = None

    with open(os.path.join(DIR, nombre), "r") as f:
        for linea in f:
            linea = linea.strip()

            if linea.startswith(">"):
                nombre_proteina = linea[1:]
                proteinas[nombre_proteina] = ""

            elif nombre_proteina is not None:
                proteinas[nombre_proteina] += linea

    return proteinas


def reverso_complementario(secuencia):
    """Obtiene el reverso complementario del ADN."""
    complemento = {
        "A": "T", "T": "A",
        "C": "G", "G": "C"
    }

    resultado = ""

    for base in reversed(secuencia):
        resultado += complemento[base]

    return resultado


def marcos_lectura(genoma):
    """Traduce el genoma en sus seis marcos de lectura."""
    marcos = {}

    for i in range(3):
        marcos[f"+{i + 1}"] = traducir_adn(genoma[i:])

    reverso = reverso_complementario(genoma)

    for i in range(3):
        marcos[f"-{i + 1}"] = traducir_adn(reverso[i:])

    return marcos


def buscar_cambio_marco(genoma, proteina):
    """Busca el cambio de reading frame -1 en ORF1ab."""

    # Inicio de QHD43415_11 en base cero
    inicio = 13441

    for cambio in range(1, len(proteina)):

        posiciones = []

        for i in range(len(proteina)):
            if i < cambio:
                posicion = inicio + i * 3
            else:
                posicion = inicio + i * 3 - 1

            posiciones.append(posicion)

        traduccion = ""

        for posicion in posiciones:
            codon = genoma[posicion:posicion + 3]

            if len(codon) != 3:
                break

            traduccion += CODONES.get(codon, "?")

        if traduccion == proteina:
            fin = posiciones[-1] + 3

            codones = ""

            for posicion in posiciones[:4]:
                codones += genoma[posicion:posicion + 3]

            return inicio + 1, fin, codones, cambio

    return None


def punto1(genoma, secuencias):
    print("Punto 1: indices de aparicion de cada gen en el genoma:")

    for nombre, gen in secuencias.items():
        idx = kmp(genoma, gen)
        rango = [(i + 1, i + len(gen)) for i in idx]

        print(f"Gen {nombre}: indices (inicio-fin, base 1) = {rango}")
        print(f"  primeros 12 caracteres: {gen[:12]}")


def punto2(secuencias):
    print("\nPunto 2: palindromo mas largo por gen:")

    salida = []

    for nombre, gen in secuencias.items():
        ini, lon = palindromo_mas_largo(gen)
        pal = gen[ini:ini + lon]

        print(
            f"Gen {nombre}: longitud = {lon}, "
            f"posicion (base 1) = {ini + 1}-{ini + lon}"
        )

        salida.append(
            f"Gen {nombre}\n"
            f"longitud: {lon}\n"
            f"posicion: {ini + 1}-{ini + lon}\n"
            f"palindromo: {pal}\n"
        )

    with open(os.path.join(DIR, "palindromos.txt"), "w") as f:
        f.write("\n".join(salida))

    print("Guardado en palindromos.txt")


def punto3(genoma):
    """Busca las proteinas en los seis marcos de lectura."""

    print("\nPunto 3: proteinas en el genoma:")

    proteinas = leer_proteinas("seq-proteins.txt")
    marcos = marcos_lectura(genoma)

    total = 0
    reverso = reverso_complementario(genoma)

    for nombre, proteina in proteinas.items():
        encontrado = False

        for marco, aminoacidos in marcos.items():

            posiciones = kmp(aminoacidos, proteina)

            if len(posiciones) > 0:
                posicion = posiciones[0]

                desplazamiento = int(marco[1]) - 1
                longitud_adn = len(proteina) * 3

                if marco[0] == "+":
                    inicio = desplazamiento + posicion * 3
                    fin = inicio + longitud_adn
                    codones = genoma[inicio:inicio + 12]

                else:
                    fin = (
                        len(genoma)
                        - desplazamiento
                        - posicion * 3
                    )

                    inicio = fin - longitud_adn

                    posicion_reverso = (
                        desplazamiento + posicion * 3
                    )

                    codones = reverso[
                        posicion_reverso:posicion_reverso + 12
                    ]

                print(f"\nProteina: {nombre}")
                print(f"Indices (base 1): {inicio + 1}-{fin}")
                print(f"Marco de lectura: {marco}")
                print(f"Primeros 4 aminoacidos: {proteina[:4]}")
                print(f"Codones: {codones}")

                encontrado = True
                total += 1
                break

        # Caso especial de ORF1ab
        if not encontrado and nombre == "QHD43415_11":

            resultado = buscar_cambio_marco(genoma, proteina)

            if resultado is not None:
                inicio, fin, codones, cambio = resultado

                print(f"\nProteina: {nombre}")
                print(f"Indices (base 1): {inicio}-{fin}")
                print("Marco de lectura: +2 a +1")
                print("Cambio de reading frame: -1")
                print(f"Cambio despues de {cambio} aminoacidos")
                print(f"Primeros 4 aminoacidos: {proteina[:4]}")
                print(f"Codones: {codones}")

                encontrado = True
                total += 1

        if not encontrado:
            print(f"\nProteina: {nombre}")
            print("No encontrada.")

    print(
        f"\nProteinas encontradas: {total}/{len(proteinas)}"
    )


if __name__ == "__main__":
    genoma = leer_secuencia(GENOMA_WUHAN)

    secuencias = {
        n: leer_secuencia(a)
        for n, a in GENES.items()
    }

    punto1(genoma, secuencias)
    punto2(secuencias)
    punto3(genoma)
