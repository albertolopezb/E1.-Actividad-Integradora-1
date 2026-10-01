"""E1. Actividad Integradora 1 - Avance: puntos 1 (genes) y 2 (palindromos).

KMP para busqueda de patrones y Manacher para el
palindromo mas largo. No se usan funciones de busqueda de la libreria estandar.
"""
import os

DIR = os.path.dirname(os.path.abspath(__file__))
GENOMA_WUHAN = "SARS-COV-2-MN908947.3.txt"
GENES = {"M": "gen-M.txt", "S": "gen-S.txt", "ORF1AB": "gen-ORF1AB.txt"}


def leer_secuencia(nombre):
    """Lee un archivo (FASTA o texto plano) y regresa solo las bases."""
    with open(os.path.join(DIR, nombre), "r") as f:
        return "".join(l.strip() for l in f if not l.startswith(">")).upper()


def calcular_lps(patron):
    """Arreglo LPS (longest prefix suffix) del patron. O(m)."""
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
    """Todos los indices (base 0) donde aparece patron en texto. O(n+m)."""
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
    """Manacher: (inicio, longitud) del palindromo mas largo. O(n)."""
    # @ al inicio, # al final y $ entre caracteres: no pertenecen al alfabeto
    texto = "@$" + "$".join(s) + "$#"
    e = 2 * len(s) + 3
    P = [0] * e
    centro = limite = 0
    for i in range(1, e - 1):
        if i < limite:
            simetrica = 2 * centro - i
            P[i] = min(limite - i, P[simetrica])
        gap = P[i] + 1
        while texto[i - gap] == texto[i + gap]:
            P[i] += 1
            gap += 1
        if i + P[i] > limite:
            limite = i + P[i]
            centro = i
    max_indice = max(range(e), key=lambda k: P[k])
    inicio = (max_indice - P[max_indice] - 1) // 2
    return inicio, P[max_indice]


def punto1(genoma, secuencias):
    print("Punto 1: indices de aparicion de cada gen en el genoma: ")
    for nombre, gen in secuencias.items():
        idx = kmp(genoma, gen)
        rango = [(i + 1, i + len(gen)) for i in idx]  # 1-indexado, como NCBI
        print(f"Gen {nombre}: indices (inicio-fin, base 1) = {rango}")
        print(f"  primeros 12 caracteres: {gen[:12]}")


def punto2(secuencias):
    print("\nPunto 2: palindromo mas largo por gen: ")
    salida = []
    for nombre, gen in secuencias.items():
        ini, lon = palindromo_mas_largo(gen)
        pal = gen[ini:ini + lon]
        print(f"Gen {nombre}: longitud = {lon}, posicion (base 1) = {ini + 1}-{ini + lon}")
        salida.append(f"Gen {nombre}\nlongitud: {lon}\nposicion: {ini + 1}-{ini + lon}\npalindromo: {pal}\n")
    with open(os.path.join(DIR, "palindromos.txt"), "w") as f:
        f.write("\n".join(salida))
    print("Guardado en palindromos.txt")


if __name__ == "__main__":
    genoma = leer_secuencia(GENOMA_WUHAN)
    secuencias = {n: leer_secuencia(a) for n, a in GENES.items()}
    punto1(genoma, secuencias)
    punto2(secuencias)
