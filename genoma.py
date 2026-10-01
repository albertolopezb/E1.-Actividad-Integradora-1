"""E1. Actividad Integradora 1 - Avance: puntos 1 (genes) y 2 (palindromos).

Implementacion propia: KMP para busqueda de patrones y Manacher para el
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


def tabla_fallo(p):
    """Tabla de prefijos de KMP. O(m)."""
    fallo = [0] * len(p)
    k = 0
    for i in range(1, len(p)):
        while k and p[i] != p[k]:
            k = fallo[k - 1]
        if p[i] == p[k]:
            k += 1
        fallo[i] = k
    return fallo


def kmp(texto, patron):
    """Todos los indices (base 0) donde aparece patron en texto. O(n+m)."""
    if not patron:
        return []
    fallo, k, res = tabla_fallo(patron), 0, []
    for i, c in enumerate(texto):
        while k and c != patron[k]:
            k = fallo[k - 1]
        if c == patron[k]:
            k += 1
        if k == len(patron):
            res.append(i - k + 1)
            k = fallo[k - 1]
    return res


def palindromo_mas_largo(s):
    """Manacher: (inicio, longitud) del palindromo mas largo. O(n)."""
    t = "#" + "#".join(s) + "#"
    n = len(t)
    p = [0] * n
    centro = derecha = 0
    for i in range(n):
        if i < derecha:
            p[i] = min(derecha - i, p[2 * centro - i])
        while i - p[i] - 1 >= 0 and i + p[i] + 1 < n and t[i - p[i] - 1] == t[i + p[i] + 1]:
            p[i] += 1
        if i + p[i] > derecha:
            centro, derecha = i, i + p[i]
    r, c = max((v, i) for i, v in enumerate(p))
    return (c - r) // 2, r


def punto1(genoma, secuencias):
    print("=== Punto 1: indices de aparicion de cada gen en el genoma ===")
    for nombre, gen in secuencias.items():
        idx = kmp(genoma, gen)
        rango = [(i + 1, i + len(gen)) for i in idx]  # 1-indexado, como NCBI
        print(f"Gen {nombre}: indices (inicio-fin, base 1) = {rango}")
        print(f"  primeros 12 caracteres: {gen[:12]}")


def punto2(secuencias):
    print("\n=== Punto 2: palindromo mas largo por gen ===")
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
