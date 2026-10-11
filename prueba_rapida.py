"""Pruebas rapidas: cada algoritmo se compara contra fuerza bruta."""
import random
from genoma import (kmp, palindromo_mas_largo, lcs_dp, diferencias,
                    subcadena_comun_alineada)

COMPLEMENTO = {"A": "T", "T": "A", "C": "G", "G": "C"}


def aleatoria(n):
    return "".join(random.choice("ACGT") for _ in range(n))


random.seed(1)

for _ in range(2000):
    s = aleatoria(random.randint(1, 40))

    # KMP contra busqueda ingenua
    p = aleatoria(random.randint(1, 3))
    esperado = [i for i in range(len(s) - len(p) + 1) if s[i:i + len(p)] == p]
    assert kmp(s, p) == esperado

    # Palindromo de ADN (complemento inverso) contra fuerza bruta
    mejor = 0
    for i in range(len(s)):
        for j in range(i + 2, len(s) + 1, 2):
            sub = s[i:j]
            rc = "".join(COMPLEMENTO[c] for c in reversed(sub))
            if sub == rc and len(sub) > mejor:
                mejor = len(sub)
    ini, lon = palindromo_mas_largo(s)
    assert lon == mejor

print("OK: KMP y Manacher de ADN coinciden con fuerza bruta (2000 casos)")

for _ in range(500):
    # LCS (programacion dinamica) contra fuerza bruta
    a, b = aleatoria(random.randint(1, 25)), aleatoria(random.randint(1, 25))
    mejor = max(
        j - i
        for i in range(len(a)) for j in range(i + 1, len(a) + 1)
        if a[i:j] in b
    ) if any(c in b for c in a) else 0
    assert lcs_dp(a, b)[0] == mejor

    # Tramo mas largo sin diferencias contra LCS, con mutaciones puntuales
    g1 = aleatoria(random.randint(30, 80))
    g2 = list(g1)
    for k in random.sample(range(len(g1)), random.randint(0, 5)):
        g2[k] = random.choice([c for c in "ACGT" if c != g1[k]])
    g2 = "".join(g2)
    inicio, longitud = subcadena_comun_alineada(g1, g2, diferencias(g1, g2))
    assert longitud == lcs_dp(g1, g2)[0]

print("OK: LCS y tramo sin diferencias coinciden con fuerza bruta (500 casos)")
