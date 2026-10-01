"""Pruebas rapidas de KMP y Manacher contra fuerza bruta."""
import random
from genoma import kmp, palindromo_mas_largo

random.seed(1)
for _ in range(2000):
    s = "".join(random.choice("ACGT") for _ in range(random.randint(1, 40)))
    p = "".join(random.choice("ACGT") for _ in range(random.randint(1, 3)))
    esperado = [i for i in range(len(s) - len(p) + 1) if s[i:i + len(p)] == p]
    assert kmp(s, p) == esperado
    mejor = max(len(s[i:j]) for i in range(len(s)) for j in range(i + 1, len(s) + 1) if s[i:j] == s[i:j][::-1])
    assert palindromo_mas_largo(s)[1] == mejor
print("OK: KMP y Manacher coinciden con fuerza bruta (2000 casos)")
