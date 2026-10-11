# E1. Actividad Integradora 1

Analisis del genoma de SARS-CoV-2 (Wuhan 2019 vs Texas 2020) con algoritmos de strings implementados desde cero (sin librerias de busqueda).

## Algoritmos utilizados
| Punto | Algoritmo | Complejidad |
|---|---|---|
| 1. Indices de cada gen en el genoma | KMP con arreglo LPS (clase 09) | O(n + m) |
| 2. Palindromo mas largo por gen | Manacher adaptado a palindromos de ADN (complemento inverso) | O(n) |
| 3. Ubicar cada proteina en el genoma | Tabla de codones + 6 marcos de lectura + KMP sobre aminoacidos; cambio de marco -1 para ORF1ab | O(n + p) por marco |
| 4a. Diferencias Wuhan vs Texas | Comparacion posicion a posicion y traduccion del codon afectado | O(n) |
| 4b. Cadena mas larga sin diferencias | Tramo mas largo entre diferencias (equivale a LCS cuando solo hay sustituciones) | O(n) |

`lcs_dp` implementa longest common substring con programacion dinamica y dos arreglos (clase 09), O(n*m) tiempo y O(m) espacio. Se usa para validar el punto 4b. Con los genomas completos da el mismo resultado: 9820 bases, indices 8783-18602.

## Resultados
Ver `capturas/salida_programa.txt` (salida completa del programa) y `palindromos.txt`.

- Punto 1: M 26523-27191, S 21563-25384, ORF1AB 266-21555.
- Punto 2: palindromo mas largo M = 10, S = 16, ORF1AB = 20.
- Punto 3: 24/24 proteinas encontradas.
- Punto 4: 7 diferencias (4 sinonimas, 3 que cambian el aminoacido); los genomas difieren ademas en 21 bases al final.

## Pruebas
`prueba_rapida.py` compara KMP, Manacher de ADN, LCS y el tramo sin diferencias contra fuerza bruta con casos aleatorios.

## Ejecucion
Colocar los archivos de Teams en esta carpeta y correr:
```
python genoma.py
python prueba_rapida.py
```
