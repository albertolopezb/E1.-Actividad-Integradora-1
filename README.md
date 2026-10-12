# E1. Actividad Integradora 1

Analisis del genoma de SARS-CoV-2 (Wuhan 2019 vs Texas 2020) con algoritmos de strings implementados desde cero (sin librerias de busqueda).

## Algoritmos utilizados
| Punto | Algoritmo | Complejidad |
|---|---|---|
| 1. Indices de cada gen en el genoma | KMP con arreglo LPS (clase 09) | O(n + m) |
| 2. Palindromo mas largo por gen | Manacher adaptado a palindromos de ADN (complemento inverso) | O(n) |
| 3. Ubicar cada proteina en el genoma | Tabla de codones + 6 marcos de lectura + KMP sobre aminoacidos; cambio de marco -1 para ORF1ab | O(n + p) por marco |
| 4a. Diferencias Wuhan vs Texas | Verificacion de alineacion con KMP + comparacion posicion a posicion + traduccion del codon afectado | O(n) |
| 4b. Cadena mas larga sin diferencias | Tramo mas largo entre diferencias (equivale a LCS cuando solo hay sustituciones) | O(n) |

`lcs_dp` implementa longest common substring con programacion dinamica y dos arreglos (clase 09), O(n*m) tiempo y O(m) espacio. Se usa para validar el punto 4b y como respaldo si los genomas no estuvieran alineados. Con los genomas completos da el mismo resultado: 9820 bases, indices 8783-18602.

Antes de comparar posicion a posicion, `genomas_alineados` busca con KMP el inicio y el final de Texas dentro de Wuhan: ambos aparecen en las mismas coordenadas, asi que no hay inserciones ni borrados que recorran los indices. La unica diferencia de longitud son 21 A's extra en la cola poly-A de Wuhan (no codificante).

## Resultados
Ver `capturas/salida_programa.txt` (salida completa del programa), `palindromos.txt` y `cadena_sin_diferencias.txt`.

- Punto 1: M 26523-27191, S 21563-25384, ORF1AB 266-21555.
- Punto 2: palindromo mas largo M = 10, S = 16, ORF1AB = 20.
- Punto 3: 24/24 proteinas encontradas.
- Punto 4: 7 sustituciones (4 sinonimas, 3 que cambian el aminoacido: QHD43415_13 D379A, QHD43422 T11I y L84S) y 21 bases extra de poly-A al final de Wuhan. Cadena mas larga sin diferencias: 9820 bases (8783-18602).

## Pruebas
`prueba_rapida.py` compara KMP, Manacher de ADN, LCS y el tramo sin diferencias contra fuerza bruta con casos aleatorios, y comprueba que la verificacion de alineacion detecta indels.

## Ejecucion
Colocar los archivos de Teams en esta carpeta y correr:
```
python genoma.py
python prueba_rapida.py
```
