# E1. Actividad Integradora 1 - Avance

Analisis del genoma de SARS-CoV-2 (Wuhan 2019 vs Texas 2020) con algoritmos de strings implementados desde cero.

## Algoritmos utilizados
| Punto | Algoritmo | Complejidad |
|---|---|---|
| 1. Indices de cada gen en el genoma | KMP (`kmp`, `tabla_fallo` en `genoma.py`) | O(n + m) |
| 2. Palindromo mas largo por gen | Manacher (`palindromo_mas_largo`) | O(n) |

`prueba_rapida.py` valida ambos contra fuerza bruta (2000 casos aleatorios): **OK**.

## Que ya esta completo
- Lectura automatica de los archivos (FASTA o texto plano).
- **Punto 1**: busqueda de los genes M, S y ORF1AB en el genoma (indices y primeros 12 caracteres).
- **Punto 2**: palindromo mas largo de cada gen, longitud mostrada y guardada en `palindromos.txt`.
- Pruebas de correctitud de KMP y Manacher.

## Que falta
- **Punto 3**: ubicar en el genoma cada proteina de `seq-proteins.txt` (tabla de codones, traduccion en los 3 marcos de lectura, primeros 4 aminoacidos y codones).
- **Punto 4**: comparar Wuhan vs Texas (indices con diferencias, codones/aminoacidos afectados, subcadena comun mas larga).
- Capturas de pantalla de resultados con los archivos reales.
- Reporte (PDF) y video.

## Ejecucion
Colocar los archivos de Teams en esta carpeta y correr:
```
python genoma.py
python prueba_rapida.py
```

## Evidencia
Capturas en la carpeta `capturas/`.
