# Proyecto de Sketches

Este proyecto implementa y evalúa diferentes algoritmos de sketches aplicados a datos de genomas.

## Ejecución

En el directorio raíz se encuentran los programas principales que permiten ejecutar cada método:

- **Count Sketch**  
  ```bash
  g++ calcular_cs.cpp -o calcular_cs
  ./calcular_cs
  ```
  Guarda en la carpeta `CSV/` un archivo `.csv` con los *heavy hitters* y sus datos.

- **Tower Sketch**  
  ```bash
  g++ calcular_ts.cpp -o calcular_ts
  ./calcular_ts
  ```
  Guarda en la carpeta `CSV/` un archivo `.csv` con los *heavy hitters* y sus datos.

- **Ground Truth**  
  ```bash
  g++ ground_truth.cpp -o ground_truth
  ./ground_truth
  ```
  Guarda en la carpeta `CSV/` un archivo `.csv` con los *heavy hitters* y sus datos.

- **Calibracion Sketches**
  ```bash
  g++ calibracion_sketchs.cpp -o calibracion_sketchs
  ./calibracion_sketchs
  ```
  Guarda en la carpeta `CSV/` un archivo `.csv` con los resultados de la calibración.

## Estructura de carpetas

- **`CSV/`**  
  Contiene los resultados en formato `.csv` generados por cada ejecución con los *heavy hitters*.

- **`sketchs/`**  
  Contiene las implementaciones de los algoritmos de sketches en formato header:
  - `countsketch.hpp`
  - `towersketch.hpp`
  - `murmurhash32.hpp`

- **`utils/`**  
  Contiene herramientas auxiliares:
  - `MetricasEvaluacion.hpp`: métricas de evaluación.  
  - `LectorGenomas.hpp`: extracción de genomas desde los archivos de entrada.

- **`results_calibracion/`**  
  Contiene los resultados de las calibraciones y sus respectivos gráficos.

## Notas

- Todos los programas deben compilarse y ejecutarse desde el directorio raíz.
- Los resultados siempre se guardan automáticamente en la carpeta `CSV/`.
