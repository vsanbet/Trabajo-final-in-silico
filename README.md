# Trabajo-final-in-silico
Reproducción y validación del análisis de perfiles proteómicos plasmáticos para la detección de sesgos en estudios de biomarcadores clínicos

Este repositorio contiene un pipeline en **R** para el preprocesamiento, control de calidad y análisis exploratorio de datos de proteómica cuantitativa generados con **MaxQuant**.  
El pipeline está diseñado para evaluar posibles **sesgos preanalíticos** en distintas fracciones sanguíneas mediante métricas de calidad y estabilidad proteómica.

---

## Estructura de los datos

El pipeline espera como entrada un directorio que contenga los archivos de salida estándar de MaxQuant, en particular:

- `proteinGroups.txt` – cuantificación a nivel de proteína
- `summary.txt` – información general de adquisición y procesamiento
- `modificationSpecificPeptides.txt` o `Oxidation (M)Sites.txt` – péptidos modificados

Cada análisis se ejecuta sobre un directorio independiente que contiene estos archivos.

---

## Análisis incluidos

El pipeline implementa los siguientes pasos:

1. **Filtrado de calidad**
   - Eliminación de contaminantes potenciales
   - Eliminación de secuencias reversas
   - Eliminación de proteínas identificadas únicamente por sitio

2. **Preprocesamiento**
   - Selección automática de intensidades LFQ (o intensidades brutas si LFQ no está disponible)
   - Transformación log₂
   - Normalización por mediana por muestra
   - Retención de proteínas cuantificadas en al menos el 70% de las muestras

3. **Análisis exploratorio**
   - Análisis de Componentes Principales (PCA)
   - Evaluación de la proporción de valores ausentes por muestra
   - Cálculo del coeficiente de variación (CV) intra-grupo
   - **Detección automática de grupos de muestras** a partir de los nombres de columna si `inspect_samples = TRUE`

4. **Control de degradación / estrés oxidativo**
   - Cuantificación de la oxidación de metionina a partir de péptidos modificados

Todos los resultados se guardan automáticamente en una carpeta `plots/` dentro del directorio analizado.

---

## Requisitos

- R ≥ 4.2
- Paquetes de R:
  - `dplyr`
  - `ggplot2`
  - `tidyverse`

---

## Uso del pipeline

1. Clona el repositorio y abre R o RStudio.
2. Carga el script principal:

```r
source("run_geyer_pipeline.R")
