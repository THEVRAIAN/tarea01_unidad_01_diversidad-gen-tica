# ============================================================
# Diversidad genética y estructura poblacional de Pedicularis dudleyi
# mediante marcadores SNP
#
# Dataset:
# Genetic data obtained from Dryad repository
#
# Objetivo:
# Caracterizar la estructura genética poblacional mediante:
#
# 1. Control y preparación de datos SNP
# 2. Análisis de Componentes Principales (PCA)
# 3. Análisis Discriminante de Componentes Principales (DAPC)
# 4. Agrupamiento genético
# 5. Diferenciación poblacional mediante FST
#

# 1. Librerías
# ============================================================
library(vcfR)
library(adegenet)
library(hierfstat)
library(StAMPP)
library(tidyverse)
library(ggplot2)
library(poppr)
library(ape)


# ============================================================
# 2. Configuración del directorio de trabajo
# ============================================================
setwd(
  "D:/2026/MAESTRIA 2026 II/DIVERSIDAD_GENÉTICA/GIT/tarea01_unidad_01_diversidad genética"
)
# ============================================================
# 3. Importación de datos SNP
# ============================================================
archivo_vcf <-
  "FINAL_cleaned_data_pedu_only_single_snp.vcf.recode.vcf"

vcf <- read.vcfR(
  archivo_vcf,
  verbose = FALSE
)

vcf

# Número de individuos

n_individuos <- ncol(vcf@gt)-1

# Número de SNP

n_snp <- nrow(vcf@fix)

cat(
  "Individuos:",
  n_individuos,
  "\nSNP:",
  n_snp
)

# ============================================================
# 4. Conversión del VCF a formato genético
# ============================================================

# Conversión del archivo VCF a objeto genlight

genlight_obj <-
  vcfR2genlight(vcf)

genlight_obj

# ============================================================
# 5. Importación de información poblacional
# ============================================================

metadata <-
  read.csv(
    "metadata_pedicularis.csv",
    stringsAsFactors = FALSE
  )

head(metadata)

# ============================================================
# 6. Integración de información poblacional
# ============================================================
# Verificación de coincidencia entre nombres

sample_geneticos <-
  indNames(genlight_obj)

sample_metadata <-
  metadata$Sample

# Ordenamiento del metadata según el VCF

metadata <-
  metadata[
    match(sample_geneticos,
          metadata$Sample),
  ]


# Asignación de población

pop(genlight_obj) <-
  metadata$Population

table(pop(genlight_obj))

# ============================================================
# 7. Evaluación inicial del dataset SNP
# ============================================================

# Porcentaje de datos faltantes

missing_data <-
  mean(is.na(as.matrix(genlight_obj)))

missing_data

# Frecuencia de individuos por población

table(pop(genlight_obj))

# ============================================================
# 8. Análisis de Componentes Principales (PCA)
# ============================================================

# El PCA fue realizado utilizando la matriz genómica
# derivada de los marcadores SNP.

pca_genetico <-
  glPca(
    genlight_obj,
    nf = 10
  )


# Varianza explicada

varianza_pca <-
  pca_genetico$eig /
  sum(pca_genetico$eig)


varianza_pca[1:5]

# ============================================================
# 9. Preparación de resultados PCA
# ============================================================


pca_resultados <-
  data.frame(
    Individuo = indNames(genlight_obj),
    PC1 = pca_genetico$scores[,1],
    PC2 = pca_genetico$scores[,2],
    Poblacion = pop(genlight_obj)
  )


head(pca_resultados)

# ============================================================
# 10. Visualización PCA
# ============================================================

grafico_pca <-
  ggplot(
    pca_resultados,
    aes(
      x = PC1,
      y = PC2,
      color = Poblacion
    )
  )+
  
  geom_point(
    size = 3
  )+
  
  theme_classic()+
  
  labs(
    title =
      "Estructura genética de Pedicularis dudleyi mediante PCA",
    
    x =
      paste0(
        "PC1 (",
        round(varianza_pca[1]*100,2),
        "%)"
      ),
    
    y =
      paste0(
        "PC2 (",
        round(varianza_pca[2]*100,2),
        "%)"
      )
    
  )


grafico_pca

# ============================================================
# 11. Análisis Discriminante de Componentes Principales (DAPC)//
# ============================================================
# El DAPC fue realizado utilizando las poblaciones definidas
# según el sitio de colecta para evaluar la diferenciación
# genética entre grupos poblacionales.


dapc_poblacion <-
  dapc(
    genlight_obj,
    pop(genlight_obj),
    n.pca = 20,
    n.da = 2
  )


scatter(
  dapc_poblacion,
  scree.da = FALSE,
  scree.pca = FALSE,
  legend = TRUE,
  pch = 19,
  cex = 1.5,
  main =
    "Estructura genética poblacional mediante DAPC"
)
# ============================================================
# 12. Agrupamiento jerárquico
# ============================================================

distancia_genetica <-
  dist(
    genlight_obj
  )



cluster_genetico <-
  hclust(
    distancia_genetica,
    method="ward.D2"
  )


plot(
  cluster_genetico,
  main =
    "Dendrograma basado en distancia genética",
  xlab="Individuos",
  sub=""
)

# ============================================================
# 13. Diferenciación genética mediante FST
# ============================================================


# Conversión del objeto genlight para análisis de diferenciación
# genética entre poblaciones mediante SNP


fst_data <-
  stamppConvert(
    genlight_obj,
    "genlight"
  )


fst_resultado <-
  stamppFst(
    fst_data,
    nboots = 100,
    percent = 95,
    nclusters = 1
  )


fst_resultado$Fsts
fst_data <-
  stamppConvert(
    genlight_obj,
    "genlight"
  )
# FST entre pares de poblaciones
fst_resultado <-
  stamppFst(
    fst_data,
    nboots = 100,
    percent = 95,
    nclusters = 1
  )
# Matriz de diferenciación genética

fst_resultado$Fsts

# ============================================================
# 14. Visualización de diferenciación genética
# ============================================================
fst_heatmap <- fst_resultado$Fsts

fst_heatmap[is.na(fst_heatmap)] <- 0


pheatmap(
  fst_heatmap,
  display_numbers = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main =
    "Matriz de diferenciación genética FST"
)

# ============================================================
# ============================================================

# DEBIDO A QUE PARA HACER EL ANÁLISIS EN NUESTRA BASE DE DATOS SE TUVO QUE FILTRAR NA Y SNP MÓRFICOS ENTRE OTROS.
# Convertir a matriz SNP
snp_matrix <- as.matrix(genlight_obj)
dim(snp_matrix)
head(snp_matrix[,1:10])
missing_snp <- colMeans(is.na(snp_matrix))

summary(missing_snp)

variacion_snp <- apply(
  snp_matrix,
  2,
  function(x) length(unique(na.omit(x)))
)

table(variacion_snp)
variacion_snp <- apply(
  snp_matrix,
  2,
  function(x) length(unique(na.omit(x)))
)

table(variacion_snp)

# ============================================================
# Filtrado de SNP informativos
# ============================================================


snp_matrix <- as.matrix(genlight_obj)


# ------------------------------------------------------------
# 1. Eliminación de SNP con alta proporción de datos faltantes
# ------------------------------------------------------------

missing_snp <- colMeans(is.na(snp_matrix))


snp_filtrado <- 
  snp_matrix[, missing_snp < 0.20]


# ------------------------------------------------------------
# 2. Eliminación de SNP monomórficos
# ------------------------------------------------------------


snp_variacion <-
  apply(
    snp_filtrado,
    2,
    function(x) length(unique(na.omit(x))) > 1
  )


snp_filtrado <-
  snp_filtrado[, snp_variacion]


# Dimensiones finales

dim(snp_filtrado)

genlight_filtrado <-
  new(
    "genlight",
    snp_filtrado
  )


# asignar nombres de individuos

indNames(genlight_filtrado) <-
  rownames(snp_filtrado)


# asignar población

pop(genlight_filtrado) <-
  pop(genlight_obj)


genlight_filtrado

class(genlight_obj)
class(snp_filtrado)

fst_data <- stamppConvert(
  genlight_obj,
  "genlight"
)
fst_resultado <- stamppFst(
  fst_data,
  nboots = 100,
  percent = 95,
  nclusters = 1
)
fst_resultado$Fsts

table(
  pop(genlight_obj)
)
table(
  clusters$grp,
  pop(genlight_obj)
)
