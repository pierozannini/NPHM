install.packages(c("dplyr", "sf", "tmap", "tmaptools", "raster", "sp", "stringr", "magrittr"))
library(magrittr)


####################################
### Importing and preparing data ###
####################################

## Importing shapefile of habitat occurences in a gridded format

# N.B. to correctly run this script the shapefile should have EPSG: 3035
# rows (grid polygons) should be characterised by habitat occurences (columns) - "value presence (1)/absence (0)"
# Internal Note: Reminding: in "habitat_italy.shp"it is not listed habitat "1420". Anyway another shape file including habitat "1420" is ready to be used
# Internal Note: Reminding: in "habitat_italy.shp"it is listed habitat "1110" that is a marine habitat. Anyway another shape file not including habitat "1110" is ready to be used
ita_grid <- sf::st_read("habitat_italy.shp")
ita_grid <- sf::st_set_crs(ita_grid, 3035)

## Converting shapefile in a table holding coordinates of each grid centroid

tab <- sf::st_set_geometry(ita_grid, NULL)
tab <- as.data.frame(tab)
rownames(tab) <- tab[,2] # assign CELLCODEs to row names
tab <- tab[,-c(1,2)]     # remove cell IDs, this and the previous row can change according to the data employed
tab <- sf::st_centroid(ita_grid) %>% 
  sf::st_coordinates() %>% 
  cbind(tab, .)

## Calculates occurences (N) and sampling size size (n) of each habitat

N <- colSums(tab)[1:(ncol(tab)-2)]
n <- ifelse(N <= 10, N, ifelse(N >= 1200, 0.1*N ,N*0.092+10))


############################################################
### Phase 1 Stage 1 - Extract 1 quadrat for each q-block ###
############################################################

source("extract_quadrats.R")

## Extraction

index_habitat <- which(colnames(tab) == "F5110")
quadrats_habitat <- extract_quadrats(x = tab, 
                                     grd = ita_grid, 
                                     col = index_habitat, 
                                     n = n[index_habitat], 
                                     seed = 42)

## Graphical output

#svg("habitat_extracted_quadrats.svg", height = 10)
#quadrats_habitat[[3]]
#dev.off()

## Tabular output
  
do.call(rbind, sf::st_geometry(quadrats_habitat[[1]])) %>% 
  dplyr::as_tibble() %>% 
  setNames(c("lon","lat")) %>% 
  cbind(., dplyr::as_tibble(quadrats_habitat[[1]])) %>% 
  dplyr::select(3, 4, 1, 2) %>% 
  write.csv("habitat_extracted_quadrats.csv")


# N.B. "n" relative to cells is a different number than "n" relative to quadrats per habitat
# this sampling design considers each 10km x 10km quadrat is made up of 25 2km x 2km c-blocks
# the c-blocks are split in 400 cells of 100m x 100m, then each quadrat is made up of 10000 cells
# here we sample 1 cell per c-block in each quadrat (total = 25 cells per quadrat)
# different sampling design should modify the functions employed accordingly

source("quadrats_to_cells_sup.R")
source("quadrats_to_cells.R")

## Load the habitat from the tabular output and convert it to shapefile

quadrats_habitat_tab <- readr::read_csv("habitat_extracted_quadrats.csv") %>%
  sf::st_as_sf(coords = c("lon", "lat"),
               crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")

## Transform the extracted quadrats to quadrat/c-bloc/cell structure

cells_within_quadrats <- quadrats_to_cells(quadrats_habitat_tab) %>% sf::st_as_sf()
cells_within_quadrats[[1]] <- rep(1:(nrow(cells_within_quadrats)/10000), each = 10000)
cells_clustered <- cells_within_quadrats %>% 
  dplyr::left_join(as.data.frame(sf::st_set_geometry(quadrats_habitat_tab, NULL)), by = c("quadrat" = "X1"))
cells <- cells_clustered %>%
  sf::st_set_crs(3035)
cells_rep <- sf::st_transform(cells, 32632)

## Import maxent

maxent <- raster::raster("5110.tif")


## Extract maxent values
## Internal note: here we extract only the values 0<x<1. In this way we mantain only the cells with a value different from "0"  

cells_max <- exactextractr::exact_extract(maxent, cells_rep, function(values, coverage_fraction) {max(values[coverage_fraction == 1])})
ind <- which(cells_max > 0 & !is.na(cells_max))
cells_val <- cells_rep
cells_val$maxent <- cells_max
cells_val$maxent[-ind] <- 0
sf::st_write(cells_val, "cells_5110.shp")

## Extract 1 quadrat per q-block

cells_sums <- cells_val %>% 
  sf::st_set_geometry(NULL) %>% 
  dplyr::group_by(quadrat) %>% 
  dplyr::summarise(maxent_sum = sum(maxent)) %>% 
  dplyr::left_join(unique(sf::st_set_geometry(cells_val, NULL)[,c(1, 5)]), by = c("quadrat" = "quadrat")) %>% 
  dplyr::select(quadrat, cluster, maxent_sum)

extracted_1_1 <- cells_sums %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::sample_n(1, weigth = maxent)

cells_1_1 <- cells_val %>% 
  dplyr::filter(quadrat %in% extracted_1_1$quadrat)

clusters <- quadrats_habitat[[1]] %>% 
  dplyr::filter(quad %in% unique(cells_1_1$quad))

svg("maps_1_1/5110_quadrats_ph1_st1.svg")
tmap::tm_shape(ita_grid) +
  tmap::tm_borders() + 
  tmap::tm_shape(quadrats_habitat[[1]]) +
  tmap::tm_symbols(size = 0.05, col = "cluster", legend.col.show = F, palette = "Set1", n = 20, shape = 19) +
  tmap::tm_shape(clusters) +
  tmap::tm_symbols(col = "red", size = 0.2, shape = 4) +
  tmap::tm_layout(title = paste(c("N = ", N[["F5110"]], "; ", "n = ", round(n[["F5110"]])), sep = "", collapse = ""), title.size = 0.8, title.position = c("right", "top"))
dev.off()

sf::st_write(cells_1_1, "intermediate_shp/5110_1_1.shp", delete_layer = TRUE)

#########################################################################################
### Phase 1 Stage 2 - Extract n cells in each quadrat extracted in the previous stage ###
#########################################################################################

cells_1_1 <- cells_1_1 %>% 
  dplyr::filter(maxent > 0)

cells_1_1$id <- 1:nrow(cells_1_1)

extracted_1_2 <- cells_1_1 %>% 
  as.data.frame() %>% 
  dplyr::group_by(quadrat, cblock) %>% 
  dplyr::sample_n(1, weight = maxent) %>% 
  dplyr::select(id)

cells_1_2 <- cells_1_1  %>% 
  dplyr::filter(id %in% extracted_1_2$id)

sf::st_write(cells_1_2, "intermediate_shp/5110_1_2.shp", delete_layer = TRUE)

###############################################################################################
### Phase 2 Stage 2 - Sampling cells by quadrat according to PPS without replacement design ###
###############################################################################################

extracted_2 <- sf::st_set_geometry(cells_1_2, NULL) %>%
  dplyr::group_by(quadrat) %>% 
  dplyr::sample_n(if(dplyr::n() < 4) dplyr::n() else 4, weight = maxent) %>% 
  dplyr::select(id)

cells_2 <- cells_1_2 %>% 
  dplyr::filter(id %in% extracted_2$id)

sf::st_write(cells_2, "output/5110.shp", delete_layer = TRUE)

cells_ext <- cells_2 %>% 
  sf::st_centroid() %>% 
  sf::st_coordinates() %>%  
  as.matrix() %>% 
  cbind(cells_2, .)

vert <- sf::st_geometry(cells_ext) %>% sf::st_coordinates()
vert <- vert[-(seq(1:(nrow(vert)/5))*5), 1:2]
vert <- matrix(as.numeric(t(vert)), ncol = 8, byrow = T)

cells_ext <- sf::st_set_geometry(cells_ext, NULL)
cells_ext <- cbind(cells_ext, vert)
names(cells_ext)[6:ncol(cells_ext)] <- c("score", "id", "x_cent", "y_cent", "x_nw", "y_nw", "x_ne", "y_ne", "x_se", "y_se", "x_sw", "y_sw") 
cells_ext <- cells_ext %>% dplyr::arrange(as.numeric(id))
write.csv(cells_ext, "output/5110.csv")
      