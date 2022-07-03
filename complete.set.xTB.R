#####################################################
###############      FULL SCRIPT      ###############
#####################################################

# Requires the folllowing packages to be installed:

list.of.packages <- c("data.table", "stringr", "tibble", "dplyr", "knitr",
                      "tools", "pracma", "igraph", "filesstrings", "diversitree", "plyr",
                      "caret", "matlib", "ggplot2", "ggrepel", "extrafont"
                      )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)


#### Aliases  - Not Necessary ####

l <- list.files()
cd <- setwd
### Support functions ###

xyz_file_generator <- function(dir) {
  options(scipen = 999)
  setwd(dir)
  #unlink(list.files(pattern = '*.xyz'))
  molecules <- list.files(full.names = F, recursive = F, pattern = "xyz_")
  for (molecule in molecules) {
    xyz <- data.frame(read.csv(molecule, header = F))
    num.atoms <- nrow(xyz[-1, ])
    m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
    m[1, 1] <- num.atoms
    m[is.na(m)] <- ""
    names(m) <- names(xyz)
    xyz <- rbind(m, xyz[-1, ])
    new_xyz <- knitr::kable(xyz, format = "simple", row.names = F, col.names = NULL, align = "r")
    new_xyz <- new_xyz[-1]
    new_xyz <- new_xyz[-length(new_xyz)]
    write(new_xyz, paste(stringr::str_remove(tools::file_path_sans_ext(molecule), 'xyz_'), ".xyz", sep = ""))
  }
  setwd("..")
}

xyz_file_generator_library <- function(dir) {
  dir.create("readable_xyz", showWarnings = F)
  setwd(dir)
  molecules <- list.files(full.names = F, recursive = F)
  for (molecule in molecules) {
    xyz_file_generator(molecule)
    setwd(molecule)
    filesstrings::file.move(list.files(pattern = ".xyz"), "../../readable_xyz", overwrite = T)
    setwd('..')
  }
  setwd("..")
}

xyz_file_generator_library_batch <- function() {
  dirs <- list.files()
  for (dir in dirs) {
    setwd(dir)
    xyz_file_generator_library(list.files(pattern = 'extracted'))
    setwd('..')
  }
}

name_changer <- function(from, to) {
  molecules <- list.files(recursive = F, full.names = F)
  for (molecule in molecules) {
    file.rename(molecule, stringr::str_replace(molecule, from, to))
  }
}

angle <- function(x,y){
  dot.prod <- x %*% y 
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}
# change coordinates system

coor.trans <- function(molecule, coor.atoms) {
  atoms <- strsplit(coor.atoms, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  xyz <- data.table::fread(molecule, header = F)
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  if (length(numeric.atoms) == 4) {
    new_origin <- (xyz[numeric.atoms[[1]], 2:4] + xyz[numeric.atoms[[2]], 2:4])/2
    new_y <- as.numeric((xyz[numeric.atoms[[3]], 2:4] - new_origin) /
                          mag(xyz[numeric.atoms[[3]], 2:4] - new_origin))
    coplane <- as.numeric((xyz[numeric.atoms[[4]], 2:4] - new_origin) /
                            mag(xyz[numeric.atoms[[4]], 2:4] - new_origin))
  } else {
    new_origin <- xyz[numeric.atoms[[1]], 2:4]
    new_y <- as.numeric((xyz[numeric.atoms[[2]], 2:4] - new_origin) /
                          mag(xyz[numeric.atoms[[2]], 2:4] - new_origin))
    coplane <- as.numeric((xyz[numeric.atoms[[3]], 2:4] - new_origin) /
                            mag(xyz[numeric.atoms[[3]], 2:4] - new_origin))
  }
  cross_y_coplane <- pracma::cross(coplane, new_y)
  coef_mat <- aperm(array(c(
    new_y,
    coplane,
    cross_y_coplane
  ),
  dim = c(3, 3)
  ))
  angle_new.y_coplane <- angle(coplane, new_y)
  x_ang_new.y <- pi / 2
  cop_ang_x <- angle_new.y_coplane - x_ang_new.y
  result_vec <- c(0, cos(cop_ang_x), 0)
  new_x <- solve(coef_mat, result_vec)
  new_z <- pracma::cross(new_x, new_y)
  new_basis <- aperm(array(c(new_x, new_y, new_z), dim = c(3, 3)))
  new_origin <- xyz[numeric.atoms[[1]], 2:4]
  new_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
  transformed_coordinates <- matrix(nrow = dim(xyz)[[1]], ncol = 3)
  for (i in 1:dim(xyz)[[1]]) {
    new_coordinates[i, ] <- as.numeric(xyz[i, 2:4] - new_origin)
    transformed_coordinates[i, ] <- aperm(new_basis %*% new_coordinates[i, ])
  }
  transformed_coordinates <- round(transformed_coordinates, 4)
  elements <- xyz[, 1]
  transformed_coordinates <- cbind(elements, transformed_coordinates)
  num.atoms <- nrow(transformed_coordinates)
  m <- as.data.frame(matrix(NA, ncol = 4, nrow = 2))
  m[1, 1] <- num.atoms
  m[is.na(m)] <- ""
  names(m) <- names(transformed_coordinates)
  transformed_coordinates <- rbind(m, transformed_coordinates)
  transformed_coordinates[transformed_coordinates == "0"] <- "0.0"
  new_xyz <- knitr::kable(transformed_coordinates,
                          format = "simple", row.names = F, col.names = NULL
  )
  new_xyz <- new_xyz[-1]
  new_xyz <- new_xyz[-length(new_xyz)]
  write(new_xyz, paste(tools::file_path_sans_ext(molecule), "_tc", ".xyz", sep = ""))
}

coor.trans.dir <- function(dir, coor.atoms){
  setwd(dir)
  molecules <- list.files(pattern = '.xyz', full.names = F, recursive = F)
  for (molecule in molecules) {
    coor.trans(molecule, coor.atoms)
  }
  setwd('..')
}

# cross validation and model selection

cv <- function(formula, data, out.col, folds) {
  predictions <- list()
  if (folds == nrow(data)) {
    split.assign <- sample(1:folds, nrow(data), replace = F)
  } else {
    split.assign <- caret::createFolds(1:dim(data)[1], folds, list = F)
  }
  cvdat <- cbind(data, split.assign)
  colnames(cvdat)[dim(data)[2] + 1] <- "assign"
  for (i in 1:folds) {
    train <- cvdat[cvdat$assign != i, ]
    test <- cvdat[cvdat$assign == i, ]
    model <- lm(formula, data = train)
    predictions[[i]] <- data.frame(predict(model, newdata = test))
  }
  pred <- dplyr::bind_rows(predictions)
  pred <- pred[match(rownames(cvdat), rownames(pred)), ]
  sq.dis.cv <- data.frame(abs(pred - data[, out.col]))
  MAE.3cv <- mean(sq.dis.cv$abs.pred...data...out.col..)
  q2.3cv <- caret::R2(pred, data[, out.col])
  return(list(MAE.3cv, q2.3cv))
}

measure_cv <- function(formula, data, out.col, folds, iterations) {
  mae.list <- list()
  q2.list <- list()
  for (i in 1:iterations) {
    tool <- cv(formula, data, out.col, folds)
    mae.list[[i]] <- tool[[1]]
    q2.list[[i]] <- tool[[2]]
  }
  MAE.validation <- Reduce(`+`, mae.list) / iterations
  q2.validation <- Reduce(`+`, q2.list) / iterations
  return(list(MAE.validation, q2.validation))
}

sub_model <- function(data, out.col = dim(data)[2],
                      min = 3, max = floor(dim(data)[1] / 5),
                      folds = nrow(data), iterations = 1,
                      cutoff = 0.85, cross.terms = F) {
  output <- stringr::str_c("`", names(data[out.col]), "`")
  vars <- names(data[, -out.col])
  for (i in 1:length(vars)) {
    vars[i] <- stringr::str_c("`", vars[i], "`")
  }
  comb.list <- list()
  ols.list <- list()
  q2.list <- list()
  mae.list <- list()
  if (cross.terms == T) {
    max <- max - 1
  }
  for (i in min:max) {
    comb.list[[i]] <- data.frame(aperm(combn(vars, i)), stringsAsFactors = F)
    comb.list[[i]][, dim(comb.list[[i]])[2] + 1] <- do.call(
      paste,
      c(comb.list[[i]][names(comb.list[[i]])],
        sep = " + "
      )
    )
    names(comb.list[[i]])[dim(comb.list[[i]])[2]] <- "formula"
    for (co in names(comb.list[[i]])[1:length(names(comb.list[[i]])) - 1]) comb.list[[i]][co] <- NULL
  }
  comb.list <- plyr::compact(comb.list)
  forms <- do.call(rbind, comb.list)
  if (cross.terms == T) {
    cross.list <- list()
    cross.list[[1]] <- data.frame(aperm(combn(vars, 2)), stringsAsFactors = F)
    cross.list[[2]] <- do.call(paste, c(cross.list[[1]][names(cross.list[[1]])], sep = " : "))
    repli_single <- do.call(rbind, replicate(length(cross.list[[2]]), forms, simplify = F))
    repli_single <- repli_single[order(repli_single), ]
    repli_cross <- do.call(c, replicate(nrow(forms), cross.list[[2]], simplify = F))
    new_comb <- list()
    for (i in 1:length(repli_cross)) {
      new_comb[[i]] <- paste(repli_single[[i]], repli_cross[[i]], sep = " + ")
    }
    orig.forms <- forms
    forms <- data.frame(do.call(rbind, plyr::compact(new_comb)), stringsAsFactors = F)
    names(forms) <- "formula"
    forms <- data.frame(rbind(orig.forms, forms))
    colnames(forms) <- "formula"
  }
  for (i in 1:dim(forms)[1]) {
    forms$formula[i] <- stringr::str_c(output, " ~ ", forms$formula[i])
    ols.list[[i]] <- summary(lm(forms[i, ], data = data))$r.squared
  }
  forms[, 2] <- do.call(rbind, ols.list)
  names(forms)[2] <- "R.sq"
  forms.cut <- forms[forms$R.sq > cutoff, ]
  while (nrow(forms.cut) == 0) {
    cutoff <- cutoff - 0.1
    forms.cut <- forms[forms$R.sq > cutoff, ]
  }
  for (i in 1:dim(forms.cut)[1]) {
    stts <- measure_cv(forms.cut[i, 1], data, out.col, folds, iterations)
    q2.list[[i]] <- stts[2]
    mae.list[[i]] <- stts[1]
  }
  forms.cut[, 3] <- data.table::transpose(do.call(rbind, q2.list))
  forms.cut[, 4] <- data.table::transpose(do.call(rbind, mae.list))
  names(forms.cut)[3:4] <- c("Q.sq", "MAE")
  return(head(dplyr::arrange(forms.cut, desc(forms.cut$Q.sq)), 2))
}

### Constructive functions ###

mol.info <- function(info_filename, vib_num_filename) {
  info <- data.frame(data.table::fread(info_filename,
                                       sep = " ", header = F, fill = T
  ))
  leave.out <- c(
    "Frequencies",
    "IR",
    "Inten",
    "--",
    "A"
  )
  for (i in leave.out) {
    info[info == i] <- NA
  }
  info.nonarow <- info[rowSums(is.na(info)) != ncol(info), ]
  row.names(info.nonarow) <- 1:dim(info.nonarow)[1]
  seq.1 <- seq(1, dim(info.nonarow)[1], 3)
  seq.2 <- seq(2, dim(info.nonarow)[1], 3)
  seq.3 <- seq(3, dim(info.nonarow)[1], 3)
  info.nonarow[seq.2, 1:3] <- info.nonarow[seq.2, 3:5]
  info.nonarow[seq.3, 1:3] <- info.nonarow[seq.3, 4:6]
  info.clean <- info.nonarow[, 1:3]
  block.list <- list()
  for (i in seq.1) {
    a <- info.clean[i, ]
    b <- info.clean[i + 1, ]
    c <- info.clean[i + 2, ]
    d <- rbind(a, b, c)
    block.list[[i]] <- assign(
      paste("block", i, sep = "."),
      data.frame(d)
    )
  }
  block.list.compact <- plyr::compact(block.list)
  flat.info <- do.call(cbind, block.list.compact)
  names(flat.info) <- flat.info[1, ]
  flat.info <- flat.info[-1, ]
  row.names(flat.info) <- c("Frequency [1/cm]", "IR intensity")
  vib.numbered <- data.table::fread(vib_num_filename)
  vib <- vib.numbered[, -c(1:2)]
  names(vib) <- c("x", "y", "z", "x", "y", "z", "x", "y", "z")
  iterations <- dim(vib)[1]
  datalist <- list()
  for (i in 1:iterations) {
    a <- vib[i, 1:3]
    b <- vib[i, 4:6]
    c <- vib[i, 7:9]
    d <- rbind(a, b, c)
    datalist[[i]] <- assign(
      paste("out", i, sep = "."),
      data.frame(d, row.names = c(1, 2, 3))
    )
  }
  output <- do.call(rbind, datalist)
  row.names(output) <- seq(1, iterations * 3, 1)
  for (i in 1:dim(output)[1]) {
    output[i, 4] <- sqrt(output[i, 1]^2 +
                           output[i, 2]^2 +
                           output[i, 3]^2)
  }
  names(output)[[4]] <- "magnitude"
  output[, 5] <- data.table::transpose(flat.info[1, ])
  names(output)[[5]] <- "frequency"
  output <- transform(output, frequency = as.numeric(frequency))
  output_finger <- output[(output$frequency > 1500), ]
  max.mag.vib <- which.max(output_finger[, 4])
  max.mag.vib <- row.names(output_finger[max.mag.vib, ])
  your.vib <- data.frame(flat.info[1:2, max.mag.vib])
  row.names(your.vib) <- c("Frequency [1/cm]", "IR intensity")
  colnames(your.vib) <- sub("\\.csv$", "", vib_num_filename)
  your.vib
  return(flat.info)
  return(your.vib)
}
## NO NEED FOR MOL.INFO - USE VIBESPECTRUM

### 2 Sterimol (paton.info & and paton.df - depricated -- in-house steRimol and steRimol.df)

steRimol <- function(mol.dir, coordinates, radii = 'CPK', only.sub = T) {
  tryCatch(
    expr = {
      origin <- as.numeric(unlist(strsplit(coordinates, " "))[[1]])
      direction <- as.numeric(unlist(strsplit(coordinates, " "))[[2]])
      bondi <- data.frame(matrix(nrow = 12, ncol = 2))
      bondi[, 1] <- c("H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I")
      bondi[, 2] <- c(1.09, 1.92, 1.70, 1.55, 1.52, 1.47, 2.10, 1.80, 1.80, 1.75, 1.83, 1.98)
      colnames(bondi) <- c("Atom", "Radius")
      Paton.Atypes <- data.frame(c("C", "C2", "C3", "C4", "C5/N5", "C6/N6", "C7", "C8", "H", "N", "C66", "N4", "O", "O2", "P", "S", 'S.O', "S1", "F", "Cl", "S4", "Br", "I"), 
                                 c(1.50,1.60,1.60,1.50,1.70,1.70,1.70,1.50,1.00,1.50,1.70,1.45,1.35,1.35,1.40,1.70,1.70,1.00,1.35,1.80,1.40,1.95,2.15))
      colnames(Paton.Atypes) <- c("Atom", "Radius")
      Paton.Atypes <- rbind(Paton.Atypes, bondi[!bondi$Atom %in% Paton.Atypes$`Atom`,])
      xyz_file_generator(mol.dir)
      setwd(mol.dir)
      bonds <- data.table::fread(list.files(pattern = 'wbo'))
      bonds <- bonds[bonds$V3 > 0.8, ]
      outsiders <- vector()
      for (i in bonds$V1) {
        if ((i %in% bonds$V1) && !(i %in% bonds$V2) && (i != 1)) {
          outsiders <- append(outsiders, i)
          bonds <- bonds[-which(bonds$V1 == i),]
        }
      }
      bonds <- bonds[, 1:2]
      bonds <- data.frame(sapply(bonds, as.numeric))
      if (is.na(bonds[bonds$V1 == direction, ][1, 2])) {
        coordinates <- paste(coordinates, as.character(bonds[bonds$V2 == direction, ][1, 1]), sep = " ")
      } else {
        coordinates <- paste(coordinates, as.character(bonds[bonds$V1 == direction, ][1, 2]), sep = " ")
      }
      if (as.numeric(unlist(strsplit(coordinates, " "))[[3]]) == origin) {
        coordinates <- as.numeric(unlist(strsplit(coordinates, " "))[1:2])
        coordinates <- paste(as.character(coordinates[1]), as.character(coordinates[2]), sep = ' ')
        if (is.na(bonds[bonds$V1 == direction, ][1, 2])) {
          coordinates <- paste(coordinates, as.character(bonds[bonds$V2 == direction, ][2, 1]), sep = " ")
        } else {
          coordinates <- paste(coordinates, as.character(bonds[bonds$V1 == direction, ][2, 2]), sep = " ")
        }
      }
      if (unlist(strsplit(coordinates, " "))[3] == 'NA') {
        coordinates <- as.numeric(unlist(strsplit(coordinates, " "))[1:2])
        coordinates <- paste(as.character(coordinates[1]), as.character(coordinates[2]), sep = ' ')
        remove.direction <- bonds[bonds$V2 != direction, ]
        if (!is.na(remove.direction[remove.direction$V1 == origin, ][1, 2])) {
          coordinates <- paste(coordinates, as.character(remove.direction[remove.direction$V1 == origin, ][1, 2]), sep = " ")
        } else {
          coordinates <- paste(coordinates, as.character(remove.direction[remove.direction$V1 == origin, ][1, 2]), sep = " ")
        }
      }
      coor.trans(list.files(pattern = ".xyz"), coordinates)
      mag <- function(vector) {
        sqrt(vector[[1]]^2 + vector[[2]]^2)
      }
      substi <- data.table::fread(list.files(pattern = "tc.xyz"))
      substi <- tibble::rownames_to_column(substi)
      if (only.sub == T) {
        G <- igraph::graph.data.frame(bonds, directed = T)
        C <- igraph::all_simple_paths(G, as.character(origin), mode = "all")
        rlev <- vector()
        for (i in 1:length(C)) {
          if (length(C[[i]]) > 2) {
            rlev[[i]] <- stringr::str_extract_all(as.character(C[i]), "`[0-9]{1,3}`")
          }
          if (length(C[[i]]) == 2) {
            two_atoms <- stringr::str_extract_all(as.character(C[i]), "[0-9]{1,3}")[[1]][1:4]
            if (as.character(origin) %in% two_atoms && as.character(direction) %in% two_atoms) {
              rlev[[i]] <- stringr::str_extract_all(as.character(C[i]), "`[0-9]{1,3}`")
            }
          }
        }
        rlev <- rlev[grep(paste("`", direction, "`", sep = ""), rlev)]
        rlev_atoms <- as.numeric(stringr::str_extract_all(unique(unlist(rlev)), "[0-9]{1,3}"))
        substi <- substi[substi$rowname %in% rlev_atoms]
      }
      from_source <- dplyr::intersect(
        bonds[bonds$V1 %in% substi$rowname, ],
        bonds[bonds$V2 %in% substi$rowname, ]
      )
      nms <- data.frame(matrix(nrow = nrow(from_source), ncol = 1))
      colnames(nms) <- "atom"
      for (i in 1:dim(from_source)[1]) {
        nms[i, ] <- substi$V1[substi$rowname == from_source$V1[i]]
      }
      nms.2 <- data.frame(matrix(nrow = nrow(from_source), ncol = 1))
      colnames(nms.2) <- "atom.2"
      for (i in 1:dim(from_source)[1]) {
        nms.2[i, ] <- substi$V1[substi$rowname == from_source$V2[i]]
      }
      from_source <- cbind(nms, from_source, nms.2)
      for (i in as.numeric(unique(as.vector(t(from_source[, -c(1, 4)]))))) {
        if (i %in% from_source$V1 && from_source$atom[from_source$V1 == i] == "H") {
          from_source <- from_source[from_source$V1 != i, ]
        }
        if (i %in% from_source$V2 && from_source$atom.2[from_source$V2 == i] == "H" && from_source$atom[from_source$V2 == i] %in% c("O", "N", "F")) {
          bonded_OH <- from_source$V1[from_source$V2 == i]
          for (j in bonded_OH) {
            bond.length <- sqrt(sum((substi[substi$rowname == i, 3:5] - substi[substi$rowname == j, 3:5])^2))
            if (bond.length > 1.6) {
              from_source <- from_source[from_source$V2 != i, ]
            }
          }
        }
      }
      bonds.clean <- from_source[, 2:3]
      substi <- substi[substi$rowname %in% bonds.clean$V1 | substi$rowname %in% bonds.clean$V2, ]
      substi <- plyr::mutate(substi, magnitude = mag(substi[, c(3, 5)]))
      substi <- plyr::mutate(substi, Radius = rep(0, nrow(substi)))
      atypes <- data.table::fread(list.files(pattern = 'atypes'), header = F)
      atypes <- atypes[row.names(atypes) %in% substi$rowname]
      trans.table <- data.frame(
        c("O.2","O.3","C.1","C.2",'C.3','C.ar','N.2','N.1','N.3','N.ar','N.am','N.pl3','N.4','S.2','S.3','S.O2'), 
        c("O2","O",'C3','C2','C','C6/N6','C6/N6','N','C6/N6','C6/N6','C6/N6','C6/N6','N',"S",'S4','S1'))
      names(trans.table) <- c('mol2', 'Verloop')
      for (i in 1:nrow(atypes)) {
        if (atypes[i, ] %in% trans.table$mol2) {
          atypes[i, ] <- trans.table$Verloop[which(trans.table$mol2 %in% atypes[i, ])]
        }
      }
      substi <- plyr::mutate(substi, atypes = atypes$V1)
      for (i in 1:dim(substi)[1]) {
        if (radii == 'bondi') {
          substi$Radius[i] <- bondi$Radius[bondi$Atom == substi$V1[i]]
        }
        if (radii == 'CPK') {
          substi$Radius[i] <- Paton.Atypes$Radius[Paton.Atypes$Atom == substi$atypes[i]]
        }
      }
      substi <- plyr::mutate(substi, Bs = substi$magnitude + substi$Radius)
      substi <- plyr::mutate(substi, L = substi$V3 + substi$Radius)
      substi$L[1] <- 0
      scans <- 90 / 5
      deg.list <- seq(scans, 90, scans)
      plane <- substi[, c(3, 5)]
      b1s <- vector()
      b1s.loc <- vector()
      for (i in deg.list) {
        rot.mat <- matrix(ncol = 2, nrow = 2)
        co.deg <- cos(i * (pi / 180))
        si.deg <- sin(i * (pi / 180))
        rot.mat[1, ] <- c(co.deg, -1 * si.deg)
        rot.mat[2, ] <- c(si.deg, co.deg)
        ind <- NULL
        tc <- matrix(nrow = dim(plane)[[1]], ncol = 2)
        for (i in 1:dim(plane)[[1]]) {
          tc[i, ] <- aperm(rot.mat %*% as.numeric(plane[i, ]))
        }
        tc <- round(tc, 3)
        avs <- abs(c(
          max(tc[, 1]),
          min(tc[, 1]),
          max(tc[, 2]),
          min(tc[, 2])
        ))
        if (min(avs) == 0) {
          if (any(which(avs == min(avs)) %in% c(1, 2))) {
            tc <- round(tc, 1)
            B1 <- max(substi$Radius[which(tc[, 1] == 0, arr.ind = T)])
            B1.loc <- substi$L[which.max(substi$Radius[which(tc[, 1] == 0, arr.ind = T)])]
            b1s.loc <- append(b1s.loc, B1.loc)
            b1s <- append(b1s, B1)
          }
          if (any(which(avs == min(avs)) %in% c(3, 4))) {
            tc <- round(tc, 1)
            B1 <- max(substi$Radius[which(tc[, 2] == 0, arr.ind = T)])
            B1.loc <- substi$L[which.max(substi$Radius[which(tc[, 2] == 0, arr.ind = T)])]
            b1s.loc <- append(b1s.loc, B1.loc)
            b1s <- append(b1s, B1)
          }
        } else {
          if (any(which(avs == min(avs)) %in% c(1, 2))) {
            ind <- (which(abs(tc[, 1]) == min(avs), arr.ind = T))[1]
            against <- vector()
            against.loc <- vector()
            if (any(tc[ind, 1] < 0)) {
              new.ind <- which(tc[, 1] == min(tc[, 1]), arr.ind = T)
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 1], tc[new.ind[1], 1], tc[new.ind[1], 1] + 1))
              tc[,1] <- -tc[,1]
            } else {
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 1], tc[ind[1], 1] - 1, tc[ind[1], 1]))
            }
            for (i in 2:dim(tc)[1]) {
              if (tc[i, 3] == T) {
                against <- append(against, tc[i, 1] + substi$Radius[i])
                against.loc <- append(against.loc, substi$L[i])
              }
              B1 <- ifelse(length(against) > 0, max(against), abs(tc[ind[1], 1]) + substi$Radius[ind[1]])
              B1.loc <- ifelse(length(against) > 0, against.loc[which.max(against)], substi$L[ind[1]])
            }
            b1s <- append(b1s, B1)
            b1s.loc <- append(b1s.loc, B1.loc)
          }
          if (any(which(avs == min(avs)) %in% c(3, 4))) {
            ind <- (which(abs(tc[, 2]) == min(avs), arr.ind = T))[1]
            against <- vector()
            against.loc <- vector()
            if (any(tc[ind, 2] < 0)) {
              new.ind <- which(tc[, 2] == min(tc[, 2]), arr.ind = T)
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 2], tc[new.ind[1], 2], tc[new.ind[1], 2] + 1))
              tc[,2] <- -tc[,2]
            } else {
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 2], tc[ind[1], 2] - 1, tc[ind[1], 2]))
            }
            for (i in 2:dim(tc)[1]) {
              if (tc[i, 3] == T) {
                against <- append(against, tc[i, 2] + substi$Radius[i])
                against.loc <- append(against.loc, substi$L[i])
              }
              B1 <- ifelse(length(against) > 0, max(against), abs(tc[ind[1], 2]) + substi$Radius[ind[1]])
              B1.loc <- ifelse(length(against) > 0, against.loc[which.max(against)], substi$L[ind[1]])
            }
            b1s <- append(b1s, B1)
            b1s.loc <- append(b1s.loc, B1.loc)
          }
        }
      }
      
      back.ang <- deg.list[which(b1s == min(b1s))][1] - scans
      front.ang <- deg.list[which(b1s == min(b1s))][1] + scans
      thin <- seq(back.ang, front.ang, 1)
      
      for (i in thin) {
        rot.mat <- matrix(ncol = 2, nrow = 2)
        co.deg <- cos(i * (pi / 180))
        si.deg <- sin(i * (pi / 180))
        rot.mat[1, ] <- c(co.deg, -1 * si.deg)
        rot.mat[2, ] <- c(si.deg, co.deg)
        ind <- NULL
        tc <- matrix(nrow = dim(plane)[[1]], ncol = 2)
        for (i in 1:dim(plane)[[1]]) {
          tc[i, ] <- aperm(rot.mat %*% as.numeric(plane[i, ]))
        }
        tc <- round(tc, 3)
        avs <- abs(c(
          max(tc[, 1]),
          min(tc[, 1]),
          max(tc[, 2]),
          min(tc[, 2])
        ))
        if (min(avs) == 0) {
          if (any(which(avs == min(avs)) %in% c(1, 2))) {
            tc <- round(tc, 1)
            B1 <- max(substi$Radius[which(tc[, 1] == 0, arr.ind = T)])
            B1.loc <- substi$L[which.max(substi$Radius[which(tc[, 1] == 0, arr.ind = T)])]
            b1s.loc <- append(b1s.loc, B1.loc)
            b1s <- append(b1s, B1)
          }
          if (any(which(avs == min(avs)) %in% c(3, 4))) {
            tc <- round(tc, 1)
            B1 <- max(substi$Radius[which(tc[, 2] == 0, arr.ind = T)])
            B1.loc <- substi$L[which.max(substi$Radius[which(tc[, 2] == 0, arr.ind = T)])]
            b1s.loc <- append(b1s.loc, B1.loc)
            b1s <- append(b1s, B1)
          }
        } else {
          if (any(which(avs == min(avs)) %in% c(1, 2))) {
            ind <- (which(abs(tc[, 1]) == min(avs), arr.ind = T))[1]
            against <- vector()
            against.loc <- vector()
            if (any(tc[ind, 1] < 0)) {
              new.ind <- which(tc[, 1] == min(tc[, 1]), arr.ind = T)
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 1], tc[new.ind[1], 1], tc[new.ind[1], 1] + 1))
              tc[,1] <- -tc[,1]
            } else {
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 1], tc[ind[1], 1] - 1, tc[ind[1], 1]))
            }
            for (i in 2:dim(tc)[1]) {
              if (tc[i, 3] == T) {
                against <- append(against, tc[i, 1] + substi$Radius[i])
                against.loc <- append(against.loc, substi$L[i])
              }
              B1 <- ifelse(length(against) > 0, max(against), abs(tc[ind[1], 1]) + substi$Radius[ind[1]])
              B1.loc <- ifelse(length(against) > 0, against.loc[which.max(against)], substi$L[ind[1]])
            }
            b1s <- append(b1s, B1)
            b1s.loc <- append(b1s.loc, B1.loc)
          }
          if (any(which(avs == min(avs)) %in% c(3, 4))) {
            ind <- (which(abs(tc[, 2]) == min(avs), arr.ind = T))[1]
            against <- vector()
            against.loc <- vector()
            if (any(tc[ind, 2] < 0)) {
              new.ind <- which(tc[, 2] == min(tc[, 2]), arr.ind = T)
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 2], tc[new.ind[1], 2], tc[new.ind[1], 2] + 1))
              tc[,2] <- -tc[,2]
            } else {
              tc <- plyr::mutate(as.data.frame(tc), "indi" = dplyr::between(tc[, 2], tc[ind[1], 2] - 1, tc[ind[1], 2]))
            }
            for (i in 2:dim(tc)[1]) {
              if (tc[i, 3] == T) {
                against <- append(against, tc[i, 2] + substi$Radius[i])
                against.loc <- append(against.loc, substi$L[i])
              }
              B1 <- ifelse(length(against) > 0, max(against), abs(tc[ind[1], 2]) + substi$Radius[ind[1]])
              B1.loc <- ifelse(length(against) > 0, against.loc[which.max(against)], substi$L[ind[1]])
            }
            b1s <- append(b1s, B1)
            b1s.loc <- append(b1s.loc, B1.loc)
          }
        }
      }
      B1 <- min(b1s[b1s >= 0])
      loc.B1 <- max(b1s.loc[which(b1s[b1s >= 0] == min(b1s[b1s >= 0], na.rm = TRUE))])
      B5 <- max(substi$Bs)
      L <- max(substi$L) + 0.4
      loc.B5 <- min(substi$V3[which(substi$Bs == B5)])
      result <- unique(round(data.frame(B1, B5, L, loc.B5,loc.B1), 2))
      unlink(list.files(pattern = ".xyz"))
      setwd("..")
      return(result)
    }, warning = function(w) {
      print(getwd())
      message("Something is wrong")
      unlink(list.files(pattern = ".xyz"))
      setwd("..")
    }, error = function(e) {
      print(basename(getwd()))
      message("Something is wrong, stopping steRimol")
      unlink(list.files(pattern = ".xyz"))
      setwd("..")
    }
  )
}

steRimol.df <- function(coor.atoms, radii = 'CPK', only.sub = T) {
  molecules <- list.dirs(full.names = F, recursive = F)
  steri.list <- list()
  for (molecule in molecules) {
    steri.list[[match(molecule, molecules)]] <- steRimol(molecule, coor.atoms, radii, only.sub)
  }
  steri.dafr <- data.frame(data.table::rbindlist(steri.list))
  row.names(steri.dafr) <- molecules
  return(steri.dafr)
}

diversitree::set.defaults(steRimol.df, getwd())

### 4 Dipole moment (dip.gaussian & dip.gaussian.df [substituted] - npa_dipole & npa_dipole.df )

Center.Of.Mass <- function(coor.atoms, mol.dir) {
  setwd('..')
  xyz_file_generator(list.files(pattern = "basic"))
  setwd(list.files(pattern = "basic"))
  coor.trans(list.files(pattern = ".xyz"),coor.atoms)
  xyz <- data.table::fread(list.files(pattern = "tc.xyz"), header = F, colClasses = c("character", "numeric", "numeric", "numeric"),
                           stringsAsFactors = F)
  suppressMessages(xyz$V1 <- plyr::mapvalues(xyz$V1,
                                             from = c(
                                               "H", "B", "C", "N", "O", "F", "Si", "P",
                                               "S", "Cl", "Br", "I"
                                             ),
                                             to = c(1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53)
  ))
  xyz <- data.frame(transform(xyz, V1 = as.numeric(V1)))
  for (i in 1:dim(xyz)[1]) ifelse(xyz[i,1] > 1, xyz[i,1] <- xyz[i,1] * 2, 1)
  M <- sum(xyz[,1])
  for (i in 1:dim(xyz)[1]) xyz[i,2:4] <- xyz[i,1] * xyz[i,2:4]
  com <- c(sum(xyz[,2]), sum(xyz[,3]), sum(xyz[,4]))
  com <- (1/M)*com
  unlink(list.files(pattern = ".xyz"))
  setwd('..')
  setwd(mol.dir)
  return(com)
}

Center.Of.Mass.Substructure <- function(coor.atoms, mol.dir, sub.atoms) {
  setwd('..')
  atoms <- strsplit(sub.atoms, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  xyz_file_generator(list.files(pattern = 'basic'))
  setwd(list.files(pattern = "basic"))
  coor.trans(list.files(pattern = ".xyz"),coor.atoms)
  xyz <- data.table::fread(list.files(pattern = "tc.xyz"), header = F, colClasses = c("character", "numeric", "numeric", "numeric"),
                           stringsAsFactors = F)
  suppressMessages(xyz$V1 <- plyr::mapvalues(xyz$V1,
                                             from = c(
                                               "H", "B", "C", "N", "O", "F", "Si", "P",
                                               "S", "Cl", "Br", "I"
                                             ),
                                             to = c(1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53)
  ))
  xyz <- data.frame(transform(xyz, V1 = as.numeric(V1)))
  xyz <- xyz[numeric.atoms,]
  for (i in 1:dim(xyz)[1]) ifelse(xyz[i,1] > 1, xyz[i,1] <- xyz[i,1] * 2, 1)
  M <- sum(xyz[,1])
  for (i in 1:dim(xyz)[1]) xyz[i,2:4] <- xyz[i,1] * xyz[i,2:4]
  com <- c(sum(xyz[,2]), sum(xyz[,3]), sum(xyz[,4]))
  com <- (1/M)*com
  unlink(list.files(pattern = '*.xyz'))
  setwd('..')
  setwd(mol.dir)
  return(com)
}

dip.xTB <- function(mol.dir, coor.atoms = '', center.of.mass = 'F', center.of.mass.substructure = 'F', sub.atoms = NULL) {
  setwd(mol.dir)
  dip <- list.files(pattern = "dipole")
  dipole <- data.frame(data.table::fread(dip,
                                         sep = " ", header = F, fill = T
  ))
  dipole <- dipole[, 2:5]
  names(dipole) <- c('dip_x', 'dip_y', 'dip_z', 'Total')
  if (!(coor.atoms == '')) {
    atoms <- strsplit(coor.atoms, " ")
    unlisted.atoms <- unlist(atoms)
    numeric.atoms <- as.numeric(unlisted.atoms)
    xyz <- data.frame(read.csv(list.files(pattern = "xyz_"), header = F, col.names = c('atom','x','y','z')))
    xyz <- xyz[-1, -1]
    xyz <- sapply(xyz, as.numeric)
    mag <- function(vector) {
      sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
    }
    
    if (length(numeric.atoms) == 4) {
      new_origin <- (xyz[numeric.atoms[[1]], ] + xyz[numeric.atoms[[2]], ])/2
      new_y <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                            mag(xyz[numeric.atoms[[3]], ] - new_origin))
      coplane <- as.numeric((xyz[numeric.atoms[[4]], ] - new_origin) /
                              mag(xyz[numeric.atoms[[4]], ] - new_origin))
    } else {
      if (center.of.mass == 'T') {
        com <- Center.Of.Mass(coor.atoms, mol.dir)
        new_origin <- com
        new_y <- as.numeric((xyz[numeric.atoms[[2]], ] - new_origin) /
                              mag(xyz[numeric.atoms[[2]], ] - new_origin))
        coplane <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                                mag(xyz[numeric.atoms[[3]], ] - new_origin))
      }
      if (center.of.mass.substructure == 'T') {
        com <- Center.Of.Mass.Substructure(coor.atoms, mol.dir, sub.atoms)
        new_origin <- com
        new_y <- as.numeric((xyz[numeric.atoms[[2]], ] - new_origin) /
                              mag(xyz[numeric.atoms[[2]], ] - new_origin))
        coplane <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                                mag(xyz[numeric.atoms[[3]], ] - new_origin))
      }
      if (center.of.mass.substructure == 'F' && center.of.mass == 'F') {
        new_origin <- xyz[numeric.atoms[[1]], ]
        new_y <- as.numeric((xyz[numeric.atoms[[2]], ] - new_origin) /
                              mag(xyz[numeric.atoms[[2]], ] - new_origin))
        coplane <- as.numeric((xyz[numeric.atoms[[3]], ] - new_origin) /
                                mag(xyz[numeric.atoms[[3]], ] - new_origin))
      }
    }
    
    cross_y_coplane <- pracma::cross(coplane, new_y)
    coef_mat <- aperm(array(c(
      new_y,
      coplane,
      cross_y_coplane), dim = c(3, 3)))
    
    angle_new.y_coplane <- angle(coplane, new_y)
    x_ang_new.y <- pi / 2
    cop_ang_x <- angle_new.y_coplane - x_ang_new.y
    result_vec <- c(0, cos(cop_ang_x), 0)
    new_x <- solve(coef_mat, result_vec)
    new_z <- pracma::cross(new_x, new_y)
    new_basis <- aperm(array(c(new_x, new_y, new_z), dim = c(3, 3)))
    dipole[1, 1:3] <- round(aperm(new_basis %*% as.numeric(dipole[1,1:3 ])), 4)
  }
  setwd('..')
  dipole <- lapply(dipole, as.numeric)
  return(dipole)
}

dip.xTB.df <- function(coor.atoms = '', center.of.mass = 'F', center.of.mass.substructure = 'F', sub.atoms = NULL) {
  molecules <- list.dirs(recursive = F, full.names = F)
  dip.list <- list()
  for (molecule in molecules) {
    dip.list[[match(molecule, molecules)]] <- dip.xTB(molecule, coor.atoms, center.of.mass, center.of.mass.substructure, sub.atoms)
  }
  dip.dafr <- data.frame(data.table::rbindlist(dip.list, fill = T))
  row.names(dip.dafr) <- molecules
  return(dip.dafr)
}

# diversitree::set.defaults(dip.gaussian.df, getwd())

# 5 Dihedrals, angles and bond lengths (bond.lengths)

angles <- function(atoms) {
  atoms.vec <- strsplit(atoms, " ")
  unlisted.atoms <- unlist(atoms.vec)
  numeric.atoms <- as.numeric(unlisted.atoms)
  ang.or.dih <- length(numeric.atoms)
  if (length(numeric.atoms) == 3) numeric.atoms <- c(numeric.atoms[1],
                                                     rep(numeric.atoms[2], 2),
                                                     numeric.atoms[3])
  molecules <- list.dirs(full.names = F,recursive = F)
  angle.df <- data.frame(matrix(ncol = 1, nrow = length(molecules)))
  colnames(angle.df)[1] <- ifelse(ang.or.dih == 3, paste('Angle(', stringr::str_replace_all(atoms, ' ', '_'),')', sep = ''),
                                  paste('Dihedral(', stringr::str_replace_all(atoms, ' ', '_'),')', sep = ''))
  for (molecule in molecules) {
    xyz_file_generator(molecule)
    setwd(molecule)
    xyz <- data.table::fread(list.files(pattern = '.xyz'))
    xyz <- xyz[, -1]
    first.bond <- xyz[numeric.atoms[1]] - xyz[numeric.atoms[2]]
    second.bond <- xyz[numeric.atoms[4]] - xyz[numeric.atoms[3]]
    angle.df[match(molecule, molecules), 1] <- angle(as.numeric(first.bond),
                                                     as.numeric(second.bond)) * (180/pi)
    unlink(list.files(pattern = '.xyz'))
    setwd('..')
  }
  return(angle.df)
}

df.angles <- function(atoms.vector) {
  ang.list <- list()
  for (i in 1:length(atoms.vector)) {
    ang.list[[match(i, 1:length(atoms.vector))]] <-  angles(atoms.vector[i])
  }
  ang.df <- do.call(cbind, ang.list)
  return(ang.df)
}

bond.lengths <- function(atom.pairs) {
  bonds.vec <- strsplit(atom.pairs, " ")
  unlisted.bvec <- unlist(bonds.vec)
  numeric.bvec <- as.numeric(unlisted.bvec)
  paired <- split(
    numeric.bvec,
    ceiling(seq_along(numeric.bvec) / 2)
  )
  molecules <- list.dirs(full.names = F,recursive = F)
  mag <- function(vector) {
    sqrt(vector[[1]]^2 + vector[[2]]^2 + vector[[3]]^2)
  }
  dist.df <- data.frame(matrix(ncol = length(paired), nrow = length(molecules)))
  names(dist.df) <- stringr::str_replace(as.character(paired),'c','Dist')
  for (molecule in molecules) {
    xyz_file_generator(molecule)
    setwd(molecule)
    xyz <- data.table::fread(list.files(pattern = '.xyz'))
    dist.list <- data.frame(matrix(ncol = length(paired), nrow = 1))
    for (i in 1:length(paired)) {
      dist.list[i] <- mag(xyz[paired[[i]][1], 2:4] - xyz[paired[[i]][2], 2:4])
    }
    names(dist.list) <- paired
    dist.df[which(molecules == stringr::str_remove(molecule,'/')), ] <- dist.list
    unlink(list.files(pattern = '.xyz'))
    setwd('..')
  }
  for (file in list.files(full.names = F, recursive = F)) {
    setwd(file)
    unlink(list.files(pattern = ".xyz"))
    setwd('..')
  }
  row.names(dist.df) <- molecules
  return(dist.df)
}

# 6 Vibrations (dot.prod.info, atomic.vec, gen.vib, df.gen.vib, ring.vibs and ring.vibs.df)

a.vec <- function(atom) {
  movec <- data.table::fread(list.files(pattern = 'movec'), fill = T)
  movec <- movec[-c(1:which(movec$V1 == 'reduced')), ]
  movec <- movec[-which(movec$V1 == 'reduced' |
                          movec$V1 == 'Frequencies' |
                          movec$V1 == 'Red.' |
                          movec$V1 == 'Frc' |
                          movec$V1 == 'IR' |
                          movec$V1 == 'Raman' |
                          movec$V1 == 'Depolar' |
                          movec$V1 == 'Atom'), ]
  movec <- movec[movec$V4 != '']
  movec <- movec[movec$V1 == atom]
  movec <- movec[, -c(1:2)]
  names(movec) <- c("x", "y", "z", "x", "y", "z", "x", "y", "z")
  iterations <- dim(movec)[1]
  datalist <- list()
  for (i in 1:iterations) {
    a <- movec[i, 1:3]
    b <- movec[i, 4:6]
    c <- movec[i, 7:9]
    d <- rbind(a, b, c)
    datalist[[i]] <- assign(
      paste("out", i, sep = "."),
      data.frame(d, row.names = c(1, 2, 3))
    )
  }
  output <- do.call(rbind, datalist)
  row.names(output) <- seq(1, iterations * 3, 1)
  output <- data.frame(apply(output, 2, as.numeric))
  return(output)
} 

prod.info <- function() {
  vibe.spec <- data.table::fread(list.files(pattern = 'vibspectrum'), fill = T)
  vibe.spec <- vibe.spec[-c(1:9, dim(vibe.spec)[1]), ]
  vibe.spec <- vibe.spec[,3:4]
  vibe.spec <- data.table::transpose(data.frame(apply(vibe.spec, 2, as.numeric)))
  colnames(vibe.spec) <- 1:ncol(vibe.spec)
  row.names(vibe.spec) <- c('Freq', 'Inten')
  return(vibe.spec)
}

gen.vib <- function(directory, atom.pairs, threshold = 1300) {
  setwd(directory)
  xyz <- data.frame(read.csv(list.files(pattern = "xyz_"), header = F, col.names = c('atom','x','y','z')))
  xyz <- xyz[-1, -1]
  xyz <- data.frame(sapply(xyz, as.numeric))
  colnames(xyz) <- c("x", "y", "z")
  bonds.vec <- strsplit(atom.pairs, " ")
  unlisted.bvec <- unlist(bonds.vec)
  numeric.bvec <- as.numeric(unlisted.bvec)
  paired <- split(
    numeric.bvec,
    ceiling(seq_along(numeric.bvec) / 2)
  )
  bond_info <- split(
    unlisted.bvec,
    ceiling(seq_along(unlisted.bvec) / 2)
  )
  bond.pairs <- list()
  for (i in 1:length(bond_info)) {
    bond.pairs[[i]] <- paste(bond_info[[i]][[1]],
                             bond_info[[i]][[2]],
                             sep = ","
    )
  }
  bonds <- data.table::fread(list.files(pattern = 'wbo'))
  bonds <- bonds[bonds$V3 > 0.8, ]
  bonds <- bonds[, 1:2]
  bonds <- data.frame(sapply(bonds, as.numeric))
  real <- bonds
  real.list <- list()
  for (i in 1:dim(real)[[1]]) {
    real.list[[i]] <- paste(as.character(real[i, 1]),
                            as.character(real[i, 2]),
                            sep = ","
    )
  }
  check <- bond.pairs %in% real.list
  if (all(check) == T) {
    vec.list <- list()
    for (atom in numeric.bvec) {
      vec.list[[match(atom, numeric.bvec)]] <- a.vec(atom)
      names(vec.list)[[match(atom, numeric.bvec)]] <- as.character(atom)
    }
    units <- list()
    for (i in 1:length(paired)) {
      units[[i]] <- cbind(
        vec.list[[as.character(paired[[i]][[1]])]],
        vec.list[[as.character(paired[[i]][[2]])]]
      )
      units[[i]][, 7:9] <- xyz[paired[[i]][[1]], ] - xyz[paired[[i]][[2]], ]
    }
    info <- prod.info()
    max.list <- list()
    for (i in 1:length(units)) {
      units[[i]][, 10] <- NA
      units[[i]][, 11] <- data.table::transpose(info[1, ])
      units[[i]] <- units[[i]][units[[i]]$V1 > threshold, ]
      rows <- dim(units[[i]])[1]
      for (j in 1:rows) {
        units[[i]][j, 10] <- abs(pracma::dot(as.matrix(units[[i]][j, 1:3]), as.matrix(units[[i]][j, 7:9]))) +
          abs(pracma::dot(as.matrix(units[[i]][j, 4:6]), as.matrix(units[[i]][j, 7:9])))
      }
    }
    for (i in 1:length(units)) {
      max.list[[i]] <- which.max(units[[i]][, 10])
    }
    max.list <- data.frame(max.list)
    uni <- data.frame(matrix(nrow = 1, ncol = dim(max.list)[2]))
    for (i in 1:dim(max.list)[2]) {
      uni[1, i] <- units[[i]][max.list[1, i], 11]
    }
    for (i in 1:dim(uni)[2]) {
      colnames(uni)[[i]] <- gsub("\\D+", "-", as.character(paired[i]))
    }
    row.names(uni) <- directory
    setwd('..')
    return(uni)
  }
  else {
    print(directory)
    print("The following bonds do not exist - check atom numbering")
    print(bond.pairs[match(F, check)])
    setwd("..")
  }
}

df.gen.vib <- function(atom.pairs, threshold = 1350) {
  molecules <- list.files()
  mol.list <- list()
  for (molecule in molecules) {
    mol.list[[match(molecule, molecules)]] <- gen.vib(molecule, atom.pairs, threshold)
  }
  vib.dafr <- do.call(rbind, mol.list)
  row.names(vib.dafr) <- molecules
  return(vib.dafr)
}

diversitree::set.defaults(df.gen.vib, getwd())

ring.vibs <- function(mol.dir, ordered.ring.atoms) {
  setwd(mol.dir)
  freq <- prod.info()
  atoms <- strsplit(ordered.ring.atoms, " ")
  unlisted.atoms <- unlist(atoms)
  numeric.atoms <- as.numeric(unlisted.atoms)
  paired.atoms <- split(
    numeric.atoms,
    ceiling(seq_along(numeric.atoms) / 2)
  )
  atoms_info <- split(
    unlisted.atoms,
    ceiling(seq_along(unlisted.atoms) / 2)
  )
  xyz <- data.frame(read.csv(list.files(pattern = "xyz_"), header = F, col.names = c('atom','x','y','z')))
  xyz <- xyz[-1, -1]
  xyz <- data.frame(sapply(xyz, as.numeric))
  xyz <- xyz[numeric.atoms, ]
  colnames(xyz) <- c("x", "y", "z")
  pa <- xyz[1, 1:3] - xyz[2, 1:3]
  atom.one <- a.vec(atoms_info[[1]][[1]])
  atom.two <- a.vec(atoms_info[[3]][[1]])
  atom.three <- a.vec(atoms_info[[2]][[1]])
  atom.four <- a.vec(atoms_info[[1]][[2]])
  atom.five <- a.vec(atoms_info[[2]][[2]])
  atom.six <- a.vec(atoms_info[[3]][[2]])
  sum.1.3.5 <- list()
  for (i in 1:dim(atom.one)[1]) {
    sum.1.3.5[[i]] <- atom.one[i, ] + atom.three[i, ] + atom.five[i, ]
  }
  vec.sum.1.3.5 <- data.frame(do.call(rbind, sum.1.3.5))
  sum.2.4.6 <- list()
  for (i in 1:dim(atom.two)[1]) {
    sum.2.4.6[[i]] <- atom.two[i, ] + atom.four[i, ] + atom.six[i, ]
  }
  vec.sum.2.4.6 <- data.frame(do.call(rbind, sum.2.4.6))
  prods <- list()
  for (i in 1:dim(vec.sum.1.3.5)[1]) {
    prods[[i]] <- pracma::dot(as.numeric(vec.sum.1.3.5[i, ]), as.numeric(vec.sum.2.4.6[i, ]))
  }
  prod.vec.sums <- data.frame(do.call(rbind, prods))
  prod.vec.sums[, 2] <- data.table::transpose(as.list(freq[1, ]))
  for (i in 1:dim(prod.vec.sums)[1]) {
    prod.vec.sums[i, 3] <- abs(sin(angle(as.numeric(vec.sum.2.4.6[i, ]), as.numeric(pa))))
  }
  pvec.filter <- dplyr::filter(prod.vec.sums, prod.vec.sums$do.call.rbind..prods. != 0)
  vec.prod.filtered <- dplyr::filter(pvec.filter, pvec.filter$V2 > 1500 &
                                       abs(pvec.filter$do.call.rbind..prods.) > 0.1)
  result <- data.frame(
    vec.prod.filtered[which.max(vec.prod.filtered[, 3]), 2],
    asin(max(vec.prod.filtered[, 3]))*(180/pi),
    vec.prod.filtered[which.min(vec.prod.filtered[, 3]), 2],
    asin(min(vec.prod.filtered[, 3]))*(180/pi)
  )
  colnames(result) <- c("cross",'cross.angle', "para", 'para.angle')
  if (dim(result)[1] == 0) {
    vec.prod.filtered <- dplyr::filter(pvec.filter, pvec.filter$V2 > 1550 &
                                         pvec.filter$V2 < 1750)
    result <- data.frame(
      vec.prod.filtered[which.max(vec.prod.filtered[, 3]), 2],
      acos(max(vec.prod.filtered[, 3]))*(180/pi),
      vec.prod.filtered[which.min(vec.prod.filtered[, 3]), 2],
      acos(min(vec.prod.filtered[, 3]))*(180/pi)
    )
    colnames(result) <- c("cross",'cross.angle', "para", 'para.angle')
    print(basename(getwd()))
    print('Dot products are lower than 0.1 - returning the default 1750 - 1550 1/cm')
  }
  setwd("..")
  return(result)
}

ring.vib.df <- function(input) {
  molecules <- list.dirs(recursive = F, full.names = F)
  ring.vib.list <- list()
  for (molecule in molecules) {
    ring.vib.list[[match(molecule, molecules)]] <- ring.vibs(molecule, input)
  }
  ring.vib.dafr <- data.frame(data.table::rbindlist(ring.vib.list))
  row.names(ring.vib.dafr) <- molecules
  return(ring.vib.dafr)
}

diversitree::set.defaults(ring.vib.df, getwd())

# 8 check for imaginary frequencies and report

ground.state <- function(path) {
  molecules <- list.dirs(path, recursive = F, full.names = F)
  Ground_State <- vector()
  for (molecule in molecules) {
    setwd(molecule)
    imaginary <- prod.info()
    Ground_State[[match(molecule, molecules)]] <- all(!(imaginary[1, ] < 0))
    setwd('..')
  }
  imag.report.dafr <- data.frame(Ground_State)
  row.names(imag.report.dafr) <- molecules
  return(imag.report.dafr)
}

# 7 Produce data set

comp.set <- function(input = NULL, atom.pairs, coor.atoms.dip, center.of.mass = 'F', center.of.mass.substructure = 'F', sub.atoms = NULL, 
                     coor.atoms.ste, atom.dist, ang.atoms) {
  molecules <- list.dirs(recursive = F, full.names = F)
  if (!is.null(input)) ring <- ring.vib.df(input)
  vib <- df.gen.vib(atom.pairs)
  dip <- dip.xTB.df(coor.atoms.dip, center.of.mass, center.of.mass.substructure, sub.atoms)
  ste <- steRimol.df(radii = 'CPK', coor.atoms.ste)
  B.L <- bond.lengths(atom.dist)
  angs <- df.angles(unlist(strsplit(ang.atoms, '  ')))
  if (!is.null(input)) {
    comp <- dplyr::bind_cols(ring, vib, B.L, angs, dip, ste)
  } else {
    comp <- dplyr::bind_cols(vib, B.L, angs, dip, ste)
  }
  row.names(comp) <- molecules
  for (i in 1:length(comp)) {
    if (class(comp[[i]]) == "character") {
      comp[[i]] <- as.numeric(comp[[i]])
    }
  }
  return(comp)
}


comp.set.hetero <- function(path, dipole.mode = 'xtb', radii = 'CPK') {
  molecules <- list.dirs(recursive = F, full.names = F)
  vib <- df.gen.vib()
  if (dipole.mode == 'compute') {
    dip <- npa_dipole.df()
  } 
  if (dipole.mode == 'xtb') {
    dip <- dip.xTB.df()
  }
  if (radii == 'CPK') {
    ste <- steRimol.df(radii = 'CPK')
  } 
  if (radii == 'bondi') {
    ste <- steRimol.df(radii = 'bondi')
  }
  atom.dist <- readline('Distances - Atom pairs: ')
  B.L <- bond.lengths(atom.dist)
  angs <- df.angles(ang.atoms)
  comp <- dplyr::bind_cols(vib, dip, B.L, angs, dip, ste)
  row.names(comp) <- molecules
  for (i in 1:length(comp)) {
    if (class(comp[[i]]) == "character") {
      comp[[i]] <- as.numeric(comp[[i]])
    }
  }
  return(comp)
}
diversitree::set.defaults(comp.set.hetero, getwd())

# Produce multiple sets including experimental observation and models, coefficients and plots

batch.comp.set <- function(input_file = NULL) {
  if (is.null(input_file)) {
    input <- readline("Ring atoms - by order -> primary axis (para first), ortho atoms and meta atoms: ")
    atom.pairs <- readline("Your atom pairs: ")
    coor.atoms.dip <- readline('Dipole coordinates atoms: ')
    center.of.mass <- readline('Use center of mass as origin (T/F)? ')
    center.of.mass.substructure <- readline('Use center of mass of a subset of atoms as origin (T/F)? ')
    if (center.of.mass.substructure == 'T') {
      sub.atoms <- readline("Indicate the subset's numbering (in the basic structure: ")
    }
    coor.atoms.ste <- readline("Primary axis along: ")
    atom.dist <- readline('Distances - Atom pairs: ')
    ang.atoms <- readline('Angles and Dihedrals - Atom Triad and Quartets: ')
  } else {
    input_file <- read.csv(input_file, header = F, fill = T, sep = ' ')
    input_file <- dplyr::mutate(input_file, last = rep(NA, dim(input_file)[1]))
    input <- do.call(paste, input_file[1, 1:(which(is.na(input_file[1, ]))[1] - 1)])
    atom.pairs <- do.call(paste, input_file[2, 1:(which(is.na(input_file[2, ]))[1] - 1)])
    coor.atoms.dip <- do.call(paste, input_file[3, 1:(which(is.na(input_file[3, ]))[1] - 1)])
    center.of.mass <- input_file[4, 1:(which(is.na(input_file[4, ]))[1] - 1)]
    center.of.mass.substructure <- input_file[5, 1]
    if (center.of.mass.substructure == 'T') {
      sub.atoms <- do.call(paste, input_file[5, 2:(which(is.na(input_file[1, ]))[1] - 1)])
    }
    coor.atoms.ste <- do.call(paste, input_file[6, 1:(which(is.na(input_file[6, ]))[1] - 1)])
    atom.dist <- do.call(paste, input_file[7, 1:(which(is.na(input_file[7, ]))[1] - 1)])
    angs <- paste(do.call(paste, input_file[8, 1:(which(is.na(input_file[8, ]))[1] - 1)]),
                          do.call(paste, input_file[8, (which(is.na(input_file[8, ]))[1] + 1):which(is.na(input_file[8, ]))[2] - 1])),
                  sep = '  ')
  }
  conformers <- list.files()
  sets <- list()
  for (i in conformers) {
    setwd(i)
    setwd('extracted_info')
    name_changer('_close','')
    name_changer('_far','')
    name_changer('_num','')
    name_changer('_tc','')
    name_changer('_fixed','')
    name_changer('_0','')
    sets[[match(i, conformers)]] <- comp.set(input, atom.pairs, coor.atoms.dip,
                                             center.of.mass, center.of.mass.substructure,
                                             sub.atoms, coor.atoms.ste, atom.dist, dih, ang.atoms)
    setwd('../..')
  }
  names(sets) <- conformers
  return(sets)
}

modelling.data <- function(datasets, output_file) {
  print('make sure you are in a neutral directory - this creates a lot of csv files!')
  output = data.table::fread(output_file)[, 2]
  dir.name <- readline('Name the datasets directory: ')
  dir.create(dir.name)
  setwd(dir.name)
  names(output) <- 'output'
  ready.for.modelling <- lapply(datasets, function(x) dplyr::mutate(x, output))
  for (i in 1:length(ready.for.modelling)) {
    write.csv(ready.for.modelling[[i]], paste(names(ready.for.modelling)[[i]], '.csv', sep = ''))
  }
  setwd('..')
}


model.plot <- function(model = models[1,1], data) {
  best.mod <- lm(model, data = data)
  pred_interval <- predict(best.mod, newdata = data, interval = 'pre', level = 0.9)
  plot.dat <- data.frame(cbind(data$output, pred_interval))
  colnames(plot.dat) <- c('Measured', 'Predicted', 'lwr', 'upr')
  rownames(plot.dat) <- row.names(data)
  
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o_",'2-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m_",'3-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"p_",'4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"o4-",'2,4-')
  row.names(plot.dat) <- stringr::str_replace(row.names(plot.dat),"m3-",'3,3-')
  
  for (i in 1:nrow(plot.dat)) {
    if (grepl('3-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'meta'
    }
    if (grepl('2-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'ortho'
    }
    if (grepl('basic',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'Ph'
    }
    if (grepl('4-',row.names(plot.dat)[i])) {
      plot.dat[i,5] <- 'para'
    }
  }
  
  plot.dat <- dplyr::mutate(plot.dat, label = row.names(plot.dat))
  if (ncol(plot.dat) > 5) {
    colnames(plot.dat)[5] <- 'Position'
  } else {
    plot.dat <- dplyr::mutate(plot.dat, Position = rep('black', nrow(plot.dat)))
  }
  
  plot <- ggplot(plot.dat, aes(x = Measured, y = Predicted)) +
    geom_point(size = 3, shape = 15,aes(color = Position)) +
    stat_smooth(aes(y = lwr), color = "cadetblue", linetype = "dashed",
                se = F, method = 'lm', fullrange = T, size = 0.8) +
    stat_smooth(aes(y = upr), color = "cadetblue", linetype = "dashed", 
                se = F, method = 'lm', fullrange = T, size = 0.8) +
    labs(x = 'Measured %',y = 'Predicted %') +
    stat_smooth(method = 'lm',se = F, formula = y~x,
                color = 'black',fullrange = T, linetype = 'dashed') +
    theme(axis.line.x = element_line(size = 1, colour = "black"),
          axis.line.y = element_line(size = 1, colour = "black"),
          axis.text.x = element_text(colour = "black", size = 12,face = 'bold'),
          axis.text.y = element_text(colour = "black", size = 12,face = 'bold'),
          axis.title.x = element_text(colour = "black", size = 12,face = 'bold'),
          axis.title.y = element_text(colour = "black", size = 12,face = 'bold'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(),
          legend.position = c(2,2)) +
    scale_color_manual(values = c(Ph = "black", meta = 'tan1',
                                  para = '#66a182',ortho = '#d1495b')) +
    xlim(min(plot.dat[,3]), max(plot.dat[,4])) +
    ylim(min(plot.dat[,3]), max(plot.dat[,4])) +
    coord_fixed(ratio = 1) +
    geom_text_repel(size = 5.2,
                    aes(label = label),
                    min.segment.length = Inf,
                    seed = 42,
                    point.padding = 0.4,
                    segment.color = 'white',
                    force_pull = 0.02,
                    nudge_x = 0.022,
                    direction = 'y') +
    theme(text = element_text(family = 'Helvetica'))
  plot
}

library(ggplot2)
library(ggrepel)
library(extrafont)

model.info <- function(dataset, min = 3, max = floor(dim(mod_data)[1] / 5), leave.out = '', predict = F) {
  print(dataset)
  mod_data <- data.frame(data.table::fread(dataset, header = T, check.names = T))
  RN <- mod_data[,1]
  mod_data <- mod_data[,-1]
  CN <- names(mod_data)
  mod_data <- data.frame(cbind(scale(mod_data[,1:dim(mod_data)[2] - 1], T, T), mod_data$output))
  names(mod_data) <- CN
  row.names(mod_data) <- RN
  pred.data <- mod_data[row.names(mod_data) == leave.out, ]
  mod_data <- mod_data[row.names(mod_data) != leave.out, ]
  models <- sub_model(mod_data, min = min, max = max)
  tab <- knitr::kable(models)
  print(tab)
  mod.sum <- summary(lm(models[1, 1], mod_data))
  k.mod <- knitr::kable(mod.sum$coefficients)
  print(k.mod)
  cv_3fold <- measure_cv(models[1,1], mod_data, dim(mod_data)[2], 3, 500)
  dt3 <- data.frame(cv_3fold[[2]], cv_3fold[[1]])
  names(dt3) <- c('Q2', 'MAE')
  tab_dt3 <- knitr::kable(dt3)
  print(tab_dt3)
  cv_5fold <- measure_cv(models[1,1], mod_data, dim(mod_data)[2], 5, 1000)
  dt5 <- data.frame(cv_5fold[[2]], cv_5fold[[1]])
  names(dt5) <- c('Q2', 'MAE')
  tab_dt5 <- knitr::kable(dt5)
  print(tab_dt5)
  if (predict == T) {
    prediction <- predict(lm(models[1, 1], mod_data), pred.data)
    real <- pred.data$output
    prd.tab <- data.frame(prediction, real)
    names(prd.tab) <- c('OOS Pred', 'OOS Measured')
    k.prd.tab <- knitr::kable(prd.tab)
    print(k.prd.tab)
  }
  model.plot(model = models[1, 1], data = mod_data)
}



