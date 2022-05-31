dt <- list(); dt$R <- paste0(Sys.getenv("OneDriveCommercial"), "/FE_Methoden/Allgemein/R_dt_project/")
beverage <- "Fanta_Orange"
version <- "01"
source(paste0(dt$R,"R/wd.R"))
source(paste0(grep(beverage, wd$bev.ccep, value = T)[1], "/", paste0("Rsource_", beverage, "_V", version, ".R")))


# Parameter ####
dt$para$para = "Acid"
dt$para$unit <- bquote("%")
dt$para$ylab <- bquote("Acid in %")
dt$para$mop.date <- "220331"

# read validated production data ####
setwd(dt$wd)
setwd("./Modelloptimierung")
setwd(paste0("./", dt$para$mop.date, "_", dt$para$model.pl[1], "_", dt$para$para))
setwd("./Produktionsdaten")

fao <- list()

fao$para$main = "Fanta"
fao$para$bev = "FaO"

fao$files$ref <- grep("ref", dir(pattern = "ref.csv$")[grep(fao$para$main, dir(pattern = "ref.csv$"))], value = T) # Background spc 
fao$files$drk <- grep("rk", dir(pattern = "rk.csv$")[grep(fao$para$main, dir(pattern = "rk.csv$"))], value = T) # Dark spc
fao$files$spc <- grep("spc", dir(pattern = "spc.csv$")[grep(fao$para$main, dir(pattern = "spc.csv$"))], value = T) # Production spc

fao$para$beverage <- substr(fao$files$spc
                            , lapply(gregexpr("_", fao$files$spc), function(x) x[[2]] + 1)
                            , lapply(gregexpr("_", fao$files$spc), function(x) x[[7]] - 1))

fao$para$beverage <- substr(fao$para$beverage, 1, nchar(fao$para$beverage) - c(2, 3, 0, 2))

for(i in 1:length(fao$para$beverage))
  if(substr(fao$para$beverage, nchar(fao$para$beverage), nchar(fao$para$beverage))[i] == "_" ) fao$para$beverage[i] <- substr(fao$para$beverage[i], 1 , nchar(fao$para$beverage[i]) - 1)

# get file info
fao$txt <- lapply(fao$files, .txt.file)

# read files
fao$raw$ref <- lapply(fao$files$ref, function(x) fread(x, dec = ",", sep = ";")) # Background spc 
fao$raw$drk <- lapply(fao$files$drk, function(x) fread(x, dec = ",", sep = ";")) # Dark spc
fao$raw$spc <- lapply(fao$files$spc, function(x) fread(x, dec = ",", sep = ";")) # Production spc

# set names
names(fao$raw$ref) <- fao$txt$ref$loc.line
names(fao$raw$drk) <- fao$txt$drk$loc.line
names(fao$raw$spc) <- fao$txt$spc$loc.line

# read wavelength columns
fao$ppp <- lapply(fao$raw, function(x) lapply(x, .transfer_csv.num.col))

fao$trs$spc <- lapply(fao$raw$spc, .transfer_csv)

# read lin
setwd(dt$wd)
setwd("./Modelloptimierung")
setwd(paste0("./", dt$para$mop.date, "_", dt$para$model.pl[1], "_", dt$para$para))
setwd("./Linearitaet")

fao$lin <- fread("2020-04-22_CCEP_Dorsten_DS_Fanta_spc.csv", sep = ";", dec = ",")
fao$lin$Acid <- c(97.76,96.78,95.8,94.9,93.82,92.8)
fao$lin$Acidp <- 100:95
fao$lin <- fao$lin[order(fao$lin$Acidp) , ]
fao$lin[,.transfer_csv.num.col(fao$lin)$numcol] <- fao$lin[, lapply(.SD, function(x) x * 2 / .31), .SDcols = .transfer_csv.num.col(fao$lin)$numcol] 
fao$lin <- .transfer_csv(fao$lin)

# Model Matrix ####
setwd(dt$wd)
setwd("./Modelloptimierung")
setwd(paste0("./", dt$para$mop.date, "_", dt$para$model.pl[1], "_", dt$para$para))
setwd("./Modellmatrix")

fao$model <- fread("FantaO_2mm_kalset.txt", dec = ",")
names(fao$model)[ grep("Zitronens", names(fao$model)) ] <- "Acid"
fao$model <- fao$model[which(!is.na(fao$model$Acid)) , ]
fao$model <- .transfer_csv(csv.file = fao$model, data.table.set = T)

# PLS para####
fao$plspara$wlr <- .wlr_function(300:500, 300:500, 5)
nrow(fao$plspara$wlr)
fao$plspara$wlm <- .wlr_function_multi(300:500, 300:500, 10)
nrow(fao$plspara$wlm)
fao$plspara$wl <- rbind.fill(fao$plspara$wlm, fao$plspara$wlr)
nrow(fao$plspara$wl)

fao$plspara$ncomp <- 4
fao$plspara$substance <- "Acid"
fao$para$unit <- bquote("%")
fao$para$ylab <- bquote("Acid in %")

# model matrix variance ####
fao$seqp <- list()

# check for parameters to vary
apply(fao$model$data[,2:5], 2, unique)

# create variance
fao$seqp$Probe <- lapply(levels(factor(fao$model$dat$Ausmsichung)), function(x) grep(x, fao$model$data$Ausmsichung, invert = T))
fao$seqp$Probe[[4]] <- 1:nrow(fao$model$data)

fao$seqp$Probe[[1]] <- fao$seqp$Probe[[4]]
fao$seqp$Probe[[2]] <- NULL
fao$seqp$Probe[[2]] <- NULL
fao$seqp$Probe[[2]] <- NULL

names(fao$seqp$Probe) <- c("original")

# create model with list
fao$seqp$model <- lapply(fao$seqp$Probe, function(x) list(data = fao$model$data[x,]
                                                          , wl = fao$model$wl
                                                          , spc = fao$model$spc[x,]
                                                          , spc1st = fao$model$spc1st[x,]
                                                          , spc2nd = fao$model$spc2nd[x,]))

