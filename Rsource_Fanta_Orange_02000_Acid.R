dt <- list(); dt$R <- paste0(Sys.getenv("OneDriveCommercial"), "/FE_Methoden/", "Allgemein/R_dt_project/")
source(paste0(dt$R,"R/source_spc_files.R"))
source(paste0(dt$R,"R/source_pls.R"))
source(paste0(dt$R,"R/source_read.R"))

# general parameters ####
dt$para$customer = "CCEP"
dt$para$beverage = "Fanta_Orange"

setwd(paste0(dt$wd <- paste0(wd$fe$CCEP$Mastermodelle, dt$para$beverage)))
setwd( print( this.path::this.dir() ) )
setwd("..")
dt$wd.git <- print( getwd() )

dt$para$location = c("Dorsten", "Moenchengladbach", "Mannheim")
dt$para$line = c("DS", "G9", "MY")
dt$para$main = paste0(dt$para$beverage, " in ", dt$para$location, ", line ", dt$para$line)
dt$para$model.date <- c("220331")
dt$para$model.pl <- c("02000")
dt$para$wl1 <- c(190)
dt$para$wl2 <- c(598)
dt$para$wl[[1]] <- seq(dt$para$wl1, dt$para$wl2, 1)

dt$para$substance <- c("Acid")
dt$para$unit <- c( bquote("%") )
dt$para$ylab <- c( bquote("Acid in %" ) )
dt$para$mop.date <- "220331"
dt$para$SOLL <- c(100)

dt$para$mva.date <- "220429"
# rename R files (run only once)
# dt$para$Rfiles <- list.files(getwd(), pattern = ".R$", recursive = T)
# file.rename(dt$para$Rfiles, gsub("beverage", dt$para$beverage, dt$para$Rfiles))

# files ####
setwd(dt$wd)
setwd("./Modellvalidierung")
setwd("./Produktionsdaten")

fao <- list()

fao$files$ref <- grep("ref", dir(pattern = "ref.csv$")[grep(dt$para$beverage, dir(pattern = "ref.csv$"))], value = T) # Background spc 
fao$files$drk <- grep("rk", dir(pattern = "rk.csv$")[grep(dt$para$beverage, dir(pattern = "drk.csv$"))], value = T) # Dark spc
fao$files$spc <- grep("spc", dir(pattern = "spc.csv$")[grep(dt$para$beverage, dir(pattern = "spc.csv$"))], value = T) # Production spc

fao$txt <- lapply(fao$files, .txt.file)

# read ####
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

# read lin ####
setwd(dt$wd)
setwd("./Modellvalidierung")
setwd(paste0("./", dt$para$mva.date, "_", dt$para$model.pl[1], "_", dt$para$substance))
setwd("./Linearitaet")

fao$lin <- fread("2020-04-22_CCEP_Dorsten_DS_Fanta_spc.csv", sep = ";", dec = ",")
fao$lin$Acid <- c(97.76,96.78,95.8,94.9,93.82,92.8)
fao$lin$Acidp <- 100:95
# fao$lin <- fao$lin[order(fao$lin$Acidp) , ]
fao$lin[,.transfer_csv.num.col(fao$lin)$numcol] <- fao$lin[, lapply(.SD, function(x) x * 2 / .31), .SDcols = .transfer_csv.num.col(fao$lin)$numcol] 
fao$lin <- .transfer_csv(fao$lin)

# Model Matrix ####
setwd(dt$wd)
setwd("./Modelloptimierung")
setwd(paste0("./", dt$para$mop.date, "_", dt$para$model.pl[1], "_", dt$para$substance))
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

# model matrix variance ####
fao$seqp <- list()

# check for parameters to vary
apply(fao$model$data[,2:5], 2, unique)

# create variance
fao$seqp$Probe <- lapply(levels(factor(fao$model$dat$Ausmsichung)), function(x) grep(x, fao$model$data$Ausmsichung, invert = T))
fao$seqp$Probe[[4]] <- 1:nrow(fao$model$data)

names(fao$seqp$Probe) <- c(as.character(levels(factor(fao$model$dat$Ausmsichung))), "original")

# create model with list
fao$seqp$model <- lapply(fao$seqp$Probe, function(x) list(data = fao$model$data[x,]
                                                          , wl = fao$model$wl
                                                          , spc = fao$model$spc[x,]
                                                          , spc1st = fao$model$spc1st[x,]
                                                          , spc2nd = fao$model$spc2nd[x,]))