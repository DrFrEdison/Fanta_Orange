dt <- list(); dt$R <- paste0(Sys.getenv("OneDriveCommercial"), "/FE_Methoden/Allgemein/R_dt_project/")
beverage <- "Fanta_Orange"
version <- "01"
source(paste0(dt$R,"R/wd.R"))
source(paste0(grep(beverage, wd$bev.ccep, value = T)[1], "/", paste0("Rsource_", beverage, "_V", version, ".R")))


# Parameter ####
dt$para$para = "Coffein"
dt$para$unit <- bquote("%")
dt$para$ylab <- bquote("Coffein in %")
dt$para$mop.date <- "220413"

# Fanta Orange 2mm Model Modelloptimierung
# 22-03-31

# set wd and source packages and functions
dt <- list(); dt$R <- paste0(Sys.getenv("OneDriveCommercial"), "/FE_Methoden/", "Allgemein/R_dt_project/")
source(paste0(dt$R,"R/source_read.R"))
source(paste0(dt$R,"R/source_pls.R"))
dt$wd <- paste0(wd$fe$CCEP$Mastermodelle,"Fanta_Orange//Modelloptimierung/220331_02000_Acid/"); setwd(dt$wd)

dir()
source("220405_FaO_source.R")

# PLS and LM ####
fao$pls <- lapply(fao$seqp$model, function(x) pls_function(csv_transfered = x
                                                , substance = fao$plspara$substance
                                                , wlr = fao$plspara$wl 
                                                , ncomp = fao$plspara$ncomp))
gc()
for(i in 1:length(fao$pls)){
  for(j in 1:length(fao$seqp$model)){
    fao$pls_lm[[i]] <- pls_lm_function(fao$pls[[i]]
                                       , csv_transfered = fao$seqp$model[[j]]
                                       , substance = fao$plspara$substance
                                       , wlr = fao$plspara$wl 
                                       , ncomp = fao$plspara$ncomp)
  }
}

# Prediction ####
for(i in 1:length(fao$pls)) fao$pred[[i]] <- lapply(fao$trs$spc,function(x) produktion_prediction(csv_transfered = x, pls_function_obj = fao$pls[[i]], ncomp = fao$plspara$ncomp))

# Find the best model ####
for(i in 1:length(fao$pls)) fao$pls_merge[[i]] <- lapply(fao$pred[[i]], function(x) .merge_pls(pls_pred = x, fao$pls_lm[[i]] ,R2=.5))
for(i in 1:length(fao$pls)) fao$pls_merge_site[[i]] <- .merge_pls_site(merge_pls_lm_predict_ls = fao$pls_merge[[i]], number = 5000, ncomp = fao$plspara$ncomp)

fao$pls_merge_site <- lapply(fao$pls_merge_site, function(x) x[x$spc == "spc" & x$ncomp < 5, ])

lapply(fao$pls_merge_site, head)
  
matplot(fao$trs$spc$Dorsten_DS$wl
        , t(fao$trs$spc$Dorsten_DS$spc2nd)
        , type = "l", lty = 1)

# calculate best model ####
mod_c <- list()
mod_c$ncomp <- 2
mod_c$wl1 <- 410
mod_c$wl2 <- 460
mod_c$wl3 <- NA
mod_c$wl4 <- NA
mod_c$spc <- "spc"
mod_c$matrix <- 4
fao$seqp$listnames[mod_c$matrix]
mod_c$pngname <- paste(mod_c$spc, mod_c$wl1, mod_c$wl2, mod_c$wl3, mod_c$wl4, "PC", mod_c$ncomp, "matrix", mod_c$matrix, sep = "_")
mod_c$modelname <- paste(mod_c$wl1, mod_c$wl2, mod_c$wl3, mod_c$wl4, mod_c$spc, mod_c$ncomp, sep = "_")

mod_c$model <- pls_function(csv_transfered = fao$seqp$model[[ mod_c$matrix ]]
                            , substance = fao$plspara$substance
                            , wlr = data.frame(mod_c$wl1, mod_c$wl2, mod_c$wl3, mod_c$wl4)
                            , ncomp = mod_c$ncomp
                            , spc = mod_c$spc
                            , validation = "none")

mod_c$pred$spc <- list()
for(i in 1:length(fao$trs$spc))
mod_c$pred$spc[[i]] <- pred_of_new_model(modell_csv_transfered = fao$seqp$model[[ mod_c$matrix ]]
                                       , substance = fao$plspara$substance
                                       , wl1 = mod_c$wl1
                                       , wl2 = mod_c$wl2
                                       , wl3 = mod_c$wl3
                                       , wl4 = mod_c$wl4
                                       , ncomp = mod_c$ncomp
                                       , derivative = mod_c$spc
                                       , csv_transfered = fao$trs$spc[[i]])

for(i in 1:length(mod_c$pred$spc)) mod_c$pred$spc[[i]] <- as.numeric(ma(mod_c$pred$spc[[i]], 5))

setwd(dt$wd)
setwd("..")
setwd("./Analyse")

png(paste0(.date(),"_FaO_Prediction_2mm.png"),xxx<-4800,xxx/16*9,"px",12,"white",res=500,"sans",T,"cairo")
par(mfrow = c(2,2), mar = c(5,4,5,1))
for(i in 1:length(fao$trs$spc)){
  x.bias <- .bias(100,0,median(mod_c$pred$spc[[i]], na.rm = T))
  plot(mod_c$pred$spc[[i]] + x.bias
       , axes = F, xlab = "", ylab = "Acid in %"
       , main = paste(fao$para$beverage[i])
       , pch = 20, cex = .5
       , ylim = c(90, 110)
       , sub = paste("Bias =", round(x.bias,2)))
.xaxisdate(fao$trs$spc[[i]]$data$datetime)
}
dev.off()

setwd("..")
fao$lin$raw <- read.csv2("180525_Fanta_Orange_Linearitaet.txt")
fao$lin$ppp <- .transfer_csv.num.col(fao$lin$raw)
fao$lin$raw[ , fao$lin$ppp$numcol] <- fao$lin$raw[ , fao$lin$ppp$numcol] * 2/.2
fao$lin$trs <- .transfer_csv(fao$lin$raw)
matplot(fao$model$wl
        , t(fao$model$spc)
        , lty = 1, type = "l", col = "red")
matplot(fao$lin$ppp$wl
        , t(fao$lin$raw[ , fao$lin$ppp$numcol])
        , lty = 1, type = "l", col = "blue", add = T)
matplot(fao$trs$spc$Mannheim_MY$wl
        , t(fao$trs$spc$Mannheim_MY$spc)[ , seq(1, nrow(fao$trs$spc$Mannheim_MY$spc), 100)]
        , lty = 1, type = "l", col = "darkgreen", add = T)

mod_c$pred$lin <- pred_of_new_model(modell_csv_transfered = fao$seqp$model[[ mod_c$matrix ]]
                                         , substance = fao$plspara$substance
                                         , wl1 = mod_c$wl1
                                         , wl2 = mod_c$wl2
                                         , wl3 = mod_c$wl3
                                         , wl4 = mod_c$wl4
                                         , ncomp = mod_c$ncomp
                                         , derivative = mod_c$spc
                                         , csv_transfered = fao$lin$trs )


# final model matrix ####
setwd(dt$wd)
setwd("./Modellmatrix")

fao$seqp$export_matrix <- data.frame(cbind(fao$seqp$model[[ mod_c$matrix ]]$data, fao$seqp$model[[ mod_c$matrix ]]$spc))

write.csv2(fao$seqp$export_matrix
           , paste0(.date(), "_", fao$para$bev,"_",mod_c$pngname,"_mop_modelmatrix_final.csv")
           , row.names = F)


