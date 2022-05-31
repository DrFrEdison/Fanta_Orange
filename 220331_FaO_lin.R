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
for(i in 1:length(fao$pls)) fao$pred.lin[[i]] <- produktion_prediction(csv_transfered = fao$lin, pls_function_obj = fao$pls[[i]], ncomp = fao$plspara$ncomp)

for(i in 1:length(fao$pls)) fao$pls_merge[[i]] <- lapply(fao$pred[[i]], function(x) .merge_pls(pls_pred = x, fao$pls_lm[[i]] ,R2=.5))
for(i in 1:length(fao$pls)) fao$pls_merge_site[[i]] <- .merge_pls_site(merge_pls_lm_predict_ls = fao$pls_merge[[i]], number = 5000, ncomp = fao$plspara$ncomp)

fao$pls_merge_site <- lapply(fao$pls_merge_site, function(x) x[x$spc == "spc" & x$ncomp < 5, ])
  
lin <- list()
lin$limit1 <- 2
lin$limit2 <- 8
lin$R2 <- .4
i = 1
for(i in 1:length(fao$pred.lin)){
  if(nrow(fao$pls_merge_site[[i]]) == 0) next
  lin$filter[[i]] <- linearitaet_filter(linearitaet_prediction = fao$pred.lin[[i]]$prediction
                                        , ncomp = fao$plspara$ncomp
                                        , linearitaet_limit_1 = lin$limit1
                                        , linearitaet_limit_2 = lin$limit2
                                        , R_2 = lin$R2
                                        , SOLL = seq(95,100,1)
                                        , pls_merge = fao$pls_merge_site[[i]])}
lin$filter

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

mod_c$pred$lin <- pred_of_new_model(modell_csv_transfered = fao$seqp$model[[ mod_c$matrix ]]
                                    , substance = fao$plspara$substance
                                    , wl1 = mod_c$wl1
                                    , wl2 = mod_c$wl2
                                    , wl3 = mod_c$wl3
                                    , wl4 = mod_c$wl4
                                    , ncomp = mod_c$ncomp
                                    , derivative = mod_c$spc
                                    , csv_transfered = fao$lin)

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

png(paste0(.date(),"_FaO_Prediction_lin_2mm.png"),xxx<-4800,xxx/16*9,"px",12,"white",res=500,"sans",T,"cairo")
par(mfrow = c(1,1), mar = c(5,4,5,1))
x.bias <- .bias(rev(fao$lin$data$Acid)[1],0,mod_c$pred$lin[1]) - 1

plot(rev(mod_c$pred$lin)  + x.bias 
     , axes = T, xlab = "", ylab = "Acid in %"
     , main = "Linearitaet in Dorsten DS am 20.04.2020,\nWerte umgedreht!"
     , pch = 20, cex = 1.5, col = "blue"
     , ylim = c(90, 100)
     , sub = paste("Bias =", round(x.bias,2)))

points(fao$lin$data$Acid, col = "red", pch = 20, cex = 1.5)

legend("topleft", c("Model", "Labor"), pch = 20, col = c("blue", "red"))
dev.off()

