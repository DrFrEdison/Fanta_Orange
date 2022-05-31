# beverage parameter ####
setwd(this.path::this.dir())
dir( pattern = "Rsource" )
source.file <- print( "Rsource_Fanta_Orange_02000_Acid.R" )
source( paste0(getwd(), "/", source.file) )

# Fanta Orange 2mm Model Modelloptimierung

# set wd and source packages and functions
# PLS and LM ####
fao$pls <- lapply(fao$seqp$model, function(x) pls_function(csv_transfered = x
                                                           , substance = dt$para$substance
                                                           , wlr = fao$plspara$wl 
                                                           , ncomp = fao$plspara$ncomp
                                                           , spc = "spc"))


fao$pls_lm <- mapply(function( x , y ) pls_lm_function(x
                                                       , csv_transfered = y
                                                       , substance = dt$para$substance
                                                       , wlr = fao$plspara$wl 
                                                       , ncomp = fao$plspara$ncomp)
                     , x = fao$pls
                     , y = fao$seqp$model
                     , SIMPLIFY = F)

# Prediction ####
gc()
memory.limit(99999)

fao$pred <- lapply(fao$pls, function(y) lapply(fao$trs$spc,function(x) produktion_prediction(csv_transfered = x, pls_function_obj = y, ncomp = fao$plspara$ncomp)))
fao$pred.lin <- lapply(fao$pls, function(y) produktion_prediction(csv_transfered = fao$lin, pls_function_obj = y, ncomp = fao$plspara$ncomp))

# Find the best model ####
for(i in 1:length(fao$pls)) fao$pls_merge[[i]] <- lapply(fao$pred[[i]], function(x) .merge_pls(pls_pred = x, pls_lm = fao$pls_lm[[i]] ,R2=.5, mean = c(75,125)))
for(i in 1:length(fao$pls)) fao$pls_merge_site[[i]] <- .merge_pls_site(merge_pls_lm_predict_ls = fao$pls_merge[[i]][ c(1,2) ], number = 500, ncomp = fao$plspara$ncomp)

lapply( fao$pls_merge_site, nrow)
lapply( fao$pls_merge_site, head )

for(i in 1:length(fao$pls)) fao$pls_merge_site[[i]] <- fao$pls_merge[[i]][[1]]


fao$pred.lin.diff <- mapply( function( x , y ) linearitaet_filter( linearitaet_prediction =  x$prediction
                                                                   , ncomp = fao$plspara$ncomp
                                                                   , linearitaet_limit_1 = 4.465 * .5
                                                                   , linearitaet_limit_2 = 5.456 * 1.5
                                                                   , R_2 = .8
                                                                   , SOLL = fao$lin$data$Acid
                                                                   , pls_merge = y )
                             , x = fao$pred.lin
                             , y = fao$pls_merge_site
                             , SIMPLIFY = F
)
fao$pred.lin.diff
lapply(fao$pred.lin.diff, head)
fao$pred.lin.diff[[3]][ order( fao$pred.lin.diff[[4]]$spcR), ]

# fao$pred.lin.diff[[1]] <- NULL

lapply(fao$pred.lin.diff, function(x) head(x[order(x$sd),]))
lapply(fao$pred.lin.diff, function(x) head(x[order(x$mad),]))


dat1 <- fao$pls_merge[[3]]$Mannheim_MY
head( dat1[ order(dat1$mad) , ] )
head( dat1[ order(dat1$sd) , ] )
# calculate best model ####
mod_c <- list()
mod_c$ncomp <-4
mod_c$wl1 <- 425
mod_c$wl2 <- 465
mod_c$wl3 <- NA
mod_c$wl4 <- NA
mod_c$spc <- "spc"
mod_c$matrix <- 1

names(fao$seqp$Probe)[mod_c$matrix]
mod_c$pngname <- paste(mod_c$spc, mod_c$wl1, mod_c$wl2, mod_c$wl3, mod_c$wl4, "PC", mod_c$ncomp, "matrix", mod_c$matrix, sep = "_")
mod_c$modelname <- paste(mod_c$wl1, mod_c$wl2, mod_c$wl3, mod_c$wl4, mod_c$spc, mod_c$ncomp, sep = "_")

mod_c$model <- pls_function(csv_transfered = fao$seqp$model[[ mod_c$matrix ]]
                            , substance = dt$para$substance
                            , wlr = data.frame(mod_c$wl1, mod_c$wl2, mod_c$wl3, mod_c$wl4)
                            , ncomp = mod_c$ncomp
                            , spc = mod_c$spc
                            , validation = "none")

mod_c$pred$lin <- pred_of_new_model(modell_csv_transfered = fao$seqp$model[[ mod_c$matrix ]]
                                    , substance = dt$para$substance
                                    , wl1 = mod_c$wl1
                                    , wl2 = mod_c$wl2
                                    , wl3 = mod_c$wl3
                                    , wl4 = mod_c$wl4
                                    , ncomp = mod_c$ncomp
                                    , derivative = mod_c$spc
                                    , csv_transfered = fao$lin)

setwd(dt$wd)
setwd("./Modelloptimierung")
setwd(paste0("./", dt$para$mop.date, "_", dt$para$model.pl[1], "_", dt$para$para))
setwd("./Analyse")

png(paste0(.date(),"_FaO_Prediction_2mm_Linearitaet.png"),xxx<-4800,xxx/16*9,"px",12,"white",res=500,"sans",T,"cairo")
par(mfrow = c(1,1), mar = c(5,4,5,5))

x.bias <- .bias(median(fao$lin$data$Acid),0,median(mod_c$pred$lin, na.rm = T))
plot(mod_c$pred$lin + x.bias
     , axes = T, xlab = "", ylab = "Acid in %"
     , main = paste("Linearitaet von Fanta Orange in Dorsten DS, Pfadlaenge umgerechnet von 0,31 auf 2 mm")
     , pch = 20, cex = 2
     , ylim = c(90, 100)
     , sub = paste("Bias =", round(x.bias,2)))
points(fao$lin$data$Acid, cex = 2, col = "red", pch = 18)
sd(mod_c$pred$lin + x.bias - fao$lin$data$Acid)

legend("bottomright", c("LG3", "Labor"), pch = c(20, 18), col = c("black", "red"), cex = 1.5)
dev.off()

mod_c$pred$spc <- list()
for(i in 1:length(fao$trs$spc))
  mod_c$pred$spc[[i]] <- pred_of_new_model(modell_csv_transfered = fao$seqp$model[[ mod_c$matrix ]]
                                           , substance = dt$para$substance
                                           , wl1 = mod_c$wl1
                                           , wl2 = mod_c$wl2
                                           , wl3 = mod_c$wl3
                                           , wl4 = mod_c$wl4
                                           , ncomp = mod_c$ncomp
                                           , derivative = mod_c$spc
                                           , csv_transfered = fao$trs$spc[[i]])

for(i in 1:length(mod_c$pred$spc)) mod_c$pred$spc[[i]] <- as.numeric(ma(mod_c$pred$spc[[i]], 5))

png(paste0(.date(),"_FaO_Prediction_2mm.png"),xxx<-4800,xxx/16*9,"px",12,"white",res=500,"sans",T,"cairo")
par(mar = c(4,4,3,4))

layout(matrix(c(1,2,3,4,5,5), byrow = T, ncol = 2), heights = c(1, 1, .125))

for(i in 1:length(fao$trs$spc)){
  x.bias <- .bias(median(mod_c$pred$spc[[i]], na.rm = T),  0, 100)
  plot(mod_c$pred$spc[[i]] - x.bias
       , axes = F, xlab = "", ylab = "Acid in %"
       , main = paste(fao$para$beverage[i])
       , pch = 20, cex = .5
       , ylim = c(95, 105)
       , sub = paste("Bias =", round(x.bias,2)))
  .xaxisdate(fao$trs$spc[[i]]$data$datetime)
  
  abline( h = 100 + c(4.5 / 299 * 100) * c( 1, -1) , col = "orange", lwd = 1)
  abline( h = 100 + c(7.6 / 299 * 100) * c( 1, -1) , col = "red", lwd = 1)
  
  par(new = T)
  plot(fao$trs$spc[[i]]$data$brix, pch = 20, col = "blue", axes = F, xlab = "", ylab = "", ylim = c(8 * .9, 8 * 1.1), cex = .5)
  axis(4)
}

par(mar = c(0,0,0,0))
plot(1,1,axes=F,xlab="",ylab="", type = "n")
legend("center", c("LG3", "Brix", "Eingriffsgrenze", "Sperrgrenze")
       , pch = c(20, 20, NA, NA)
       , lty = c(NA, NA, 1, 1)
       , col = c("black", "blue", "orange", "red"), ncol = 4, cex = 1.5, bty = "n", lwd = 2)
dev.off()

.keep.out.unsb(fao$model, mod_c$wl1, mod_c$wl2, mod_c$wl3, mod_c$wl4)

# final model matrix ####
setwd(dt$wd)
setwd("./Modelloptimierung")
setwd(paste0("./", dt$para$mop.date, "_", dt$para$model.pl[1], "_", dt$para$para))
setwd("./Modellmatrix")

fao$seqp$export_matrix <- data.frame(cbind(fao$seqp$model[[ mod_c$matrix ]]$data, fao$seqp$model[[ mod_c$matrix ]]$spc))

write.csv2(fao$seqp$export_matrix
           , paste0(.date(), "_", fao$para$bev,"_",mod_c$pngname,"_mop_modelmatrix_final.csv")
           , row.names = F)


