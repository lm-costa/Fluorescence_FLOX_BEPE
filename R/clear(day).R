namefile <- list.files("data-raw/clear/", pattern = "_F_")
ori_tab <- readr::read_rds(paste0("data-raw/clear/",namefile))
names(ori_tab) <- c("DOY",
                    round(
                      as.numeric(names(ori_tab)[2:length(ori_tab)]), 2)
)

ind_file <- readr::read_csv("data-raw/clear/table_all_index.csv") |>
  dplyr::filter(DOYdayfrac > 127.9)

time_table <- cbind(ind_file,ori_tab)
time_table <- time_table |>
  dplyr::filter(DOY >= 144 & DOY <=193) |>
  dplyr::filter(lubridate::hour(UTC_datetime)>8) |>
  dplyr::mutate(day = DOY%/%1) |>
  dplyr::filter(day == 181)


min_int<-"30 min" #step
sttime <- time_table$UTC_datetime[1]
endtime <- time_table$UTC_datetime[(nrow(time_table))]

timetoagg<-seq.POSIXt(from = as.POSIXct(sttime,tz = "UTC"),
                      to = as.POSIXct(endtime,tz = "UTC"),
                      by = min_int)

time_agg <- aggregate(x = time_table,
                      by = list(cut(as.POSIXct(time_table$UTC_datetime,tz = "UTC"),
                                    timetoagg)), FUN = mean)


gppFile <- list.files("data-raw/clear/", pattern = "aa_")
gpp_table <- readr::read_csv(paste0("data-raw/clear/", gppFile))

time_agg <- time_agg |> dplyr::mutate(Hr1 = lubridate::hour(as.character(Group.1)),
                                      Mn1 = lubridate::minute(as.character(Group.1)))

time_agg <- time_agg[,-7]

gpp_table <- gpp_table |> dplyr::mutate(day=DOY%/%1)

comp_f_gpp <- dplyr::semi_join(time_agg,gpp_table)
comp_gpp_f <- dplyr::semi_join(gpp_table,time_agg)
comp_gpp_f <- comp_gpp_f[,-c(6,7,176)]

df <- cbind(comp_f_gpp,comp_gpp_f)


spectra_Table <- df[,-c(689:(ncol(df)))]
spectra_Table <- as.data.frame(t(spectra_Table))
c_spectra <- paste0(lubridate::hour(as.character(spectra_Table[1,])),":",
                    lubridate::minute(as.character(spectra_Table[1, ])))
colnames(spectra_Table) <- c_spectra
spectra_Table <- spectra_Table[-c(1:6), ]
rownames(spectra_Table) <- seq(1:682)
spectra <- names(ori_tab[0, c(2:683)])
spectra_Table <- cbind(spectra,spectra_Table)
colnames(spectra_Table) <- make.unique(names(spectra_Table))
spectra_Table[,1] <-as.numeric(spectra_Table[,1])

integral <- function(x,y){
  approxfun(x,y)
}

result <- vector("list", length = ncol(spectra_Table))

for (i in 1:(ncol(spectra_Table))){
  int <- integral(spectra_Table$spectra, spectra_Table[,i])
  result[[i]] <- integrate(int, min(spectra_Table$spectra),
                           max(spectra_Table$spectra))
}

for(i in 1:length(result)){
  if(i==1){
    fint <- result[[i]]$value
  }else{
    fint_a <- result[[i]]$value
    fint <- rbind(fint,fint_a)
    fint <- as.data.frame(fint)
  }
}

row.names(fint) <- seq(1:nrow(fint))
fint <- as.data.frame(fint[-1,])
colnames(fint) <- "Fint"

df <- cbind(df, fint)

fint_t <- vector("numeric", length = ncol(spectra_Table))
fint_t[2:(ncol(spectra_Table))] <- t(fint)

spectra_Table[683, ] <- fint_t
gpp <- vector("numeric", length = ncol(spectra_Table))
gpp[2:ncol(spectra_Table)] <- comp_gpp_f$GPP_DT_U95
spectra_Table[684, ] <- gpp

for(i in 2:ncol(spectra_Table)){
  if(i==2){
    tab <- as.numeric(spectra_Table[,i])
  }else{
    tab_a <- as.numeric(spectra_Table[,i])
    tab <- cbind(tab,tab_a)
  }
}
tab <- as.data.frame(tab)
colnames(tab) <- comp_gpp_f$TIMESTAMP
tab <- as.data.frame(t(tab))
colnames(tab) <- spectra_Table$spectra
colnames(tab)[683] <- "fint"
colnames(tab)[684] <- "gpp"

pear_resul <- vector("list", length = ncol(tab))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(tab$gpp, tab[,i])
}

for(i in 1:length(pear_resul)){
  if(i==1){
    correl <- pear_resul[[i]]$estimate
    lower <- pear_resul[[i]]$conf.int[1]
    uper <- pear_resul[[i]]$conf.int[2]
    correl <- cbind(correl, lower,uper)
  }else{
    correl_a <- pear_resul[[i]]$estimate
    lower_a <- pear_resul[[i]]$conf.int[1]
    uper_a <- pear_resul[[i]]$conf.int[2]
    correl_a <- cbind(correl_a,lower_a,uper_a)
    correl <- rbind(correl,correl_a)
  }
}

correl <- data.frame(correl)

rownames(correl) <- colnames(tab)
correl[,4] <- spectra_Table$spectra

correl[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4, y= correl, ymax= uper, ymin= lower))+
  ggplot2::geom_ribbon(alpha=.1,linetype="dashed",color="black")+
  ggplot2::geom_jitter(color="black")+
  ggplot2::ylim(0,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (LUE vs Fy ('~lambda~'))' ))

df |>
  tidyr::pivot_longer(
    cols = "669.98":"779.96",
    names_to = "wavelength",
    values_to = "Fluorescence") |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(wavelength), y=Fluorescence, color=lubridate::hour(UTC_datetime))) +
  ggplot2::geom_point() +
  ggplot2::geom_line()+
  ggplot2::facet_wrap(~day)+
  ggplot2::labs(x="Wavelength", color='Hour')

rfile <- list.files("data-raw/clear/", pattern = "_R_")
rtab <- readr::read_rds(paste0("data-raw/clear/",rfile))

rspec <- as.data.frame(t(rtab))
c_name <- rspec[1,]
colnames(rspec) <- c_name
rspec <- rspec[-1, ]

rnir <- rspec[c(616:682),]
rnir_spec <- as.numeric(rownames(rnir))
rnir <- cbind(rnir_spec,rnir)

rred <- rspec[c(1:49),]
rred_spec <- as.numeric(rownames(rred))
rred <- cbind(rred_spec,rred)

redge <- rspec[c(207:267),]
redge_spec <- as.numeric(rownames(redge))
redge <- cbind(redge_spec,redge)

rnir_result <- vector("list", length = ncol(rnir))
rred_result <- vector("list", length = ncol(rred))
redge_result <- vector("list", length = ncol(redge))

for (i in 1:(ncol(rnir))){
  int <- integral(rnir$rnir_spec, rnir[,i])
  rnir_result[[i]] <- integrate(int, min(rnir$rnir_spec), max(rnir$rnir_spec))
}

for(i in 1:length(rnir_result)){
  if(i==1){
    rnir_int <- rnir_result[[i]]$value
  }else{
    rnir_int_a <- rnir_result[[i]]$value
    rnir_int <- rbind(rnir_int,rnir_int_a)
    rnir_int <- as.data.frame(rnir_int)
  }
}

rnir_int <- rnir_int[-1,]
rnir_int <- as.data.frame(rnir_int)

for (i in 1:(ncol(rred))){
  int <- integral(rred$rred_spec, rred[,i])
  rred_result[[i]] <- integrate(int, min(rred$rred_spec), max(rred$rred_spec))
}

for(i in 1:length(rred_result)){
  if(i==1){
    rred_int <- rred_result[[i]]$value
  }else{
    rred_int_a <- rred_result[[i]]$value
    rred_int <- rbind(rred_int,rred_int_a)
    rred_int <- as.data.frame(rred_int)
  }
}
rred_int <- rred_int[-1,]
rred_int <- as.data.frame(rred_int)

for (i in 1:(ncol(redge))){
  int <- integral(redge$redge_spec, redge[,i])
  redge_result[[i]] <- integrate(int, min(redge$redge_spec), max(redge$redge_spec))
}

for(i in 1:length(redge_result)){
  if(i==1){
    redge_int <- redge_result[[i]]$value
  }else{
    redge_int_a <- redge_result[[i]]$value
    redge_int <- rbind(redge_int,redge_int_a)
    redge_int <- as.data.frame(redge_int)
  }
}

redge_int <- redge_int[-1,]
redge_int <- as.data.frame(redge_int)

ndvi <- (rnir_int - rred_int)/(rnir_int+rred_int)
colnames(ndvi) <- "NDVI"

ndvi_r <- (rnir_int- redge_int)/(rnir_int + redge_int)
colnames(ndvi_r) <- "NDVI Red"

rtab <- cbind(rtab,ndvi,ndvi_r)

rtab <- cbind(ind_file,rtab)
rtab <- rtab |>
  dplyr::filter(DOY >= 144 & DOY <=193) |>
  dplyr::filter(lubridate::hour(UTC_datetime)>8) |>
  dplyr::mutate(day = DOY%/%1) |>
  dplyr::filter(day == 181)

ndvi_agg <- aggregate(x = rtab,
                      by = list(cut(as.POSIXct(rtab$UTC_datetime,tz = "UTC"),
                                    timetoagg)), FUN = mean)

ndvi_agg <- ndvi_agg |> dplyr::mutate(Hr1 = lubridate::hour(as.character(Group.1)),
                                      Mn1 = lubridate::minute(as.character(Group.1)))

ndvi_agg <- ndvi_agg[,-7]

comp_ndvi_gpp <- dplyr::semi_join(ndvi_agg,gpp_table)

t_rnir <- data.frame(t(rnir))
names(t_rnir) <- round(as.numeric(t_rnir[1,]),2)
t_rnir <- t_rnir[-1,]

t_rnir <- cbind(ind_file,t_rnir)

t_rnir <- t_rnir |>
  dplyr::filter(DOYdayfrac >= 144 & DOYdayfrac <=193) |>
  dplyr::filter(lubridate::hour(UTC_datetime)>8) |>
  dplyr::mutate(day = DOYdayfrac%/%1) |>
  dplyr::filter(day == 181)


timetoagg<-seq.POSIXt(from = as.POSIXct(sttime,tz = "UTC"),
                      to = as.POSIXct(endtime,tz = "UTC"),
                      by = min_int)

nir_agg <- aggregate(x = t_rnir,
                     by = list(cut(as.POSIXct(t_rnir$UTC_datetime,tz = "UTC"),
                                   timetoagg)), FUN = mean)
nir_agg <- nir_agg |>
  dplyr::select(
    Group.1,
    `775`
  )
names(nir_agg)[2] <- "R775"

nir_agg <- dplyr::semi_join(nir_agg,comp_ndvi_gpp)
NDVI <- comp_ndvi_gpp$NDVI
NDVI_R <- comp_ndvi_gpp$`NDVI Red`
NIR = NDVI*nir_agg$R775
ndvi_df <- data.frame(NDVI,NDVI_R, NIR)

df <- df |>
  dplyr::mutate(
    k1=584*sin(circular::rad(SZA)) - 88,
    k2= 3.55424 - 1.15937*cos(circular::rad(SZA)),
    aPAR= k1*sin(k2*NIR),
    Fy= Fint/aPAR,
    LUE=GPP_DT_U95/aPAR
  )


pear_resul <- vector("list", length = ncol(tab))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(tab$gpp/df$aPAR, tab[,i]/df$aPAR)
}

for(i in 1:length(pear_resul)){
  if(i==1){
    correl_fy <- pear_resul[[i]]$estimate
    lower <- pear_resul[[i]]$conf.int[1]
    uper <- pear_resul[[i]]$conf.int[2]
    correl_fy <- cbind(correl_fy, lower,uper)
  }else{
    correl_afy <- pear_resul[[i]]$estimate
    lower_a <- pear_resul[[i]]$conf.int[1]
    uper_a <- pear_resul[[i]]$conf.int[2]
    correl_afy <- cbind(correl_afy,lower_a,uper_a)
    correl_fy <- rbind(correl_fy,correl_afy)
  }
}

correl_fy <- data.frame(correl_fy)

rownames(correl_fy) <- colnames(tab)
correl_fy[,4] <- spectra_Table$spectra

correl_fy[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4, y= correl_fy, ymax= uper, ymin= lower))+
  ggplot2::geom_point()+
  ggplot2::geom_ribbon(alpha=.2)


###

fobs <- time_table

a = dplyr::case_when(
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(22))~ 81.36 - 13.45*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(24)) ~ 78.61 -10.49*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(26)) ~ 77.05 -8.78*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(28)) ~ 77 -8.727*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(30)) ~ 75.157 -6.637*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(32)) ~ 73.526 -4.753*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(34)) ~ 72.087 -3.057*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(36)) ~ 70.665 -1.342*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(38)) ~ 68.9615 +0.7643*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(40)) ~ 67.24 +2.949*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(42)) ~ 65.552 +5.152*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(44)) ~ 63.412 +8.033*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) <=cos(circular::rad(46)) ~ 61.39 +10.84*cos(circular::rad(fobs$SZA)),
  abs(cos(circular::rad(fobs$SZA))) > cos(circular::rad(46)) ~ 59.58 +13.45*cos(circular::rad(fobs$SZA)),
)

b = -0.23*cos(circular::rad(fobs$SZA))^2 + 0.58*cos(circular::rad(fobs$SZA)) + 0.43


a1 = a*fobs$`760.05`^b
a2 = 0.2546*a1
b1 = 682.2
b2=706.7
c1=11.49
c2=54.47


for(i in 7:688){
  if(i==7){
    Frc <- a1*exp(-(((as.numeric(names(fobs[i]))-b1)/c1)^2)) + a2*exp(-(((as.numeric(names(fobs[i]))-b2)/c2)^2))
  }else{
    frc <- a1*exp(-(((as.numeric(names(fobs[i]))-b1)/c1)^2)) + a2*exp(-(((as.numeric(names(fobs[i]))-b2)/c2)^2))
    Frc <- cbind(Frc, frc)
  }
}

Frc <- data.frame(Frc)
names(Frc) <- round(as.numeric(names(fobs)[(7:688)]),2)

Frc <- cbind(fobs$DOYdayfrac, fobs$UTC_datetime, Frc)
names(Frc)[c(1,2)] <- c("DOY", "UTC")
Frc <- Frc |> na.omit()


time_agg_frc <- aggregate(x = Frc,
                          by = list(cut(as.POSIXct(Frc$UTC,tz = "UTC"),
                                        timetoagg)), FUN = mean)

time_agg_frc <- time_agg_frc |> dplyr::mutate(day = time_agg_frc$DOY%/%1,
                                              Hr1 = lubridate::hour(as.character(Group.1)),
                                              Mn1 = lubridate::minute(as.character(Group.1)))

time_agg_frc <- time_agg_frc[,-2]


comp_frc_gpp <- dplyr::semi_join(time_agg_frc,gpp_table)
comp_gpp_frc <- dplyr::semi_join(gpp_table,time_agg_frc)
comp_gpp_frc <- comp_gpp_frc[,-c(6,7,176)]

df_frc <- cbind(comp_frc_gpp,comp_gpp_frc)


spectra_Table <- df_frc[,-c(685:(ncol(df_frc)))]
spectra_Table <- as.data.frame(t(spectra_Table))
c_spectra <- paste0(lubridate::hour(as.character(spectra_Table[1,])),":",
                    lubridate::minute(as.character(spectra_Table[1, ])))
colnames(spectra_Table) <- c_spectra
spectra_Table <- spectra_Table[-c(1:2), ]
rownames(spectra_Table) <- seq(1:682)
spectra <- names(ori_tab[0, c(2:683)])
spectra_Table <- cbind(spectra,spectra_Table)
colnames(spectra_Table) <- make.unique(names(spectra_Table))
spectra_Table[,1] <-as.numeric(spectra_Table[,1])

result <- vector("list", length = ncol(spectra_Table))

for (i in 1:(ncol(spectra_Table))){
  int <- integral(spectra_Table$spectra, spectra_Table[,i])
  result[[i]] <- integrate(int, min(spectra_Table$spectra),
                           max(spectra_Table$spectra))
}

for(i in 1:length(result)){
  if(i==1){
    fint <- result[[i]]$value
  }else{
    fint_a <- result[[i]]$value
    fint <- rbind(fint,fint_a)
    fint <- as.data.frame(fint)
  }
}

row.names(fint) <- seq(1:nrow(fint))
fint <- as.data.frame(fint[-1,])
colnames(fint) <- "Fint"

df_frc <- cbind(df_frc, fint)

df_frc <- df_frc |>
  dplyr::mutate(
    LUE = df_frc$GPP_DT_U95/df$aPAR,
    fy = fint/df$aPAR
  )


fint_t <- vector("numeric", length = ncol(spectra_Table))
fint_t[2:(ncol(spectra_Table))] <- t(fint)

spectra_Table[683, ] <- fint_t
gpp <- vector("numeric", length = ncol(spectra_Table))
gpp[2:ncol(spectra_Table)] <- comp_gpp_f$GPP_DT_U95
spectra_Table[684, ] <- gpp

for(i in 2:ncol(spectra_Table)){
  if(i==2){
    tab <- as.numeric(spectra_Table[,i])
  }else{
    tab_a <- as.numeric(spectra_Table[,i])
    tab <- cbind(tab,tab_a)
  }
}
tab <- as.data.frame(tab)
colnames(tab) <- comp_gpp_f$TIMESTAMP
tab <- as.data.frame(t(tab))
colnames(tab) <- spectra_Table$spectra
colnames(tab)[683] <- "fint"
colnames(tab)[684] <- "gpp"

pear_resul <- vector("list", length = ncol(tab))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(tab$gpp, tab[,i])
}

for(i in 1:length(pear_resul)){
  if(i==1){
    correl_fr <- round(pear_resul[[i]]$estimate,2)
    lower <- pear_resul[[i]]$conf.int[1]
    uper <- pear_resul[[i]]$conf.int[2]
    correl_fr <- cbind(correl_fr, lower,uper)
  }else{
    correl_afr <- round(pear_resul[[i]]$estimate,2)
    lower_a <- pear_resul[[i]]$conf.int[1]
    uper_a <- pear_resul[[i]]$conf.int[2]
    correl_afr <- cbind(correl_afr,lower_a,uper_a)
    correl_fr <- rbind(correl_fr,correl_afr)
  }
}

correl_fr <- data.frame(correl_fr)

rownames(correl_fr) <- colnames(tab)
correl_fr[,4] <- spectra_Table$spectra



correl[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4, y=correl, ymax=uper, ymin= lower
  ))+
  ggplot2::geom_ribbon(alpha=.1, color="black",linetype= "dashed")+
  ggplot2::geom_point(color="black")+
  ggplot2::geom_ribbon(alpha=.1, color="red",linetype= "dashed",ggplot2::aes(ymin=correl_fr$lower[-c(683,684)],
                                                                             ymax=correl_fr$uper[-c(683,684)]))+
  ggplot2::geom_point(ggplot2::aes(x=correl_fr$V4[-c(683,684)],
                                   y=correl_fr$correl_fr[-c(683,684)]),color="red")+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (GPP vs F ('~lambda~'))' ))


#######


pear_resul <- vector("list", length = ncol(tab))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(df_frc$LUE, tab[,i]/df$aPAR)
}

for(i in 1:length(pear_resul)){
  if(i==1){
    correl_fyc <- round(pear_resul[[i]]$estimate,2)
    lower <- pear_resul[[i]]$conf.int[1]
    uper <- pear_resul[[i]]$conf.int[2]
    correl_fyc <- cbind(correl_fyc, lower,uper)
  }else{
    correl_afyc <- round(pear_resul[[i]]$estimate,2)
    lower_a <- pear_resul[[i]]$conf.int[1]
    uper_a <- pear_resul[[i]]$conf.int[2]
    correl_afyc <- cbind(correl_afyc,lower_a,uper_a)
    correl_fyc <- rbind(correl_fyc,correl_afyc)
  }
}

correl_fyc <- data.frame(correl_fyc)

rownames(correl_fyc) <- colnames(tab)
correl_fyc[,4] <- spectra_Table$spectra



correl[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4, y=correl, ymax=uper, ymin= lower
  ))+
  ggplot2::geom_ribbon(alpha=.1, color="black",linetype= "dashed")+
  ggplot2::geom_point(color="black")+
  ggplot2::geom_ribbon(alpha=.1, color="red",linetype= "dashed",ggplot2::aes(ymin=correl_fr$lower[-c(683,684)],
                                                                             ymax=correl_fr$uper[-c(683,684)]))+
  ggplot2::geom_point(ggplot2::aes(x=correl_fr$V4[-c(683,684)],
                                   y=correl_fr$correl_fr[-c(683,684)]),color="red")+
  ggplot2::ylim(0,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (GPP vs F('~lambda~'))' ))+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=14,color='black'))
