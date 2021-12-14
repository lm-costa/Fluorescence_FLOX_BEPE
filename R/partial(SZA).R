nfile <- list.files("data-raw/partial/", pattern = "F_alf")
F_Spec <- readr::read_csv(paste0("data-raw/partial/",nfile))

rfile <- list.files("data-raw/partial/", pattern = "R_alf")
R_Spec <- readr::read_csv(paste0("data-raw/partial/",nfile))

ind_tab <- readr::read_csv("data-raw/partial/all_INDEX_table_alfalfa.csv")
metri <- readr::read_csv("data-raw/partial/table_metrics_alfalfa.csv")
metri <- metri[,-c(1:2)]



spect <- list.files("data-raw/clear/", pattern = "_F_")
spect <- readr::read_rds(paste0("data-raw/clear/",spect))
spect <- round(as.numeric(names(spect)[2:683]),2)



F_time <- t(F_Spec)
F_time <- data.frame(F_time)
names(F_time) <- spect

R_Spec <- data.frame(t(R_Spec))
names(R_Spec) <- spect
Rnir <- R_Spec$`770`


tab <- cbind(ind_tab, metri,Rnir, F_time)

tab <- tab |>
  dplyr::mutate(
    day = DOYdayfrac%/%1
  ) |>
  dplyr::filter(
    day != 144 & day!=151  & day!= 156 & day!=158 & day!=166
    & day!=168 & day!=192  # criteria i
    & day!= 147 & day != 148 & day!=149  & day!=157 & day!=163
    & day!=172 & day!=175 & day!=176 & day!=179 & day!=184 & day!=189
    & day!=193 & day!= 199 & day!=200 & day!=201
    & day!=161 & day!=162 & day!=164 & day!=165 & day!=167 & day!=169
    & day!= 170 & day!=171 & day!=173 & day!=174 & day!=177 & day!=178 & day!=185
  ) |>
  dplyr::filter(
    lubridate::hour(UTC_datetime)>=9 & lubridate::hour(UTC_datetime) <=16
  ) |>
  dplyr::filter(SZA <50 & SZA >20)

min_int<-"30 min" #step
sttime <- tab$UTC_datetime[1]
endtime <- tab$UTC_datetime[(nrow(tab))]


timetoagg<-seq.POSIXt(from = as.POSIXct(sttime,tz = "UTC"),
                      to = as.POSIXct(endtime,tz = "UTC"),
                      by = min_int)

time_agg <- aggregate(x = tab,
                      by = list(cut(as.POSIXct(tab$UTC_datetime,tz = "UTC"),
                                    timetoagg)), FUN = mean)


gppFile <- list.files("data-raw/clear/", pattern = "aa_")
gpp_table <- readr::read_csv(paste0("data-raw/clear/", gppFile))


time_agg <- time_agg |> dplyr::mutate(Hr1 = lubridate::hour(as.character(Group.1)),
                                      Mn1 = lubridate::minute(as.character(Group.1)))

gpp_table <- gpp_table |> dplyr::mutate(day=DOY%/%1)

comp_f_gpp <- dplyr::semi_join(time_agg,gpp_table)
comp_gpp_f <- dplyr::semi_join(gpp_table,time_agg)
comp_gpp_f <- comp_gpp_f[,-c(6,7,176)]

df <- cbind(comp_f_gpp,comp_gpp_f)


spect_tab <- df[,-c(1:85,768:943)]
spect_tab <- as.data.frame(t(spect_tab))

c_spec <- as.character(df$Group.1)
colnames(spect_tab) <- c_spec
spect <- names(df)[86:767]
spect_tab <- cbind(spect,spect_tab)

integral <- function(x,y){
  approxfun(x,y)
}

result <- vector("list", length = ncol(spect_tab))

for (i in 2:(ncol(spect_tab))){
  int <- integral(spect_tab$spect, spect_tab[,i])
  result[[i]] <- integrate(int, as.numeric(min(spect_tab$spect)),
                           as.numeric(max(spect_tab$spect)))
}

for(i in 2:length(result)){
  if(i==2){
    fint <- result[[i]]$value
  }else{
    fint_a <- result[[i]]$value
    fint <- rbind(fint,fint_a)
    fint <- as.data.frame(fint)
  }
}

row.names(fint) <- seq(1:(nrow(fint)))
colnames(fint) <- "Fint"

## Now we are gonna add the integral in the final dataset

df <- cbind(df,fint)

######
fint_t <- vector("numeric", length = ncol(spect_tab))

fint_t[2:ncol(spect_tab)] <- t(fint)

spect_tab[683, ] <- fint_t

gpp <- vector("numeric", length = ncol(spect_tab))
gpp[2:ncol(spect_tab)] <- df$GPP_DT_U95
spect_tab[684, ] <- gpp

for(i in 2:ncol(spect_tab)){
  if(i==2){
    tab2 <- as.numeric(spect_tab[,i])
  }else{
    tab2_a <- as.numeric(spect_tab[,i])
    tab2 <- cbind(tab2,tab2_a)
  }
}
tab2 <- data.frame(tab2)
colnames(tab2) <- df$Group.1
tab2 <- data.frame(t(tab2))
colnames(tab2) <- spect_tab$spect
colnames(tab2)[683] <- "fint"
colnames(tab2)[684] <- "gpp"

##
pear_resul <- vector("list", length = ncol(tab2))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(tab2$gpp, tab2[,i])
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

rownames(correl) <- colnames(tab2)
correl[,4] <- spect_tab$spect

df <- df |>
  dplyr::mutate(
    k1=584*sin(circular::rad(SZA)) - 88,
    k2= 3.55424 - 1.15937*cos(circular::rad(SZA)),
    aPAR= k1*sin(k2*Rnir),
    Fy= Fint/aPAR,
    LUE=GPP_DT_U95/aPAR
  )

pear_resul <- vector("list", length = ncol(tab2))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(df$LUE, tab2[,i]/df$aPAR)
}

for(i in 1:length(pear_resul)){
  if(i==1){
    correl_fy <- pear_resul[[i]]$estimate
    lower <- pear_resul[[i]]$conf.int[1]
    uper <- pear_resul[[i]]$conf.int[2]
    correl_fy <- cbind(correl_fy, lower,uper)
  }else{
    correl_fya <- pear_resul[[i]]$estimate
    lower_a <- pear_resul[[i]]$conf.int[1]
    uper_a <- pear_resul[[i]]$conf.int[2]
    correl_fya <- cbind(correl_fya,lower_a,uper_a)
    correl_fy <- rbind(correl_fy,correl_fya)
  }
}

correl_fy <- data.frame(correl_fy)

rownames(correl_fy) <- colnames(tab2)
correl_fy[,4] <- spect_tab$spect




fobs <- tab[,-c(767)]

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


for(i in 85:length(fobs)){
  if(i==85){
    Frc <- a1*exp(-(((as.numeric(names(fobs[i]))-b1)/c1)^2)) + a2*exp(-(((as.numeric(names(fobs[i]))-b2)/c2)^2))
  }else{
    frc <- a1*exp(-(((as.numeric(names(fobs[i]))-b1)/c1)^2)) + a2*exp(-(((as.numeric(names(fobs[i]))-b2)/c2)^2))
    Frc <- cbind(Frc, frc)
  }
}

Frc <- data.frame(Frc)
names(Frc) <- round(as.numeric(names(fobs)[(85:length(fobs))]),2)

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
spectra <- names(Frc[0, c(3:684)])
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
    LUE = GPP_DT_U95/df$aPAR,
    fy = fint/df$aPAR
  )

fint_t <- vector("numeric", length = ncol(spectra_Table))

fint_t[2:ncol(spectra_Table)] <- t(fint)

spectra_Table[683, ] <- fint_t

gpp <- vector("numeric", length = ncol(spectra_Table))
gpp[2:ncol(spectra_Table)] <- df$GPP_DT_U95
spectra_Table[684, ] <- gpp

for(i in 2:ncol(spectra_Table)){
  if(i==2){
    tab2 <- as.numeric(spectra_Table[,i])
  }else{
    tab2_a <- as.numeric(spectra_Table[,i])
    tab2 <- cbind(tab2,tab2_a)
  }
}
tab2 <- data.frame(tab2)
colnames(tab2) <- df$Group.1
tab2 <- data.frame(t(tab2))
colnames(tab2) <- spectra_Table$spect
colnames(tab2)[683] <- "fint"
colnames(tab2)[684] <- "gpp"


pear_resul <- vector("list", length = ncol(tab2))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(tab2$gpp, tab2[,i])
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

rownames(correl_fr) <- colnames(tab2)
correl_fr[,4] <- spectra_Table$spectra

pear_resul <- vector("list", length = ncol(tab2))
for(i in 1:length(pear_resul)){
  pear_resul[[i]] <- cor.test(df$LUE, tab2[,i]/df$aPAR)
}

for(i in 1:length(pear_resul)){
  if(i==1){
    correl_fyc <- pear_resul[[i]]$estimate
    lower <- pear_resul[[i]]$conf.int[1]
    uper <- pear_resul[[i]]$conf.int[2]
    correl_fyc <- cbind(correl_fyc, lower,uper)
  }else{
    correl_fyac <- pear_resul[[i]]$estimate
    lower_a <- pear_resul[[i]]$conf.int[1]
    uper_a <- pear_resul[[i]]$conf.int[2]
    correl_fyac <- cbind(correl_fyac,lower_a,uper_a)
    correl_fyc <- rbind(correl_fyc,correl_fyac)
  }
}

correl_fyc <- data.frame(correl_fyc)

rownames(correl_fyc) <- colnames(tab2)
correl_fyc[,4] <- spect_tab$spect




correl[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(V4), y=correl, ymax=uper, ymin= lower
  ))+
  ggplot2::geom_ribbon(alpha=.1, color="black",linetype= "dashed")+
  ggplot2::geom_point(color="black")+
  ggplot2::geom_ribbon(alpha=.1, color="red",linetype= "dashed",ggplot2::aes(ymin=correl_fr$lower[-c(683,684)],
                                                                             ymax=correl_fr$uper[-c(683,684)]))+
  ggplot2::geom_point(ggplot2::aes(x=correl_fr$V4[-c(683,684)],
                                   y=correl_fr$correl_fr[-c(683,684)]),color="red")+
  ggplot2::ylim(-0.1,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (GPP vs F('~lambda~'))' ))+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=14,color='black'))




correl_fy[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(V4), y=correl_fy, ymax=uper, ymin= lower
  ))+
  ggplot2::geom_ribbon(alpha=.1, color="black",linetype= "dashed")+
  ggplot2::geom_point(color="black")+
  ggplot2::geom_ribbon(alpha=.1, color="red",linetype= "dashed",ggplot2::aes(ymin=correl_fyc$lower[-c(683,684)],
                                                                             ymax=correl_fyc$uper[-c(683,684)]))+
  ggplot2::geom_point(ggplot2::aes(x=as.numeric(correl_fyc$V4[-c(683,684)]),
                                   y=correl_fyc$correl_fyc[-c(683,684)]),color="red")+
  ggplot2::ylim(.999,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (LUE vs Fy ('~lambda~'))' ))+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=14,color='black'))


df_frc|>
  tidyr::pivot_longer(
    cols = "669.98":"779.96",
    names_to = "wavelength",
    values_to = "Fluorescence") |>
  dplyr::filter(day!=155) |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(wavelength), y=Fluorescence, color=lubridate::hour(UTC))) +
  ggplot2::geom_point() +
  ggplot2::geom_line()+
  ggplot2::facet_wrap(~day)+
  ggplot2::labs(x="Wavelength",color="Hour")+
  ggplot2::labs(x=expression("Wavelength ("~ lambda~")"), color="Hour")+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12))
