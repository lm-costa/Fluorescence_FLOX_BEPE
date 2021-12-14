correl_partial <- readr::read_rds("data/partial/correl.rds")
correl_frc_partial <- readr::read_rds("data/partial/correl_frc.rds")
correl_fy_partial <- readr::read_rds("data/partial/correl_fy.rds")
correl_fyc_partial <- readr::read_rds("data/partial/correl_fyc.rds")
df_partial <- readr::read_rds("data/partial/df.rds")
df_frc_partial <- readr::read_rds("data/partial/df_frc.rds")


correl_partial[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4,y=correl,ymax=uper,ymin=lower))+
  ggplot2::geom_ribbon(alpha=.1,linetype="dashed",color="black")+
  ggplot2::geom_jitter(color="black")+
  ggplot2::geom_ribbon(ggplot2::aes(ymax=correl_frc_partial$uper[-c(683,684)],
                                    ymin=correl_frc_partial$lower[-c(683,684)]),
                       color='red', linetype='dashed',alpha=.1)+
  ggplot2::geom_point(ggplot2::aes(x=correl_frc_partial$V4[-c(683,684)],
                                    y=correl_frc_partial$correl_fr[-c(683,684)]),
                      color='red')+
  ggplot2::ylim(0,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (GPP vs F('~lambda~'))' ))+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=14,color='black'))



correl_fy_partial[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4,y=correl_fy,ymax=uper,ymin=lower))+
  ggplot2::geom_ribbon(alpha=.1,linetype="dashed",color="black")+
  ggplot2::geom_jitter(color="black")+
  ggplot2::geom_ribbon(ggplot2::aes(ymax=correl_fyc_partial$uper[-c(683,684)],
                                    ymin=correl_fyc_partial$lower[-c(683,684)]),
                       color='red', linetype='dashed',alpha=.1)+
  ggplot2::geom_point(ggplot2::aes(x=correl_fyc_partial$V4[-c(683,684)],
                                   y=correl_fyc_partial$correl_fyc[-c(683,684)]),
                      color='red')+
  ggplot2::ylim(.9,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (LUE vs Fy ('~lambda~'))' ))+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=14,color='black'))

df_frc_partial |>
  tidyr::pivot_longer(
    cols = "669.98":"779.96",
    names_to = "wavelength",
    values_to = "Fluorescence") |>
  dplyr::filter(day!=155) |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(wavelength), y=Fluorescence, color=lubridate::hour(UTC))) +
  ggplot2::geom_point() +
  ggplot2::geom_line()+
  ggplot2::facet_wrap(~day)+
  ggplot2::labs(x=expression("Wavelength ("~ lambda~")"), color="Hour")+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12))

df_partial |>
  tidyr::pivot_longer(
    cols = "669.98":"779.96",
    names_to = "wavelength",
    values_to = "Fluorescence") |>
  dplyr::filter(day!=155) |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(wavelength), y=Fluorescence, color=lubridate::hour(UTC_datetime))) +
  ggplot2::geom_point() +
  ggplot2::geom_line()+
  ggplot2::facet_wrap(~day)+
  ggplot2::xlab(label="Wavelength")+
  ggplot2::labs(col='Hour', fontsize=14)+
  ggplot2::labs(x=expression("Wavelength ("~ lambda~")"), color="Hour")+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12))

df_partial |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=`685.09`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[RED]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('GPP ('~mu*mol,' ',m^-2,' ',s^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_partial |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=`740.02`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[FAR-RED]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('GPP ('~mu*mol,' ',m^-2,' ',s^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_partial |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=Fint))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[Int]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('GPP ('~mu*mol,' ',m^-2,' ',s^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_partial |>
  ggplot2::ggplot(ggplot2::aes(y=df_frc_partial$Fint,x=Fint))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[Int-Obs]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('F'[Int-RC]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_partial |>
  ggplot2::ggplot(ggplot2::aes(y=df_frc_partial$`685.09`,x=`685.09`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[Red-Obs]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('F'[Red-RC]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_partial |>
  ggplot2::ggplot(ggplot2::aes(y=df_frc_partial$`740.02`,x=`740.02`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[FR-Obs]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('F'[FR-RC]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )


df_frc_partial |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=`685.09`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression('F'[Red-RC]~Wm^-2*sr^-1*mu*m^-1),
                y=expression(paste('GPP'~mu*mol,' ',m^-2,' ',s^-1)),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )


df_partial |>
  ggplot2::ggplot(ggplot2::aes(x=LUE,y=Fy))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x='LUE', y='Fy', col='DOY')+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_frc_partial |>
  ggplot2::ggplot(ggplot2::aes(x=LUE,y=Fint/df_partial$aPAR))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x='LUE', y=expression('Fy '[RC]), col='DOY')+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



correl_clear <- readr::read_rds("data/clear/correl.rds")
correl_frc_clear <- readr::read_rds("data/clear/correl_frc.rds")
correl_fy_clear <- readr::read_rds("data/clear/correl_fy.rds")
correl_fyc_clear <- readr::read_rds("data/clear/correl_fyc.rds")
df_clear <- readr::read_rds("data/clear/df.rds")
df_frc_clear <- readr::read_rds("data/clear/df_frc.rds")


correl_clear[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4,y=correl,ymax=uper,ymin=lower))+
  ggplot2::geom_ribbon(alpha=.1,linetype="dashed",color="black")+
  ggplot2::geom_jitter(color="black")+
  ggplot2::geom_ribbon(ggplot2::aes(ymax=correl_frc_clear$uper[-c(683,684)],
                                    ymin=correl_frc_clear$lower[-c(683,684)]),
                       color='red', linetype='dashed',alpha=.1)+
  ggplot2::geom_point(ggplot2::aes(x=correl_frc_clear$V4[-c(683,684)],
                                   y=correl_frc_clear$correl_fr[-c(683,684)]),
                      color='red')+
  ggplot2::ylim(0,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (GPP vs F('~lambda~'))' ))+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=14,color='black'))



correl_fy_clear[-c(683,684),] |>
  ggplot2::ggplot(ggplot2::aes(x=V4,y=correl_fy,ymax=uper,ymin=lower))+
  ggplot2::geom_ribbon(alpha=.1,linetype="dashed",color="black")+
  ggplot2::geom_jitter(color="black")+
  ggplot2::geom_ribbon(ggplot2::aes(ymax=correl_fyc_clear$uper[-c(683,684)],
                                    ymin=correl_fyc_clear$lower[-c(683,684)]),
                       color='red', linetype='dashed',alpha=.1)+
  ggplot2::geom_point(ggplot2::aes(x=correl_fyc_clear$V4[-c(683,684)],
                                   y=correl_fyc_clear$correl_fyc[-c(683,684)]),
                      color='red')+
  ggplot2::ylim(.89,1)+
  ggplot2::labs(x=expression("Wavelenght ("~ lambda ~")"),
                y=expression('Correlation (LUE vs Fy ('~lambda~'))' ))+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=14,color='black'))

df_frc_clear |>
  tidyr::pivot_longer(
    cols = "669.98":"779.96",
    names_to = "wavelength",
    values_to = "Fluorescence") |>
  dplyr::filter(day!=155) |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(wavelength), y=Fluorescence, color=lubridate::hour(UTC))) +
  ggplot2::geom_point() +
  ggplot2::geom_line()+
  ggplot2::facet_wrap(~day)+
  ggplot2::labs(x=expression("Wavelength ("~ lambda~")"), color="Hour")+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12))

df_clear |>
  tidyr::pivot_longer(
    cols = "669.98":"779.96",
    names_to = "wavelength",
    values_to = "Fluorescence") |>
  dplyr::filter(day!=155) |>
  ggplot2::ggplot(ggplot2::aes(x=as.numeric(wavelength), y=Fluorescence, color=lubridate::hour(UTC_datetime))) +
  ggplot2::geom_point() +
  ggplot2::geom_line()+
  ggplot2::facet_wrap(~day)+
  ggplot2::xlab(label="Wavelength")+
  ggplot2::labs(col='Hour', fontsize=14)+
  ggplot2::labs(x=expression("Wavelength ("~ lambda~")"), color="Hour")+
  ggplot2::theme(axis.title= ggplot2::element_text(size=16),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12))

df_clear |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=`685.09`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[RED]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('GPP ('~mu*mol,' ',m^-2,' ',s^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_clear |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=`740.02`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[FAR-RED]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('GPP ('~mu*mol,' ',m^-2,' ',s^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_clear |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=Fint))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[Int]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('GPP ('~mu*mol,' ',m^-2,' ',s^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_clear |>
  ggplot2::ggplot(ggplot2::aes(y=df_frc_clear$Fint,x=Fint))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[Int-Obs]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('F'[Int-RC]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_clear |>
  ggplot2::ggplot(ggplot2::aes(y=df_frc_clear$`685.09`,x=`685.09`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[Red-Obs]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('F'[Red-RC]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_clear |>
  ggplot2::ggplot(ggplot2::aes(y=df_frc_clear$`740.02`,x=`740.02`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression(paste('F'[FR-Obs]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                y=expression(paste('F'[FR-RC]~' (',Wm^-2*sr^-1*mu*m^-1,')')),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )


df_frc_clear |>
  ggplot2::ggplot(ggplot2::aes(y=GPP_DT_U95,x=`685.09`))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x=expression('F'[Red-RC]~Wm^-2*sr^-1*mu*m^-1),
                y=expression(paste('GPP'~mu*mol,' ',m^-2,' ',s^-1)),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )


df_clear |>
  ggplot2::ggplot(ggplot2::aes(x=LUE,y=Fy))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x='LUE', y='Fy', col='DOY')+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



df_frc_clear |>
  ggplot2::ggplot(ggplot2::aes(x=LUE,y=Fint/df_clear$aPAR))+
  ggplot2::geom_jitter(ggplot2::aes(color=as.factor(day)))+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::labs(x='LUE', y=expression('Fy '[RC]), col='DOY')+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )
