df <- readr::read_rds("data/clear/df.rds")


df <- df |>
  dplyr::select(GPP_DT_U95,Fint,`669.98`:`779.96`)

df_funtion <- function(data.f,name="") {
  for(i in 1:(ncol(data.f))){
    if(i==1){
      m <- mean(data.f[,i])
    }else{
      ma <- mean(data.f[,i])
      m <- rbind(m,ma)
    }
  }
  m <- data.frame(m)
  names(m) <- name
  row.names(m) <- seq(1:nrow(m))
  return(m)
}


fred <- df |>
  dplyr::select(`683.04`:`687.13`) |>
  t() |>
  data.frame() |>
  df_funtion(name = "Fred")


fre <- df |>
  dplyr::select(`705.1`:`720.07`)|>
  t() |>
  data.frame() |>
  df_funtion(name = "Fred-edge")

ffr <- df |>
  dplyr::select(`742.08`:`748.07`)|>
  t() |>
  data.frame() |>
  df_funtion(name="Ffar-red")

f757 <- df |>
  dplyr::select(`755.09`:`759.12`) |>
  t() |>
  data.frame() |>
  df_funtion(name = "F757")

f771 <- df |>
  dplyr::select(`769.09`:`772.13`)|>
  t() |>
  data.frame() |>
  df_funtion(name="F771")


df <- data.frame(
  GPP=df$GPP_DT_U95,
  Fint=df$Fint,
  Fred=fred,
  FFr=ffr,
  F771=f771,
  F757=f757,
  Fredeg=fre
)

df[,-1] <- df[,-1] |> scale()

set.seed(123)

inTrain <- caret::createDataPartition(df$GPP,p=0.6,list = FALSE)
train <- df[inTrain,]
teste <- df[-inTrain,]

RandomForest <- randomForest::randomForest(GPP~.,
                                 data=train,
                                 ntree=500,
                                 replace=T,
                                 nodesize=10,
                                 maxnode=15)

RandomForest

par(mfrow=c(1,2))
plot(RandomForest,cex.axis=1.3, cex.lab=1.5,main = "")
randomForest::varImpPlot(RandomForest, main="", cex.lab=1.5, cex.axis=1.3)

pred_1 <- predict(RandomForest,newdata = teste)

cv <- data.frame(est=pred_1,
                 obs=teste$GPP)

rmse <- Metrics::rmse(cv$obs,cv$est)
mape <- Metrics::mape(cv$obs,cv$est)*100

cv |> ggplot2::ggplot(ggplot2::aes(x=obs,y=est))+
  ggplot2::geom_point()+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::annotate('text',x=16.4,y=26.7,label=paste0('RMSE = ',round(rmse,2),', MAPE = '
                                                    ,round(mape,2),'%'),size=5)+
  ggplot2::labs(x=expression(paste('GPP'[OBS]~mu*mol,' ',m^-2,' ',s^-1)),
                y=expression(paste('GPP'[EST]~mu*mol,' ',m^-2,' ',s^-1)),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



####
fitControl <- caret::trainControl(method = "repeatedcv",
                                  number=3,
                                  repeats = 500)

step <- caret::train(GPP~., data=train,
                         method= "leapSeq",
                         tuneGrid=data.frame(nvmax=1:6),
                         trControl=fitControl
)

summary(step)
step |>
  ggplot()+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )


pred_2 <- predict(step,newdata = teste)

cv_2 <- data.frame(obs=teste$GPP,est=pred_2)

rmse <- Metrics::rmse(cv_2$obs,cv_2$est)
mape <- Metrics::mape(cv_2$obs,cv_2$est)*100

cv_2 |> ggplot2::ggplot(ggplot2::aes(x=obs,y=est))+
  ggplot2::geom_point()+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")),size=5)+
  ggplot2::annotate('text',x=16.4,y=27.7,label=paste0('RMSE = ',round(rmse,2),', MAPE = '
                                                      ,round(mape,2),'%'),size=5)+
  ggplot2::labs(x=expression(paste('GPP'[OBS]~mu*mol,' ',m^-2,' ',s^-1)),
                y=expression(paste('GPP'[EST]~mu*mol,' ',m^-2,' ',s^-1)),
                col="DOY")+
  ggplot2::theme_bw()+
  ggplot2::theme(axis.title= ggplot2::element_text(size=14),
                 axis.text = ggplot2::element_text(size=12,color='black'),
                 legend.text = ggplot2::element_text(size=12)
  )



pred_3 <- predict(step,newdata = df)
cv_3 <- data.frame(obs=df$GPP,est=pred_3)
cv_3 |> ggplot2::ggplot(ggplot2::aes(x=obs,y=est))+
  ggplot2::geom_point()+
  ggplot2::geom_smooth(method = "lm")+
  ggpubr::stat_regline_equation(ggplot2::aes(
    label =  paste(..eq.label.., ..rr.label.., sep = "*plain(\",\")~~")))
Metrics::rmse(cv_3$obs,cv_3$est)
