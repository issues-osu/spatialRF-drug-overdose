library(tigris)
library(sf)
library(dplyr)
library(readr)
library(GWmodel)      ## GW models
library(plyr)         ## Data management
library(sp)           ## Spatial Data management
library(spdep)        ## Spatial autocorrelation
library(RColorBrewer) ## Visualization
library(classInt)     ## Class intervals
library(raster)       ## spatial data
library(grid)         ## plot
library(gridExtra)    ## Multiple plot
library(ggplot2)      #  plotting
library(tidyverse)    # data 
library(SpatialML) 
library(pdp)       # for partial dependence plots (and ICE curves)
library(vip)       # for variable importance plots
library(fmsb)
library(scales)

setwd("C:/Users/barboza-salerno.1/Downloads/spatialRF-opioid/")
#sink("model_output_log.txt")

atlantic_states <- c(
  "NY", 
  "NJ", 
  "NC", 
  "SC", 
  "GA", 
  "DE", 
  "MD", 
  "PA", 
  "VA",
  "WV"
  )

options(tigris_use_cache = TRUE)
counties_sf <- counties(cb = TRUE, year = 2022) %>%
  st_as_sf() %>%
  filter(STUSPS %in% atlantic_states)

cdc_data <- read_csv("cdc-drug-deaths-atlantic-all.csv") %>%
  dplyr::rename(FIPS = `County Code`) %>%
  mutate(FIPS = sprintf("%05s", FIPS))  # Ensure FIPS is character with 5 digits

merged <- counties_sf %>%
  left_join(cdc_data, by = c("GEOID" = "FIPS"))

merged_5070 <- st_transform(merged, crs = 5070)

merged_5070 <- merged_5070 %>%
  mutate(centroid = st_centroid(geometry),
         X = st_coordinates(centroid)[,1],
         Y = st_coordinates(centroid)[,2])

head(merged_5070 %>% dplyr::select(GEOID, NAME, STUSPS, X, Y, everything()))

df <- merged_5070 %>% dplyr::select(GEOID, X, Y, `Crude Rate`, Year)

#st_write(df, "atlantic_counties_drug_deaths_5070.gpkg")

wide_df <- df %>%
  dplyr::rename(Rate = `Crude Rate`) %>%
  st_drop_geometry() %>%
  pivot_wider(names_from = Year, values_from = Rate, names_prefix = "Rate_")

# Correlation matrix between year columns
(cor_matrix <- wide_df %>%
  dplyr::select(starts_with("Rate_")) %>%
  as.data.frame() %>%
  cor(use = "pairwise.complete.obs"))

ejiindicators <- read.csv("eji-indicators.csv",  colClasses = c("GEOID" = "character"))

df <- df %>% left_join(ejiindicators)

county<- as_Spatial(counties_sf)
state.bd<-shapefile("STATE_ATLANTIC.shp")

state <- list("sp.lines", as(state.bd, "SpatialLines"), col="grey50", lwd=.7,lty=3) 

df$Rate <- log(df$`Crude Rate`)
test.df<-df %>% 
  dplyr::select(GEOID, X, Y, Year, UNEMP, UNINSUR, DISABL, PARK, MHLTH, MOBILE, PM, TRI, WLKIND, Rate ) %>%
  filter(Year == 2022) %>% st_drop_geometry()

valid.df<-df %>% 
  dplyr::select(GEOID, X, Y, Year, UNEMP, UNINSUR, DISABL, PARK, MHLTH, MOBILE, PM, TRI, WLKIND, Rate ) %>%
  filter(Year == 2020) %>% st_drop_geometry()

train.df<-df %>% 
  dplyr::select(GEOID, X, Y, Year, UNEMP, UNINSUR, DISABL, PARK, MHLTH, MOBILE, PM, TRI, WLKIND, Rate) %>%
  filter(Year == 2018) %>% st_drop_geometry()


test.df[, 5:13] = scale(test.df[, 5:13])
valid.df[, 5:13] = scale(valid.df[, 5:13])
train.df[, 5:13] = scale(train.df[, 5:13])


library(h2o)
h2o.init(nthreads = -1,max_mem_size ="48g",enable_assertions = FALSE)

test.mf<-test.df[, 5:14] 
valid.mf<-valid.df[, 5:14]
train.mf<-train.df[, 5:14]

test.hex<-  as.h2o(test.mf)
valid.hex<-  as.h2o(valid.mf)
train.hex<-  as.h2o(train.mf)

response <- "Rate"
predictors <- setdiff(names(train.hex), c("Rate", response))

# Hyper-parameter
#drf_hyper_params <-list(
#  ntrees  = seq(10, 5000, by = 10),
##  max_depth=c(10,20,30,40,50),
##  sample_rate=c(0.7, 0.8, 0.9, 1.0)
#)

drf_hyper_params <- list(
  ntrees = c(100, 250, 500, 750, 1000),  # Strategic points
  max_depth = c(10, 20, 30, 40),
  sample_rate = c(0.7, 0.8, 0.9, 1.0)
)


#  serach criteria
drf_search_criteria <- list(
  strategy = "RandomDiscrete", 
  max_models = 200,
  max_runtime_secs = 900,
  stopping_tolerance = 0.001,
  stopping_rounds = 2,
  seed = 1345767
  )

# Grid Search
drf_grid <- h2o.grid(
  algorithm="randomForest",
  grid_id = "drf_grid_IDx",
  x= predictors,
  y = response,
  training_frame = train.hex,
  validation_frame = valid.hex,
  stopping_metric = "RMSE",
  nfolds=10,
  keep_cross_validation_predictions = TRUE,
  hyper_params = drf_hyper_params,
  search_criteria = drf_search_criteria,
  seed = 42)


# RF Grid parameters
drf_get_grid <- h2o.getGrid("drf_grid_IDx",sort_by="RMSE",decreasing=FALSE)
drf_get_grid@summary_table[1,]

best_drf <- h2o.getModel(drf_get_grid@model_ids[[1]]) 
#capture.output(print(summary(best_drf)),file =  "DRF_summary_N_RY.txt")
best_drf

cv.drf<-best_drf@model$cross_validation_metrics_summary%>%.[,c(1,2)]
cv.drf

FIPS.xy<-test.df[,1:3]
FIPS.xy$Obs_2011<-valid.df[,14]
FIPS.xy$Obs_2012<-test.df[,14]
# validation data
pred.valid<-as.data.frame(h2o.predict(object = best_drf, newdata = valid.hex))

FIPS.xy$RF_2011<-pred.valid$predict
# test data
pred.test<-as.data.frame(h2o.predict(object = best_drf, newdata = test.hex))
FIPS.xy$RF_2012<-pred.test$predict

cat('RF Validation RMSE:', round(sqrt(mean((FIPS.xy$RF_2011-FIPS.xy$Obs_2011)^2 , na.rm = TRUE)), digits=3), '\n')
cat('RF Validation MAE:', round(mean(abs(FIPS.xy$RF_2011-FIPS.xy$Obs_2011) , na.rm = TRUE ), digits=3), '\n')
cat('RF Validation R2:', round(summary(lm(Obs_2011~RF_2011,FIPS.xy))$r.squared, digits=3), '\n')
cat('RF Test RMSE:', round(sqrt(mean((FIPS.xy$RF_2012-FIPS.xy$Obs_2012)^2 , na.rm = TRUE)), digits=3), '\n')
cat('RF Test MAE:', round(mean(abs(FIPS.xy$RF_2012-FIPS.xy$Obs_2012) , na.rm = TRUE ), digits=3), '\n')
cat('RF Test R2:', round(summary(lm(Obs_2012~RF_2012,FIPS.xy))$r.squared, digits=3), '\n')

shared_x_range <- range(
  c(FIPS.xy$Obs_2011, FIPS.xy$RF_2011, FIPS.xy$Obs_2012, FIPS.xy$RF_2012),
  na.rm = TRUE
)

shared_y_range <- range(
  c(FIPS.xy$Obs_2011, FIPS.xy$RF_2011, FIPS.xy$Obs_2012, FIPS.xy$RF_2012),
  na.rm = TRUE
)


p.valid<-ggplot(data=FIPS.xy, aes(x=Obs_2011, y=RF_2011))+ 
  geom_point(size = 1.0)+
  geom_smooth(method = "lm", se = FALSE, colour="black",size=0.5)+
  scale_x_continuous(limits = shared_x_range, breaks = pretty(shared_x_range)) +
  scale_y_continuous(limits = shared_y_range, breaks = pretty(shared_y_range)) +
   theme_light() +
  ggtitle("Validation data (2020-2021)")+
  theme(
    plot.title = element_text(color="black", size=12,hjust = 0.5),
    panel.background = element_rect(fill = "white",colour = "gray75",size = 0.5, linetype = "solid"),
    axis.line = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size=10, colour="black"),
    axis.text.y=element_text(size=10,angle = 90,vjust = 0.5, hjust=0.5, colour='black'))+
  geom_abline(slope=1, intercept=0,linetype="dashed",size=0.5)+
  labs(x="Observed", y = "Predicted") 

p.test<-ggplot(data=FIPS.xy, aes(x=Obs_2012, y=RF_2012))+ 
  geom_point(size = 1.0)+
  geom_smooth(method = "lm", se = FALSE, colour="black",size=0.5)+
  scale_x_continuous(limits = shared_x_range, breaks = pretty(shared_x_range)) +
  scale_y_continuous(limits = shared_y_range, breaks = pretty(shared_y_range)) +
  theme_light() +
  ggtitle("Test data (2018-2019)")+
  theme(
    plot.title = element_text(color="black", size=12,hjust = 0.5),
    panel.background = element_rect(fill = "white",colour = "gray75",size = 0.5, linetype = "solid"),
    axis.line = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size=10, colour="black"),
    axis.text.y=element_text(size=10,angle = 90,vjust = 0.5, hjust=0.5, colour='black'))+
  geom_abline(slope=1, intercept=0,linetype="dashed",size=0.5)+
  labs(x="Observed", y = "Predicted")

grid.arrange(p.valid,p.test, 
             ncol= 2, 
             heights = c(30,6), 
             top = textGrob("Random Forest Predicted Drug-Related Mortality Rate\n(death per 100,000)",gp=gpar(fontsize=18)))

rf.SPDF<-merge(county, FIPS.xy, by.y= "GEOID")
names(rf.SPDF)

my.breaks <- seq(5.5, 7.8, by = 0.1)
myPalette <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(length(my.breaks) - 1)

spplot(rf.SPDF, 
       c("RF_2011", "RF_2012", "Obs_2011", "Obs_2012" ),         
       layout=c(2, 2),
       par.strip.text=list(cex=1, lines=1, col="black"),
       main = list( "Random Forest Predicted Drug-Related Mortality Rate\n(death per 100,000)", adj = 0.5),
       sp.layout=list(state),
       par.settings=list(strip.background=list(col="transparent"),
                         strip.border = list(col = 'grey95',lwd=0.5),
                         axis.line=list(col="grey95",lwd=0.5)),
       col="transparent",
       colorkey=list(space="right",height=1, width=2,labels=list(cex=.75)),
       col.regions=myPalette)


features <- as.data.frame(train.hex) %>%  dplyr::select(-Rate)

response <- as.data.frame(train.hex) %>% pull(Rate)

pred <- function(object, newdata)  {
  results <- as.vector(h2o.predict(object, as.h2o(newdata)))
  return(results)
}

per.var.imp<-vip(
  best_drf,
  train = as.data.frame(train.hex),
  method = "permute",
  target = "Rate",
  metric = "RMSE",
  nsim = 5,
  sample_frac = 0.5,
  pred_wrapper = pred
)
per.var.imp<-per.var.imp$data
# Scale the importance
per.var.imp$Rank<-(per.var.imp$Importance-min(per.var.imp$Importance))*100/(max(per.var.imp$Importance)-min(per.var.imp$Importance))
head(per.var.imp)
ggplot(per.var.imp, aes(y=Rank, x=reorder(Variable, +Rank))) +
  ylab('Scale-Importance') +
  xlab('')+
  geom_bar(stat="identity", width = 0.45,fill = "darkblue")+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  coord_flip()+
  ggtitle("Permutation-based Feature Importance")

# Custom prediction function wrapper
pdp_pred <- function(object, newdata)  {
  results <- mean(as.vector(h2o.predict(object, as.h2o(newdata))))
  return(results)
}

var.names=per.var.imp[1]
var_1st=var.names[1,]$Variable
var_1st

var_2nd=var.names[2,]$Variable
var_2nd

pd_var_1st<- partial(
  best_drf,
  train = as.data.frame(train.hex),
  pred.var = var_1st,
  pred.fun = pdp_pred,
  parallel = TRUE,
  grid.resolution = 10
)

pd_1<-autoplot(pd_var_1st,
               rug = TRUE,
               train=as.data.frame(train.hex)) +
  theme(text = element_text(size=15))

pd_var_2nd<- partial(
  best_drf,
  train = as.data.frame(train.hex),
  pred.var = var_2nd,
  pred.fun = pdp_pred,
  parallel = TRUE,
  grid.resolution = 10
)

pd_2<-autoplot(pd_var_2nd,
               rug = TRUE,
               train=as.data.frame(train.hex)) +
  theme(text = element_text(size=15))

pd_var_1st_2nd<- partial(
  best_drf,
  train = as.data.frame(train.hex),
  pred.var = c(var_1st,var_2nd),
  pred.fun = pdp_pred,
  parallel = TRUE,
  grid.resolution = 10
)

pd_1_2<-autoplot(pd_var_1st_2nd,
                 contour = TRUE) +
  theme(text = element_text(size=15))

grid.arrange(
  pd_1,pd_2,pd_1_2, 
  ncol= 3,
  heights = c(40,8), 
  top = textGrob("Partial Dependence Plot", gp=gpar(fontsize=25))
  )

top5_vars <- per.var.imp$Variable[1:5]

pd_plots <- list()

for (i in 1:5) {
  pd_obj <- partial(
    object = best_drf,
    train = as.data.frame(train.hex),
    pred.var = top5_vars[i],
    pred.fun = pdp_pred,
    parallel = TRUE,
    grid.resolution = 10
  )
  
  pd_plot <- autoplot(pd_obj, rug = TRUE, train = as.data.frame(train.hex)) +
    ggtitle(paste0("Partial Dependence: ", top5_vars[i])) +
    theme_minimal(base_size = 14)
  
  pd_plots[[i]] <- pd_plot
}

grid.arrange(
   grobs = pd_plots,
   ncol = 3,
   top = textGrob(
     "Top 5 Predictor Partial Dependence Plots", 
     gp = gpar(
       fontsize = 20, fontface = "bold"
       )
     )
   )


Coords<-train.df[ ,2:3]
complete_idx <- complete.cases(train.df[, c("Rate", "UNEMP", "UNINSUR", "DISABL", "PARK", "MHLTH")])

train_clean <- train.df[complete_idx, ]
coords_clean <- Coords[complete_idx, ]
coords_clean <- as.matrix(coords_clean)

grf.model <- grf(Rate ~ UNEMP+UNINSUR+DISABL+PARK+MHLTH+MOBILE+PM+TRI+WLKIND, 
                 dframe=train_clean, 
                 bw=162,             
                 # a positive number, in the case of an "adaptive kernel" or a real in the case of a "fixed kernel".
                 ntree=   500, 
                 mtry = 2,           
                 # n integer referring to the number of trees to grow for each of the local random forests.
                 kernel="adaptive",  
                 # yhe kernel to be used in the regression. Options are "adaptive" or "fixed".
                 forests = TRUE,     
                 coords=Coords)      

grf.model$Global.Model$variable.importance 
mean(grf.model$Global.Model$prediction.error)
grf.model$Global.Model$r.squared 
grf.model$LocalModelSummary

counties_sf$incMSE.SMOK=grf.model$Local.Variable.Importance$UNEMP
counties_sf$incMSE.POV=grf.model$Local.Variable.Importance$UNINSUR
counties_sf$incMSE.PM25=grf.model$Local.Variable.Importance$DISABL
counties_sf$incMSE.NO2=grf.model$Local.Variable.Importance$PARK
counties_sf$incMSE.MOB=grf.model$Local.Variable.Importance$MOBILE
counties_sf$incMSE.PM=grf.model$Local.Variable.Importance$PM
counties_sf$incMSE.TRI=grf.model$Local.Variable.Importance$TRI
counties_sf$incMSE.WALK=grf.model$Local.Variable.Importance$WLKIND


gwrf.SPDF <- as_Spatial(counties_sf)
col.palette.t<-colorRampPalette(c("blue",  "sky blue", "green","yellow","pink", "red"),space="rgb",interpolate = "linear") 
state_proj <- spTransform(state[[2]], CRSobj = CRS(proj4string(gwrf.SPDF)))
state <- list("sp.lines", state_proj, col = "black", lwd = 1, lty = 3)

gwrf.SPDF@data$loc_R2=grf.model$LGofFit$LM_Rsq100
myPaletteRes <- colorRampPalette(c("lightseagreen","lightsteelblue1", "moccasin","hotpink", "red"))
r2 <- spplot(gwrf.SPDF,"loc_R2", main = "Local R2 (%)", 
       sp.layout=list(state),
       col="transparent",
       col.regions=myPaletteRes(100))

smok<-spplot(gwrf.SPDF,"incMSE.SMOK", main = "Unemployment", 
             sp.layout=list(state),
             col="transparent",
             col.regions=rev(col.palette.t(100)))

pov<-spplot(gwrf.SPDF,"incMSE.POV", main = "Uninsurance", 
            sp.layout=list(state),
            col="transparent",
            col.regions=rev(col.palette.t(100)))

pm25<-spplot(gwrf.SPDF,"incMSE.PM25", main = "Disability", 
             sp.layout=list(state),
             col="transparent",
             col.regions=rev(col.palette.t(100)))

no2<-spplot(gwrf.SPDF,"incMSE.NO2", main = "Park", 
            sp.layout=list(state),
            col="transparent",
            col.regions=rev(col.palette.t(100)))

so2<-spplot(gwrf.SPDF,"incMSE.NO2", main = "Mental Health", 
            sp.layout=list(state),
            col="transparent",
            col.regions=rev(col.palette.t(100)))

mo<-spplot(gwrf.SPDF,"incMSE.MOB", main = "Mobile Homes", 
            sp.layout=list(state),
            col="transparent",
            col.regions=rev(col.palette.t(100)))

tri<-spplot(gwrf.SPDF,"incMSE.TRI", main = "TRI", 
            sp.layout=list(state),
            col="transparent",
            col.regions=rev(col.palette.t(100)))

pm2<-spplot(gwrf.SPDF,"incMSE.PM", main = "PM", 
            sp.layout=list(state),
            col="transparent",
            col.regions=rev(col.palette.t(100)))

walk<-spplot(gwrf.SPDF,"incMSE.WALK", main = "Walkability", 
            sp.layout=list(state),
            col="transparent",
            col.regions=rev(col.palette.t(100)))

grid.arrange(smok, pov, pm25, no2, so2, mo, tri, pm2, walk, r2,
             ncol = 5, nrow = 2,
             heights = rep(1, 2),  
             top = textGrob("Local Feature Importance (IncMSE)", gp = gpar(fontsize = 25)))

# Validation data
FIPS.xy$GWRF_2011<-predict.grf(grf.model, valid.df, x.var.name="X", y.var.name="Y", local.w=1, global.w=0)
# Test data
FIPS.xy$GWRF_2012<-predict.grf(grf.model, test.df, x.var.name="X", y.var.name="Y", local.w=1, global.w=0)

cat('GWRF Validation RMSE:', round(sqrt(mean((FIPS.xy$GWRF_2011-FIPS.xy$Obs_2011)^2 , na.rm = TRUE)), digits=3), '\n')
cat('GWRF Validation MAE:', round(mean(abs(FIPS.xy$GWRF_2011-FIPS.xy$Obs_2011) , na.rm = TRUE ), digits=3), '\n')
cat('GWRF Validation R2:', round(summary(lm(Obs_2011~GWRF_2011,FIPS.xy))$r.squared, digits=3), '\n')
cat('GWRF Test RMSE:', round(sqrt(mean((FIPS.xy$GWRF_2012-FIPS.xy$Obs_2012)^2 , na.rm = TRUE)), digits=3), '\n')
cat('GWRF Test MAE:', round(mean(abs(FIPS.xy$GWRF_2012-FIPS.xy$Obs_2012) , na.rm = TRUE ), digits=3), '\n')
cat('GWRF Test R2:', round(summary(lm(Obs_2012~GWRF_2012,FIPS.xy))$r.squared, digits=3), '\n')

p1.valid<-ggplot(data=FIPS.xy, aes(x=Obs_2011, y=GWRF_2011))+ 
  geom_point(size = 1.0)+
  geom_smooth(method = "lm", se = FALSE, colour="black",size=0.5)+
  #scale_x_continuous(limits=c(20,120), breaks=seq(20, 120, 20))+ 
 # scale_y_continuous(limits=c(20,120), breaks=seq(20, 120, 20)) +
  theme_light() +
  ggtitle("Validation data (2020-2021)")+
  theme(
    plot.title = element_text(color="black", size=12,hjust = 0.5),
    panel.background = element_rect(fill = "white",colour = "gray75",size = 0.5, linetype = "solid"),
    axis.line = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size=10, colour="black"),
    axis.text.y=element_text(size=10,angle = 90,vjust = 0.5, hjust=0.5, colour='black'))+
  geom_abline(slope=1, intercept=0,linetype="dashed",size=0.5)+
  labs(x="Observed", y = "Predicted") 

p1.test<-ggplot(data=FIPS.xy, aes(x=Obs_2012, y=GWRF_2012))+ 
  geom_point(size = 1.0)+
  geom_smooth(method = "lm", se = FALSE, colour="black",size=0.5)+
  #scale_x_continuous(limits=c(20,120), breaks=seq(20, 120, 20))+ 
  #scale_y_continuous(limits=c(20,120), breaks=seq(20, 120, 20)) +
  theme_light() +
  ggtitle("Test data (2018-2019)")+
  theme(
    plot.title = element_text(color="black", size=12,hjust = 0.5),
    panel.background = element_rect(fill = "white",colour = "gray75",size = 0.5, linetype = "solid"),
    axis.line = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size=10, colour="black"),
    axis.text.y=element_text(size=10,angle = 90,vjust = 0.5, hjust=0.5, colour='black'))+
  geom_abline(slope=1, intercept=0,linetype="dashed",size=0.5)+
  labs(x="Observed", y = "Predicted")

grid.arrange(p1.valid,p1.test, 
             ncol= 2, 
             heights = c(30,6), 
             top = textGrob("GW Random Forest Predicted Substance-Related Mortality Rate\n(death per 100,000)",gp=gpar(fontsize=18)))

gwrf.SPDF<-merge(gwrf.SPDF, FIPS.xy,by = "GEOID")
names(gwrf.SPDF)


#state <- list("sp.lines", as(state, "SpatialLines"), col="grey50", lwd=.7,lty=3) 
my.breaks <- seq(5.5, 7.8, by = 0.1)
myPalette <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(length(my.breaks) - 1)

spplot(gwrf.SPDF, 
       c("GWRF_2011", "GWRF_2012", "Obs_2011", "Obs_2012"),         
       layout = c(2, 2),
       par.strip.text = list(cex = 1, lines = 1, col = "black"),
       main = list("GW Random Forest Predicted Substance-Related Mortality Rate\n(death per 100,000)", adj = 0.5),
       sp.layout = list(state),
       par.settings = list(
         strip.background = list(col = "transparent"),
         strip.border = list(col = 'grey95', lwd = 0.5),
         axis.line = list(col = "grey95", lwd = 0.5)
       ),
       col = "transparent",
       at = my.breaks,
       col.regions = myPalette,
       colorkey = list(
         space = "right",
         height = 1,
         width = 2,
         labels = list(cex = .75),
         at = my.breaks
       )
)

train_raw <- train.df

train_raw <- train_raw %>%
  mutate(State = substr(GEOID, 1, 2))

unscaled_summary_by_state <- train_raw %>%
  group_by(State) %>%
  dplyr::summarise(across(
    c(UNEMP, UNINSUR, DISABL, PARK, MHLTH, TRI, PM, WLKIND, MOBILE),
    list(
      Mean = ~mean(.x, na.rm = TRUE),
      SD   = ~sd(.x, na.rm = TRUE),
      Min  = ~min(.x, na.rm = TRUE),
      Max  = ~max(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(
    cols = -State,
    names_to = c("Variable", ".value"),
    names_sep = "_"
  ) %>%
  dplyr::select(State, Variable, Mean, SD, Min, Max)

print(unscaled_summary_by_state)

knitr::kable(unscaled_summary_by_state, format = "markdown", digits = 2)

state_lookup <- tibble::tibble(
  State = c("10", "13", "24", "34", "36", "37", "42", "45", "51", "54"),
  State_Name = c("Delaware", "Georgia", "Maryland", "New Jersey", "New York",
                 "North Carolina", "Pennsylvania", "South Carolina", "Virginia", "West Virginia")
)

radar_df <- unscaled_summary_by_state %>%
  dplyr::select(State, Variable, Mean) %>%
  left_join(state_lookup, by = "State") %>%
  dplyr::select(State_Name, Variable, Mean) %>%
  pivot_wider(names_from = Variable, values_from = Mean)

radar_mat <- radar_df %>%
  column_to_rownames("State_Name") %>%
  as.data.frame()

radar_scaled <- as.data.frame(rbind(
  apply(radar_mat, 2, max, na.rm = TRUE),
  apply(radar_mat, 2, min, na.rm = TRUE),
  radar_mat
))
rownames(radar_scaled)[1:2] <- c("Max", "Min")

state_labels <- rownames(radar_scaled)[-c(1, 2)]
n_states <- length(state_labels)
ncol <- 4
nrow <- ceiling(n_states / ncol)

par(mfrow = c(nrow, ncol), mar = c(1, 1, 2, 1))

for (state in state_labels) {
  radarchart(radar_scaled[c("Max", "Min", state), ],
             axistype = 1,
             pcol = "#1f77b4",
             pfcol = alpha("#1f77b4", 0.3),
             plwd = 2,
             title = state,
             cglcol = "grey", cglty = 1,
             axislabcol = "grey30",
             caxislabels = seq(-2, 2, 1),
             vlcex = 0.6)
}
