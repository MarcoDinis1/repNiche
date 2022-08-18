library(ggplot2)

input <- read.csv2(".\\input_sub7.txt")
sub.input <- subset(input, input$mode != "Larv")

# Change color by groups
dp <- ggplot(sub.input, aes(x=mode, y=vpd_1, fill=mode)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.05, fill="white")+
  labs(title="Vapor pressure deficit (VPD)",x="Mode", y = "VPD") + scale_fill_manual(values=c("#FFAA00","#0AF5F5", "#0A74F5", "#F50A0A")) #color codes derived from filling in function rgb() with rgb values from ArcGIS (for color homogeneity across figures). The order corresponds to :Alps, lycia, algira, ssal 
dp


###--------------###
##Paired histograms (for slope)
library(plyr)
mu <- ddply(input, "rep", summarise, grp.mean=mean(slope_1))
head(mu)


dp.slope <- ggplot(input, aes(x=slope_1, color=rep, fill=rep)) +
  geom_density(alpha=0.2)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=rep),
             linetype="dashed") +
  scale_color_manual(values=c("#E3E38D","#5331c7"))+
  scale_fill_manual(values=c("#E3E38D","#5331c7"))+
  labs(title="Slope",x="Slope (dd)", y = "Density")+
  theme_classic()
dp.slope


###paired histograms### 
#for all variables/comparisons among peuriparous groups, in order to show where perfect separation occurs. Edit manually for variable of interest
input.hist <- subset(input, input$mode == "Viv_Lycia" | input$mode == "Viv_Salg")

dp.hist <- ggplot(input.hist, aes(x=Kbin_1, color=mode, fill=mode)) +
  geom_density(alpha=0.2)+
  scale_color_manual(values=c("#0AF5F5","#0A74F5"))+
  scale_fill_manual(values=c("#0AF5F5","#0A74F5"))+
  labs(title="Karst presence/absence (KBIN)",x="KBIN", y = "Density")+
  theme_classic()
dp.hist


#Viv_Alp Viv_Lycia Viv_Salg Viv_Ssal
#bio01_1 Annual Mean Temperature (AT) AT (ºC)
#forestCo_1 Forest Cover (FORESTCO) FORESTCO (%)
#slope_1 Slope (SLOPE) SLOPE (dd)
#soil_1 Soil Moisture (SOIL) SOIL (mm)
#vpd_1 Vapor Pressure Deficit (VPD) VPD (kPa)
#Kbin_1 Karst presence/absence (KBIN) KBIN
