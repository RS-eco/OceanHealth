#' "A global assessment of ocean health using the trophic status index (TRIX)"
#' R-Code from RS-eco

#' <!--
#' Remove points and replace with geom_errorbar and geom_errorbarh, done for Fig.5
#' Update Halpern Analysis with 2013 data
#' -->

## ---- markdown_options ----

# Set global chunk setttings
knitr::opts_chunk$set(cache=TRUE, eval=TRUE, warning=FALSE, results = FALSE,
                      message=FALSE, comment=NA, echo=FALSE,
                      tidy=TRUE, fig.width=8, fig.height=4, 
                      fig.path='figures/', dev="png")

## ---- load_packages ----

# Install required packages not yet in library
packages <- c("dplyr", "ggplot2", "ggpmisc", "grid", "ggsignif", "mapdata", "maptools", 
              "pander", "parallel", "raster", "rgdal", "rgeos", "sf", "ncdf4")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load required packages
l <- sapply(packages, require, character.only = TRUE); rm(packages, l)

## ---- global_options ----

# Specify path of working directory and file directory
workdir <- "/home/matt/Documents/OceanHealth/"
filedir <- "/media/matt/Data/Documents/Wissenschaft/Data/"

# Set working directory 
setwd(workdir)

# Define spatial projection
crs.wgs84 <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Define colour theme for plotting
colourtheme <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                  "#FF7F00", "red", "#7F0000"))(255)
trix_colours <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", 
                                   "#FF7F00", "red", "#7F0000"))

# Set plotting theme
theme_set(theme_bw() + theme(
  axis.title.x = element_text(size=16),
  axis.title.y = element_text(size=16, angle=90),
  axis.text.x = element_text(size=14),
  axis.text.y = element_text(size=14),
  legend.text=element_text(size=14),
  legend.title=element_text(size=16),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "right"))

# Load world data for mapping
data(countries10, package="rnaturalearthhires")
countries10 <- sf::st_as_sf(countries10)

# Open shapefile of MEOW boundaries as Spatial Polygons
# File was downloaded from http://www.marineregions.org/
data(meow, package="geodat")

# Create realm shapefile 
sp_realm <- rgeos::gUnaryUnion(as(meow, "Spatial"), meow$REALM)
crs(sp_realm) <- crs.wgs84

# Run this
coastline <- sf::st_union(countries10)
coastline <- as(coastline, "Spatial")

#Crop realm by coastline
crs(coastline) <- crs.wgs84
marinerealms <- rgeos::gDifference(sp_realm, coastline, byid=TRUE)

# Create realm outline for nice plots
df_realm <- fortify(marinerealms, region="id")

# Define resolution
resolution <- "9km" # Alternatively: 4km

# Define time frame
years <- 2002:2015

## ---- functions ----

# Load trix function
source("R/trix_function.R")

## ---- chlorophyll ----

if(!file.exists(paste0(workdir, "data/chl_", resolution, "_realm.rds"))){
  # Path specification of directory for OceanColor files
  chl_files <- list.files(paste0(filedir, "OceanColor/Annual_", resolution, "_Chlorophyll_OCI/"), 
                          pattern=paste0("*.L3m_YR_CHL_chlor_a_", resolution, ".nc"), full.names=T)
  
  #Load nc files as raster, chl in mg/m3
  chlorophyll <- stack(chl_files[1:14], varname="chlor_a"); rm(chl_files)
  names(chlorophyll) <- c(2002:2015)
  projection(chlorophyll)  <- crs.wgs84
  
  # Mask chlorophyll by marine ecoregions
  chlorophyll <- mask(chlorophyll, marinerealms)
  chlorophyll <- as.data.frame(rasterToPoints(chlorophyll))
  saveRDS(chlorophyll, file=paste0(workdir, "data/chl_", resolution, "_realm.rds"), compress="xz")
} else{
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  chlorophyll <- rasterFromXYZ(chlorophyll)
}

## ---- nit_pho_d02 ----

# What about temporal variability in TRIX from WOD?
if(!file.exists(paste0(workdir, "data/wod_dO2_2015.rds")) | 
   !file.exists(paste0(workdir, "data/wod_nit_2015.rds")) | 
   !file.exists(paste0(workdir, "data/wod_pho_2015.rds"))){
  
  if(!file.exists(paste0(workdir, "data/wod_2002_2015.rds"))){
    # Read in exported dataset from Ocean Data View - Dowloaded from World Ocean Database
    files <- c("WOD_2015/data_from_WOD_2015_03.txt", "WOD_2015/data_from_WOD_2015_04.txt", 
               "WOD_2015/data_from_WOD_2015_05.txt")
    # WOD_2015_01 only goes from 1878 to 1977, WOD_2015_02 only goes from 1977 - 1994 and 
    # WOD_2015_03 goes from 1994 to 2011.
    wod <- lapply(files, FUN=function(x){
      wod <- read.delim(paste0(filedir, x),sep="\t",
                        header=FALSE, skip=32, comment.char="/", na.strings = c("", "NA"))
      colnames(wod) <- c("Cruise", "Station", "Type", "Date",
                         "Longitude",	"Latitude",	"Depth_m", "Temperature_deg_C", "Salinity_psu",
                         "Oxygen_mll",	"Phosphate_mmoll", "Nitrate_mmoll", "Nitrite_mmoll", 
                         "Chlorophyll_mgl")
      
      # Remove readings that are not the first cast from a cruise
      #wod <- subset(wod, !is.na(wod$Cruise))
      
      # Remove depths that are not close to the surface
      wod <- subset(wod, wod$Depth_m<20)
      
      # Clean year field
      wod$year <- as.numeric(substr(wod$Date,1,4))
      min(wod$year, na.rm=TRUE); max(wod$year, na.rm=TRUE)
      
      # Only select data from 2002 - 2015
      wod <- subset(wod, wod$year >= 2002 & wod$year <= 2015)
      return(wod)
    })
    
    # Merge list into one dataframe
    wod <- do.call("rbind", wod)
    
    # Set zeros to NAs, as they are likely erroneous
    wod$Oxygen_mll[wod$Oxygen_mll<0.00000001] <- NA
    wod$Phosphate_mmoll[wod$Phosphate_mmoll<0.00000001] <- NA
    wod$Nitrate_mmoll[wod$Nitrate_mmoll<0.00000001] <- NA
    wod$Nitrite_mmoll[wod$Nitrite_mmoll<0.00000001] <- NA
    wod$Chlorophyll_mgl[wod$Chlorophyll_mgl<0.00000001] <- NA
    
    # Save wod to file
    saveRDS(wod, "data/wod_2002_2015.rds")
  }
  
  # Read WOD data
  wod <- readRDS("data/wod_2002_2015.rds")
  
  # Extract Nitrate data
  wod_nit <- wod[,c("Longitude", "Latitude", "Date", "Nitrate_mmoll")]
  
  # Remove NAs
  wod_nit <- na.omit(wod_nit)
  
  #Convert nitrate from micro mol/l in mg/m3 (Nitrate (NO3): 61.9887 g/mol)
  wod_nit$nitrate <- (0.0619887*wod_nit$Nitrate_mmoll)
  
  # Save to file
  saveRDS(wod_nit, "data/wod_nit_2015.rds")
  
  # Extract Phosphate data
  wod_pho <- wod[,c("Longitude", "Latitude", "Date", "Phosphate_mmoll")]
  
  # Remove NAs
  wod_pho <- na.omit(wod_pho)
  
  #Convert phosphate from micro mol/l in mg/m3 (Phosphate (PO4): 94.9716 g/mol)
  wod_pho$phosphate <- (0.0949716*wod_pho$Phosphate_mmoll)
  
  # Save to file
  saveRDS(wod_pho, "data/wod_pho_2015.rds")
  
  # Remove rows that do not share all measurements
  wod_O2 <- subset(wod, !is.na(wod$Oxygen_mll))
  wod_O2 <- subset(wod_O2, !is.na(wod_O2$Temperature_deg_C))
  wod_O2 <- subset(wod_O2, !is.na(wod_O2$Salinity_psu))
  wod_O2 <- wod_O2[,c("Longitude", "Latitude", "Date", "Oxygen_mll", "Temperature_deg_C", "Salinity_psu")]
  
  #Calculate oxygen solubility
  abs_sst <- (wod_O2$Temperature_deg_C+273.15) #Absolute mean sst in Kelvin
  sst_f <- (abs_sst/100) #sst/100 = sst_f
  f_sst <- (100/abs_sst) #100/sst = f_sst
  a <- (-173.4292) + (249.6339*f_sst) + (143.3483*log(sst_f)) + ((-21.8492)*sst_f) #A1 + A2*(100/sst) + A3*Ln(sst/100) + A4*(sst/100)
  b <- (-0.033096) + (0.014259*sst_f) + ((-0.0017)*sqrt(sst_f)) #B1 + B2*(sst/100) + B3*square(sst/100)
  O2_sol <- exp(a + (wod_O2$Salinity_psu*b)) #Ln of a + Salinity*b
  
  #Remove unused rasters
  rm(abs_sst, f_sst, sst_f, a, b)
  
  #Calculate absolute deviation of % oxygen saturation
  perO2_sat <- (wod_O2$Oxygen_mll/O2_sol)*100 #Percentage oxygen saturation
  f_O2sat <- (100 - perO2_sat) #Calculate 100 - % Oxygen saturation
  wod_O2$dO2sat <- abs(f_O2sat) #Calculate absolute deviation of oxygen saturation
  rm(perO2_sat, f_O2sat)
  
  # Save to file
  saveRDS(wod_O2, "data/wod_dO2_2015.rds")
}

if(!file.exists(paste0(workdir, "data/nit_pho_dO2_", resolution, "_realm.rds"))){
  # Convert spatial points data frame into raster stack
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  chlorophyll <- raster::rasterFromXYZ(chlorophyll)
  r <- raster(resolution=c(1, 1), extent(chlorophyll), crs=crs.wgs84); rm(chlorophyll)
  
  # Inverse distance interpolation with inverse distance power set to .5:
  if(!file.exists(paste0(workdir, "data/nitrate_", resolution, "_realm.rds"))){
    wod_nit <- readRDS("data/wod_nit_2015.rds")
    colnames(wod_nit) <- c("x", "y", "date", "nitrate_mmoll", "nitrate")
    idw <- gstat::gstat(id = "nitrate", formula = nitrate~1, locations = ~x+y, data=wod_nit,
                        nmax=7, set=list(idp = .5))
    nitrate <- interpolate(object=r, model=idw)
    nitrate <- mask(nitrate, marinerealms); rm(idw)
    nitrate <- as.data.frame(rasterToPoints(nitrate))
    saveRDS(nitrate, file=paste0(workdir, "data/nitrate_", resolution, "_realm.rds"), compress="xz")
  } else{
    nitrate <- readRDS(paste0(workdir, "data/nitrate_", resolution, "_realm.rds"))
    nitrate <- raster::rasterFromXYZ(nitrate)
  }
  
  if(!file.exists(paste0(workdir, "data/phosphate_", resolution, "_realm.rds"))){
    wod_pho <- readRDS("data/wod_pho_2015.rds")
    colnames(wod_pho) <- c("x", "y", "date", "phosphate_mmoll", "phosphate")
    idw <- gstat::gstat(id = "phosphate", formula = phosphate~1, locations = ~x+y, data=wod_pho, 
                        nmax=7, set=list(idp = .5))
    phosphate <- interpolate(object=r, model=idw)
    phosphate <- mask(phosphate, marinerealms); rm(idw)
    phosphate <- as.data.frame(rasterToPoints(phosphate))
    saveRDS(phosphate, file=paste0(workdir, "data/phosphate_", resolution, "_realm.rds"), compress="xz")
  } else{
    phosphate <- readRDS(paste0(workdir, "data/phosphate_", resolution, "_realm.rds"))
    phosphate <- raster::rasterFromXYZ(phosphate)
  }
  
  if(!file.exists(paste0(workdir, "data/dO2sat_", resolution, "_realm.rds"))){
    wod_dO2 <- readRDS("data/wod_dO2_2015.rds")
    colnames(wod_dO2) <- c("x", "y", "date", "oxygen_mll", "temp_deg", "salinity_psu", "dO2sat")
    idw <- gstat::gstat(id = "dO2sat", formula = dO2sat~1, locations = ~x+y, data=wod_dO2, 
                        nmax=7, set=list(idp = .5))
    abs_dO2sat <- interpolate(object=r, model=idw)
    abs_dO2sat <- mask(abs_dO2sat, marinerealms); rm(idw)
    abs_dO2sat <- as.data.frame(rasterToPoints(abs_dO2sat))
    saveRDS(abs_dO2sat, file=paste0(workdir, "data/abs_dO2sat_", resolution, "_realm.rds"), compress="xz")
  } else{
    abs_dO2sat <- readRDS(paste0(workdir, "data/dO2sat_", resolution, "_realm.rds"))
    abs_dO2sat <- raster::rasterFromXYZ(abs_dO2sat)
  }
  
  #Stack and save to file
  nit_pho_dO2 <- stack(nitrate, phosphate, abs_dO2sat)
  #Resample WOA data to higher resolution of chlorophyll data
  chlorophyll <- stack(paste0(workdir, "data/chl_", resolution, "_realm.tif"))
  beginCluster()
  nit_pho_dO2 <- resample(nit_pho_dO2, chlorophyll[[1]])
  endCluster(); rm(chlorophyll)
  nit_pho_dO2 <- as.data.frame(rasterToPoints(nit_pho_dO2))
  saveRDS(nit_pho_dO2, file=paste0(workdir, "data/nit_pho_dO2_", resolution, "_realm.rds"), compress="xz")
} else{
  nit_pho_dO2 <- readRDS(paste0(workdir, "data/nit_pho_dO2_", resolution, "_realm.rds"))
  nit_pho_dO2 <- raster::rasterFromXYZ(nit_pho_dO2)
}

## ---- trix_realm ----

#Calculate trix for realm and province
if(!file.exists(paste0(workdir, "data/trix_", resolution, "_realm.rds"))){
  # Read chlorophyll data
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  chlorophyll <- raster::rasterFromXYZ(chlorophyll)
  
  # Read nutrient data
  nit_pho_dO2 <- readRDS(paste0(workdir, "data/nit_pho_dO2_", resolution, "_realm.rds"))
  nit_pho_dO2 <- raster::rasterFromXYZ(nit_pho_dO2)
  
  #Calculate mean chlorophyll over time (2002 - 2015)
  beginCluster(n=detectCores()-1)
  mean_chl <- clusterR(chlorophyll, fun = calc, args = list(fun = mean))
  endCluster(); rm(chlorophyll)
  
  #Stack environmental data together
  env_data <- stack(mean_chl, nit_pho_dO2); rm(mean_chl, nit_pho_dO2)
  names(env_data) <- c("chl", "nit", "pho", "dO2sat")
  
  trix_realm <- trix(x=env_data, y=marinerealms)
  trix_realm <- as.data.frame(rasterToPoints(trix_realm))
  saveRDS(trix_realm, file=paste0(workdir, "data/trix_", resolution, "_realm.rds"),
          compress="xz")
} else{
  trix_realm <- readRDS(paste0(workdir, "data/trix_", resolution, "_realm.rds"))
  trix_realm <- rasterFromXYZ(trix_realm)
}

# Convert raster to dataframe for plotting with ggplot2
trix_realm_df <- readRDS(paste0(workdir, "data/trix_", resolution, "_realm.rds"))
colnames(trix_realm_df) <- c("x","y","trix")

#ggplot(data = trix_realm_df, aes(trix)) + geom_histogram(bins=50, colour="white", fill="grey40") + 
#  labs(x = "TRIX", y = "Frequency")

## ---- trophic_state ----

# Load trix from file
trix_realm_df <- readRDS(paste0(workdir, "data/trix_", resolution, "_realm.rds"))
colnames(trix_realm_df) <- c("x","y","trix")

# Classify TRIX into trophic state categories
trix_realm_df$health <- NA 
trix_realm_df$health[trix_realm_df$trix >= 8] <- NA 
trix_realm_df$health[trix_realm_df$trix < 2] <- NA 
trix_realm_df$health[trix_realm_df$trix >= 2 & trix_realm_df$trix < 4] <- 1 
trix_realm_df$health[trix_realm_df$trix >= 4 & trix_realm_df$trix < 5] <- 2 
trix_realm_df$health[trix_realm_df$trix >= 5 & trix_realm_df$trix < 6] <- 3 
trix_realm_df$health[trix_realm_df$trix >= 6 & trix_realm_df$trix < 8] <- 4

# Trophic state by marine realm
if(!file.exists(paste0(workdir, "data/ts_realm.rds"))){
  trix_realm <- rasterFromXYZ(trix_realm_df)
  ts_realm <- extract(trix_realm, sp_realm, fun=mean, df=TRUE, na.rm=TRUE)
  saveRDS(ts_realm, "data/ts_realm.rds", compress="xz")
} else{
  ts_realm <- readRDS("data/ts_realm.rds")
}

# Classify TRIX into trophic state categories
ts_realm$health <- NA 
ts_realm$health[ts_realm$trix_REALM >= 8] <- NA 
ts_realm$health[ts_realm$trix_REALM < 2] <- NA 
ts_realm$health[ts_realm$trix_REALM >= 2 & ts_realm$trix_REALM < 4] <- 1 
ts_realm$health[ts_realm$trix_REALM >= 4 & ts_realm$trix_REALM < 5] <- 2 
ts_realm$health[ts_realm$trix_REALM >= 5 & ts_realm$trix_REALM < 6] <- 3 
ts_realm$health[ts_realm$trix_REALM >= 6 & ts_realm$trix_REALM < 8] <- 4

ts_realm$ID <- factor(ts_realm$ID, 
                      labels=c("Arctic", "Central Indo-Pacific", "Eastern Indo-Pacific", 
                               "Southern Ocean", "Temperate Australasia", 
                               "Temperate Northern Atlantic", "Temperate Northern Pacific", 
                               "Temperate South America", "Temperate Southern Africa", 
                               "Tropical Atlantic", "Tropical Eastern Pacific", "Western Indo-Pacific"))

## ---- fig1 ----

if(!file.exists(paste0(workdir, "figures/Figure_01.png"))){
  # Two panel plot including TRIX and categorised TRIX
  p1 <- ggplot() + 
    geom_raster(data=trix_realm_df, aes(x,y, fill=trix)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), fill=NA, colour="black") + 
    coord_sf(ndiscr=0, expand=F) + labs(x = "Longitude", y = "Latitude") +ggtitle("a)") + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("TRIX", limits=c(0,10), breaks=seq(0,10,2), 
                         colours=colourtheme, na.value="transparent") +
    theme(legend.key = element_blank(), plot.title=element_text(size=18, hjust=-0.10))
  
  # Create plot of Trophic state
  p2 <- ggplot() + 
    geom_raster(data=trix_realm_df, aes(x,y,fill=factor(health))) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill=NA, colour="black") + coord_sf(ndiscr=0, expand=F) +
    labs(x = "Longitude", y = "Latitude") + ggtitle("b)") +
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_manual("Status", breaks=c(4,3,2,1), labels = c("Bad", "Poor","Moderate","Good and above"), 
                      values=c("blue", "cyan", "yellow", "red"), na.value="transparent") +
    theme(legend.key = element_blank(), plot.title=element_text(size=18, hjust=-0.10))
  
  # Issues with Legends require Grob workaround for Multiplot following: 
  # http://stackoverflow.com/questions/16255579/how-can-i-make-consistent-width-plots-in-ggplot-with-legends
  gA <- ggplotGrob(p1)
  gB <- ggplotGrob(p2)
  # The parts that differs in width
  leg1 <- convertX(sum(with(gA$grobs[[15]], grobs[[1]]$widths)), "mm")
  leg2 <- convertX(sum(with(gB$grobs[[15]], grobs[[1]]$widths)), "mm")
  # Add an empty column of "abs(diff(widths)) mm" width on the right of 
  # legend box for gA (the smaller legend box)
  gA$grobs[[15]] <- gtable::gtable_add_cols(gA$grobs[[15]], unit(abs(diff(c(leg1, leg2))), "mm"))
  # Combine the plots
  g = gridExtra::gtable_rbind(gA, gB, size = "max")
  grid.newpage()
  png(file=paste0(getwd(),"/figures/Figure_01.png"), 
      width=10, height=8, units="in", res=300)
  grid.draw(g)
  dev.off()
  rm(g,p1,p2,gA,gB,leg1,leg2)
}

if(!file.exists(paste0(workdir, "figures/Figure_01b.png"))){
  # Pattern of TRIX with Latitude
  p <- ggplot(data = trix_realm_df, aes(x = abs(y), y = trix)) + 
    labs(x = "Latitude", y = "Trophic index") +
    scale_x_continuous(expand=c(0.01,0.01), limits=c(0,90), breaks=seq(0,90,15)) + 
    scale_y_continuous(expand=c(0.01,0.01), limits=c(0,10), breaks=seq(0,10,2)) +
    geom_point(shape=1, alpha = 1/25) + geom_smooth() + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 colour="black", label.x.npc = "right", label.y.npc = "bottom", parse=TRUE)
  ggsave(paste0("figures/Figure_01b.png"), p, width=8, height=6, units="in", dpi=300)
}; invisible(gc())

## ---- env_realm ----
if(!file.exists(paste0(workdir, "data/ma_env_realm.rds"))){
  # Read chlorophyll data
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  chlorophyll <- raster::rasterFromXYZ(chlorophyll)
  
  # Read nutrient data
  nit_pho_dO2 <- readRDS(paste0(workdir, "data/nit_pho_dO2_", resolution, "_realm.rds"))
  nit_pho_dO2 <- raster::rasterFromXYZ(nit_pho_dO2)
  
  #Calculate mean chlorophyll over time (2002 - 2015)
  beginCluster(n=detectCores()-1)
  mean_chl <- clusterR(chlorophyll, fun = calc, args = list(fun = mean))
  endCluster(); rm(chlorophyll)
  
  #Stack environmental data together
  ma_env_realm <- stack(mean_chl, nit_pho_dO2); rm(mean_chl, nit_pho_dO2)
  names(ma_env_realm) <- c("chl", "nit", "pho", "dO2sat")
  ma_env_realm <- as.data.frame(rasterToPoints(ma_env_realm))
  saveRDS(ma_env_realm, file="data/ma_env_realm.rds", compress="xz")
} else{
  ma_env_realm <- readRDS("data/ma_env_realm.rds")
  ma_env_realm <- rasterFromXYZ(ma_env_realm)
}

# Read realm data
r_realm <- rasterize(sp_realm, ma_env_realm)
ma_env_realm <- stack(ma_env_realm, r_realm); rm(r_realm)
ma_env_realm_df <- data.frame(rasterToPoints(ma_env_realm))
colnames(ma_env_realm_df) <- c("x", "y", "chl", "nit", "pho", "dO2sat", "realm")
ma_env_realm_df$realm <- factor(ma_env_realm_df$realm, 
                                labels=c("Arctic", "Central Indo-Pacific", "Eastern Indo-Pacific", 
                                         "Southern Ocean", "Temperate Australasia", 
                                         "Temperate Northern Atlantic", "Temperate Northern Pacific", 
                                         "Temperate South America", "Temperate Southern Africa", 
                                         "Tropical Atlantic", "Tropical Eastern Pacific", "Western Indo-Pacific"))

env_realm_sum <- ma_env_realm_df %>% group_by(realm) %>% 
  summarise(mn_chl=mean(chl, na.rm=TRUE), mn_nit=mean(nit, na.rm=TRUE),
            mn_pho=mean(pho, na.rm=TRUE), mn_dO2=mean(dO2sat, na.rm=TRUE),
            sd_chl=sd(chl, na.rm=TRUE), sd_nit=sd(nit, na.rm=TRUE),
            sd_pho=sd(pho, na.rm=TRUE), sd_dO2=sd(dO2sat, na.rm=TRUE))

## ---- figS1 ----

if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_01.png"))){
  ma_env_realm <- readRDS("data/ma_env_realm.rds")
  ma_env_realm <- rasterFromXYZ(ma_env_realm)
  names(ma_env_realm) <- c("chl", "nit", "pho", "dO2sat")
  
  # Convert raster to dataframe for plotting with ggplot
  # you can only do this with the RStoolbox package loaded
  env_data_df <- data.frame(rasterToPoints(ma_env_realm))
  colnames(env_data_df) <- c("x","y","mn_chl", "nit", "pho", "dO2")
  
  #Create plot of environmental conditions, MEOW and country outline
  p1 <- ggplot() + 
    geom_raster(data=env_data_df, aes(x,y,fill=mn_chl)) + 
    geom_sf(data=countries10, fill = "gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + coord_sf() +
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("Chlorophyll", colours=colourtheme, 
                         trans = "log10", na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") + ggtitle("a)") +
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.14))
  p2 <- ggplot() + 
    geom_raster(data=env_data_df, aes(x,y,fill=nit)) + 
    geom_sf(data=countries10, fill = "gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + coord_sf() +
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("Phosphate", colours=colourtheme, 
                         trans = "log10", na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") + ggtitle("b)") +
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.14))
  p3 <- ggplot() + 
    geom_raster(data=env_data_df, aes(x,y,fill=pho)) + 
    geom_sf(data=countries10, fill = "gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + coord_sf() +
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("Nitrate", colours=colourtheme, 
                         trans = "log10", na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") + ggtitle("c)") +
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.14))
  p4 <- ggplot() + 
    geom_raster(data=env_data_df, aes(x,y,fill=dO2)) + 
    geom_sf(data=countries10, fill = "gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + coord_sf() +
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("aD%O", colours=colourtheme, 
                         trans = "log10", na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") + ggtitle("d)") +
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.14))
  
  #Convert to plottable grobs for GridExtra
  gA <- ggplotGrob(p1)
  gB <- ggplotGrob(p2)
  gC <- ggplotGrob(p3)
  gD <- ggplotGrob(p4)
  
  # Issues with Legends require Grob workaround for Multiplot following: http://stackoverflow.com/questions/17462504/align-edges-of-ggplot-chloropleth-legend-title-varies?rq=1
  grid.newpage()
  png(file=paste(getwd(),"/figures/SupplementaryFigure_01.png",sep=""), 
      width=6, height=12, units="in", res=300)
  g <- rbind(gA,gB,gC,gD, size = "first")
  for (i in which(g$layout$name == "guide-box")) {
    g$grobs[[i]] <- g$grobs[[i]]$grobs[[1]]
  }
  grid.draw(g)
  dev.off()
  rm(g,p1,p2,p3,p4,gA,gB,gC,gD)
}; invisible(gc())

## ---- figS2 ----

if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_02.png"))){
  ma_env_realm <- readRDS("data/ma_env_realm.rds")
  colnames(ma_env_realm) <- c("x", "y", "mn_chl", "nit", "pho", "dO2")
  
  trix_realm <- readRDS(paste0(workdir, "data/trix_", resolution, "_realm.rds"))
  colnames(trix_realm) <- c("x", "y", "trix")
  
  env_trix_df <- left_join(ma_env_realm, trix_realm) %>% dplyr::select(-c(x,y)) %>% 
    tidyr::gather("var", "value", -c(trix)) %>% tidyr::drop_na() %>%
    mutate(var = factor(var, labels=c("Log chl-a", "Log nitrate", "Log phosphate", "Log aD%O")))
  
  #Plot histogram of scaled environmental variables and trix
  #ggplot(data = env_trix_df, aes(value)) + 
  #  geom_histogram(bins=50, col="white", fill="gray40") + 
  #  labs(x = "Chlorophyll a", y = "Frequency") + 
  #  facet_wrap(~ var, ncol=2, scales="free_x", strip.position="bottom") + 
  #  theme(strip.background= element_blank(), strip.placement="outside", 
  #        strip.text.x = element_text(size=16))
  
  p <- ggplot(data = env_trix_df, aes(x = value, y = trix)) + 
    geom_point(shape=1, alpha = 1/25) + geom_smooth() + 
    scale_x_log10() + scale_y_continuous(expand=c(0.01,0.01), limits=c(0,10), breaks=seq(0,10,2)) + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 parse = TRUE, label.x.npc = "left", label.y.npc = "top", colour="black") + 
    facet_wrap(.~ var, ncol=2, scales="free_x", strip.position="bottom") + 
    labs(x="", y="TRIX") + 
    theme(strip.background= element_blank(), strip.placement="outside", 
          strip.text.x = element_text(size=16))
  ggsave("figures/SupplementaryFigure_02.png", p, width=9, height=8, units="in", dpi=300)
}; invisible(gc())

## ---- ref_cond ----   

if(!file.exists(paste0(workdir, "data/ref_cond_table.rds"))){
  # Load trix from file
  trix_realm <- readRDS(paste0(workdir, "data/trix_", resolution, "_realm.rds"))
  trix_realm <- rasterFromXYZ(trix_realm)
  
  #Load reference conditions
  ref_cond <- readRDS("data/ref_cond_REALM.rds")[,-c(6)]
  colnames(ref_cond) <- c("zone", "chl_min", "nit_min", "pho_min", "dO2sat_min",
                          "chl_max", "nit_max", "pho_max", "dO2sat_max")
  
  # Calculate mean trix per realm
  ref_cond$mean_trix <- extract(trix_realm, sp_realm, fun=mean, na.rm=TRUE)
  
  # Define names of each marine realm
  ref_cond$zone <- levels(meow$REALM)
  
  # Calculate SD of each marine realm
  ref_cond$sd_trix <- extract(trix_realm, sp_realm, fun=sd, na.rm=TRUE)
  
  # Save Ref cond for Table
  saveRDS(ref_cond, "data/ref_cond_table.rds", compress="xz")
  file.remove("data/ref_cond_REALM.rds")
} else{
  ref_cond <- readRDS("data/ref_cond_table.rds")
}; invisible(gc())

## ---- table1 ----

## Table 1. Reference conditions and mean trophic index (TRIX) for each marine realm
ref_cond_table <- readRDS("data/ref_cond_table.rds")

# Define names of dataframe
colnames(ref_cond_table) <- c("Marine realm", "Min Chl-a", "Min Nit", 
                              "Min Pho", "Min aD%O", "Max Chl-a", "Max Nit", 
                              "Max Pho", "Max aD%O", "Mean TRIX")

# Produce table of reference conditions
pander::pandoc.table(ref_cond_table, round=2, split.tables=Inf)

## ---- lin_reg_env ----

# Read data
ma_env_realm <- readRDS(paste0(workdir, "data/ma_env_realm.rds"))
trix_realm <- readRDS(paste0(workdir, "data/trix_", resolution, "_realm.rds"))
colnames(ma_env_realm) <- c("x", "y", "mn_chl", "nit", "pho", "dO2")
colnames(trix_realm) <- c("x", "y", "trix")

# Open Env-TRIX data
env_trix_df <- ma_env_realm %>% group_by(x,y) %>% mutate_at(vars(-group_cols()), log10) %>%
  left_join(trix_realm) %>% tidyr::drop_na(trix) %>% as.data.frame()
rm(ma_env_realm, trix_realm); invisible(gc())

# Create data frame for results of linear model
lin_reg_env <- data.frame(matrix(nrow=4, ncol=4))

# Run linear regression between each explanatory variable 
# and the different environmental impacts 
for (i in 1:4){
  # Run individual model
  m <- lm(env_trix_df$trix ~ env_trix_df[,2+i], na.action="na.omit")
  
  # Extract R2 value and save to dataframe
  lin_reg_env[i,1] <- c("Chl-a", "Nitrate", "Phosphate", "aD%O2")[i]
  lin_reg_env[i,2] <- summary(m)$adj.r.squared
  lin_reg_env[i,3] <- summary(m)$coefficients[2]
  lin_reg_env[i,4] <- summary(m)$coefficients[8]
}; rm(env_trix_df); invisible(gc())

#Define names of dataframe
colnames(lin_reg_env) <- c("Variable", "AdjR2", "Intercept", "p-value")

#Produce table
pander::pandoc.table(lin_reg_env, round=4, split.tables=Inf)

## ---- figS3_4 ----

if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_03.png"))){
  # Load annual chlorophyll files
  an_chl <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  
  # Convert chlorophyll data into long df format for plotting
  years <- 2002:2015
  colnames(an_chl) <- c("x", "y", years)
  an_chl1 <- tidyr::gather(an_chl[,c(1:10)], "year", "chl", -c(x,y))
  an_chl2 <- tidyr::gather(an_chl[,c(1,2,11:ncol(an_chl))], "year", "chl",-c(x,y)); rm(an_chl)
  
  # Plot annual chlorophyll maps
  p1 <- ggplot() + 
    geom_raster(data=an_chl1, aes(x,y,fill=chl)) + 
    scale_fill_gradientn(name="Chl-a",
                         colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines")); rm(an_chl1)
  ggsave(paste0("figures/SupplementaryFigure_03.png"), width=10, height=12, units="in", dpi=300)
  p2 <- ggplot() + 
    geom_raster(data=an_chl2, aes(x,y,fill=chl)) + 
    scale_fill_gradientn(name="Chl-a", colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines")); rm(an_chl2)
  ggsave(paste0("figures/SupplementaryFigure_04.png"), width=10, height=9, dpi=300)
}

## ---- an_mn_chl ----

# Calculate mean and sd of chlorophyll for each year
if(!file.exists(paste0(workdir, "data/summary_an_mn_chl.rds"))){
  # Load annual chlorophyll files
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  chlorophyll <- rasterFromXYZ(chlorophyll)
  
  an_mn_chl <- data.frame(mean=cellStats(chlorophyll, stat=mean, na.rm=TRUE), 
                          sd=cellStats(chlorophyll, stat=sd, na.rm=TRUE), 
                          year=2002:2015)
  
  saveRDS(an_mn_chl, file="data/summary_an_mn_chl.rds", compress="xz"); rm(an_mn_chl) # Write data to file
}

# Calculate mean and sd of chlorophyll for each realm and each year
if(!file.exists(paste0(workdir, "data/summary_an_chl_realm.rds"))){
  # Load annual chlorophyll files
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  chlorophyll <-  rasterFromXYZ(chlorophyll)
  
  mn_chl <- extract(chlorophyll, sp_realm, fun=mean, df=TRUE, na.rm=TRUE)
  colnames(mn_chl) <- c("ID", 2002:2015)
  mn_chl <- tidyr::gather(mn_chl, "year", "mn_chl", -c(ID))
  sd_chl <- extract(chlorophyll, sp_realm, fun=sd, df=TRUE, na.rm=TRUE)
  colnames(sd_chl) <- c("ID", 2002:2015)
  sd_chl <- tidyr::gather(sd_chl, "year", "sd_chl", -c(ID))
  sum_chl_realm <- merge(mn_chl, sd_chl, by=c("ID", "year")); rm(mn_chl, sd_chl)
  sum_chl_realm <- plyr::join(sum_chl_realm, data.frame(ID=c(1:length(names(sp_realm))), 
                                                        Realm = names(sp_realm)), by="ID")
  colnames(sum_chl_realm) <- c("ID", "year", "mean", "sd", "realm")
  saveRDS(sum_chl_realm, "data/summary_an_chl_realm.rds", compress="xz"); rm(sum_chl_realm)
}

#if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_X.png"))){
#  # Read data
#  sum_chl_realm <- readRDS("data/summary_an_chl_realm.rds")
#  
#  # Time series graph with line of mean Chlorophyll for each realm
#  ggplot(sum_chl_realm, aes(year, mean)) + 
#    labs(list(x = "Year", y = "Chlorophyll a")) + 
#    geom_point(shape=1) + 
#    scale_x_continuous(breaks=seq(2002, 2015, 2)) +
#    scale_y_continuous(breaks=seq(-2.5, 5.5, 0.5)) + 
#    geom_smooth(method="lm", se=TRUE, formula = y ~ x, colour="black") + 
#    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                 parse = TRUE, colour="black") + 
#    facet_wrap(~ realm, ncol=3) + 
#    theme(legend.position = "none", strip.text.x = element_text(size=12, face="bold"), 
#          strip.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
#  ggsave(paste0("figures/SupplementaryFigure_X.png"), width=9, height=12, units="in", dpi=300)
#}

## ---- an_mn_chl_lin_reg ----

# Load annual mean chlorophyll data
an_mn_chl <- readRDS("data/summary_an_mn_chl.rds")

# Run linear regression
m_chl <- lm(mean ~ year, data=an_mn_chl)

# Load annual mean chlorophyll data per marine realm
sum_chl_realm <- readRDS("data/summary_an_chl_realm.rds")

# Create data frame for results of linear model
lin_reg_chl <- data.frame(matrix(nrow=12, ncol=4))

# Run linear regression between each explanatory variable and the different environmental impacts
for (i in 1:12){
  chl_sub <- sum_chl_realm %>% filter(realm == unique(sum_chl_realm$realm)[i])
  # Run individual model
  m <- lm(chl_sub$mean ~ chl_sub$year)
  
  # Extract R2 value and save to dataframe
  lin_reg_chl[i,1] <- unique(sum_chl_realm$realm)[i]
  lin_reg_chl[i,2] <- summary(m)$adj.r.squared
  lin_reg_chl[i,3] <- summary(m)$coefficients[2]
  lin_reg_chl[i,4] <- summary(m)$coefficients[8]
}; rm(chl_sub, m)

#Define names of dataframe
colnames(lin_reg_chl) <- c("Realm", "Adj. R2", "Intercept", "p-value")

#Produce table
pander::pandoc.table(lin_reg_chl, round=4, split.tables=Inf)

## ---- an_wod2015 ----

# Create overall layer
if(!file.exists(paste0(workdir, "data/an_nit_", resolution, "_realm.rds")) | 
   !file.exists(paste0(workdir, "data/an_pho_", resolution, "_realm.rds")) | 
   !file.exists(paste0(workdir, "data/an_dO2_", resolution, "_realm.rds"))){
  
  # Convert spatial points data frame into raster stack
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  chlorophyll <- rasterFromXYZ(chlorophyll)[[1]]
  r <- raster(resolution=c(1, 1), extent(chlorophyll), crs=crs.wgs84); rm(chlorophyll)
  
  # Subset WOD over multiple time periods
  years <- 2002:2015
  if(!file.exists(paste0(workdir, "data/an_nit_", resolution, "_realm.tif"))){
    wod_nit <- readRDS("data/wod_nit_2015.rds")
    colnames(wod_nit) <- c("x", "y", "date", "nitrate_mmoll", "nitrate")
    wod_nit$year <- as.numeric(substr(wod_nit$date,1,4))
    an_nitrate<- lapply(years, FUN=function(x){
      wod <- wod_nit[wod_nit$year == x,]
      # Inverse distance interpolation with inverse distance power set to .5:
      idw <- gstat::gstat(id = "nitrate", formula = nitrate~1, locations = ~x+y, data=wod,
                          nmax=7, set=list(idp = .5))
      nitrate <- interpolate(object=r, model=idw)
      nitrate <- mask(nitrate, marinerealms); rm(idw)
      return(nitrate)
    })
    
    #Stack data
    an_nitrate <- stack(an_nitrate)
    
    #Resample WOA data to higher resolution of chlorophyll data and save to file
    chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
    chlorophyll <- rasterFromXYZ(chlorophyll)[[1]]
    beginCluster()
    an_nitrate <- resample(an_nitrate, chlorophyll)
    endCluster(); rm(chlorophyll)
    an_nitrate <- as.data.frame(rasterToPoints( an_nitrate))
    saveRDS(an_nitrate, file=paste0(workdir, "data/ an_nitrate_", resolution, "_realm.rds"),
            compress="xz")
  }
  
  if(!file.exists(paste0(workdir, "data/an_pho_", resolution, "_realm.rds"))){
    wod_pho <- readRDS("data/wod_pho_2015.rds")
    colnames(wod_pho) <- c("x", "y", "date", "phosphate_mmoll", "phosphate")
    wod_pho$year <- as.numeric(substr(wod_pho$date,1,4))
    an_phosphate <- lapply(years, FUN=function(x){
      wod <- wod_pho[wod_pho$year == x,]
      # Inverse distance interpolation with inverse distance power set to .5:
      idw <- gstat::gstat(id = "phosphate", formula = phosphate~1, locations = ~x+y, data=wod, 
                          nmax=7, set=list(idp = .5))
      phosphate <- mask(interpolate(object=r, model=idw), marinerealms); rm(idw)
      return(phosphate)
    })
    an_phosphate <- stack(an_phosphate)
    chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
    chlorophyll <- rasterFromXYZ(chlorophyll)[[1]]
    beginCluster()
    an_phosphate <- resample(an_phosphate, chlorophyll)
    endCluster(); rm(chlorophyll)
    an_phosphate <- as.data.frame(rasterToPoints(an_phosphate))
    saveRDS(an_phosphate, file=paste0(workdir, "data/an_phosphate_", resolution, "_realm.rds"),
            compress="xz")
  }
  
  if(!file.exists(paste0(workdir, "data/an_dO2_", resolution, "_realm.rds"))){
    wod_dO2 <- readRDS("data/wod_dO2_2015.rds")
    colnames(wod_dO2) <- c("x", "y", "date", "oxygen_mll", "temperature_deg", "salinity_psu", "dO2sat")
    wod_dO2$year <- as.numeric(substr(wod_dO2$date,1,4))
    an_abs_dO2sat <- lapply(years, FUN=function(x){
      wod <- wod_dO2[wod_dO2$year == x,]
      idw <- gstat::gstat(id = "dO2sat", formula = dO2sat~1, locations = ~x+y, data=wod, 
                          nmax=7, set=list(idp = .5))
      abs_dO2sat <- mask(interpolate(object=r, model=idw), marinerealms); rm(idw)
      return(abs_dO2sat)
    })
    an_abs_dO2sat <- stack(an_abs_dO2sat)
    chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
    chlorophyll <- rasterFromXYZ(chlorophyll)[[1]]
    beginCluster()
    an_abs_dO2sat <- resample(an_abs_dO2sat, chlorophyll)
    endCluster(); rm(chlorophyll)
    an_abs_dO2sat <- as.data.frame(rasterToPoints(an_abs_dO2sat))
    saveRDS(an_abs_dO2sat, file=paste0(workdir, "data/an_dO2_", resolution, "_realm.rds"),
            compress="xz")
  }
}

## ---- an_mn_wod ----

# Calculate mean and sd of WOD data for each year
if(!file.exists(paste0(workdir, "data/summary_an_mn_dO2.rds"))){
  # Load annual WOD files
  nitrate <- readRDS(paste0(workdir, "data/an_nit_", resolution, "_realm.rds"))
  nitrate <- rasterFromXYZ(nitrate)
  phosphate <- readRDS(paste0(workdir, "data/an_pho_", resolution, "_realm.rds"))
  phosphate <- rasterFromXYZ(phosphate)
  absdO2 <- readRDS(paste0(workdir, "data/an_dO2_", resolution, "_realm.rds"))
  absdO2 <- rasterFromXYZ(absdO2)
  
  # Calculate mean and sd per year
  an_mn_nit <- data.frame(mean=cellStats(nitrate, stat=mean, na.rm=TRUE), 
                          sd=cellStats(nitrate, stat=sd, na.rm=TRUE),
                          year=2002:2015)
  an_mn_pho <- data.frame(mean=cellStats(phosphate, stat=mean, na.rm=TRUE), 
                          sd=cellStats(phosphate, stat=sd, na.rm=TRUE),
                          year=2002:2015)
  an_mn_dO2 <- data.frame(mean=cellStats(absdO2, stat=mean, na.rm=TRUE), 
                          sd=cellStats(absdO2, stat=sd, na.rm=TRUE), 
                          year=2002:2015)
  saveRDS(an_mn_nit, file="data/summary_an_mn_nit.rds", compress="xz"); rm(an_mn_nit) # Write data to file
  saveRDS(an_mn_pho, file="data/summary_an_mn_pho.rds", compress="xz"); rm(an_mn_pho) # Write data to file
  saveRDS(an_mn_dO2, file="data/summary_an_mn_dO2.rds", compress="xz"); rm(an_mn_dO2) # Write data to file
}

# Calculate mean and sd of chlorophyll for each realm and each year
if(!file.exists(paste0(workdir, "data/summary_an_dO2_realm.rds"))){
  # Load annual WOD files
  nitrate <- readRDS(paste0(workdir, "data/an_nit_", resolution, "_realm.rds"))
  nitrate <- rasterFromXYZ(nitrate)
  phosphate <- readRDS(paste0(workdir, "data/an_pho_", resolution, "_realm.rds"))
  phosphate <- rasterFromXYZ(phosphate)
  absdO2 <- readRDS(paste0(workdir, "data/an_dO2_", resolution, "_realm.rds"))
  absdO2 <- rasterFromXYZ(absdO2)
  
  mn_nit <- extract(nitrate, sp_realm, fun=mean, df=TRUE, na.rm=TRUE)
  colnames(mn_nit) <- c("ID", 2002:2015)
  mn_nit <- tidyr::gather(mn_nit, "year", "mn_nit", -c(ID))
  sd_nit <- extract(nitrate, sp_realm, fun=sd, df=TRUE, na.rm=TRUE); rm(nitrate)
  colnames(sd_nit) <- c("ID", 2002:2015)
  sd_nit <- tidyr::gather(sd_nit, "year", "sd_nit", -c(ID))
  sum_nit_realm <- merge(mn_nit, sd_nit, by=c("ID", "year")); rm(mn_nit, sd_nit)
  sum_nit_realm <- plyr::join(sum_nit_realm, data.frame(ID=c(1:length(names(sp_realm))), 
                                                        Realm = names(sp_realm)), by="ID")
  colnames(sum_nit_realm) <- c("ID", "year", "mean", "sd", "realm")
  saveRDS(sum_nit_realm, "data/summary_an_nit_realm.rds", compress="xz"); rm(sum_nit_realm)
  
  mn_pho <- extract(phosphate, sp_realm, fun=mean, df=TRUE, na.rm=TRUE)
  colnames(mn_pho) <- c("ID", 2002:2015)
  mn_pho <- tidyr::gather(mn_pho, "year", "mn_pho", -c(ID))
  sd_pho <- extract(phosphate, sp_realm, fun=sd, df=TRUE, na.rm=TRUE); rm(phosphate)
  colnames(sd_pho) <- c("ID", 2002:2015)
  sd_pho <- tidyr::gather(sd_pho, "year", "sd_pho", -c(ID))
  sum_pho_realm <- merge(mn_pho, sd_pho, by=c("ID", "year")); rm(mn_pho, sd_pho)
  sum_pho_realm <- plyr::join(sum_pho_realm, data.frame(ID=c(1:length(names(sp_realm))), 
                                                        Realm = names(sp_realm)), by="ID")
  colnames(sum_pho_realm) <- c("ID", "year", "mean", "sd", "realm")
  saveRDS(sum_pho_realm, "data/summary_an_pho_realm.rds", compress="xz"); rm(sum_pho_realm)
  
  mn_dO2 <- extract(absdO2, sp_realm, fun=mean, df=TRUE, na.rm=TRUE)
  colnames(mn_dO2) <- c("ID", 2002:2015)
  mn_dO2 <- tidyr::gather(mn_dO2, "year", "mn_dO2", -c(ID))
  sd_dO2 <- extract(absdO2, sp_realm, fun=sd, df=TRUE, na.rm=TRUE)
  colnames(sd_dO2) <- c("ID", 2002:2015)
  sd_dO2 <- tidyr::gather(sd_dO2, "year", "sd_dO2", -c(ID))
  sum_dO2_realm <- merge(mn_dO2, sd_dO2, by=c("ID", "year")); rm(mn_dO2, sd_dO2)
  sum_dO2_realm <- plyr::join(sum_dO2_realm, data.frame(ID=c(1:length(names(sp_realm))), 
                                                        Realm = names(sp_realm)), by="ID")
  colnames(sum_dO2_realm) <- c("ID", "year", "mean", "sd", "realm")
  saveRDS(sum_dO2_realm, "data/summary_an_dO2_realm.rds", compress="xz"); rm(sum_dO2_realm)
}

## ---- an_mn_wod_lin_reg ----

# Load annual mean env data
an_mn_nit <- readRDS("data/summary_an_mn_nit.rds")
an_mn_pho <- readRDS("data/summary_an_mn_pho.rds")
an_mn_dO2 <- readRDS("data/summary_an_mn_dO2.rds")

# Run linear regression
m_nit <- lm(mean ~ year, data=an_mn_nit)
m_pho <- lm(mean ~ year, data=an_mn_pho)
m_dO2 <- lm(mean ~ year, data=an_mn_dO2)

# Load annual nitrate data per marine realm
sum_nit_realm <- readRDS("data/summary_an_nit_realm.rds")

# Create data frame for results of linear model
lin_reg_nit <- data.frame()

# Run linear regression between each explanatory variable and the different environmental impacts
for (i in 1:length(unique(sum_nit_realm$realm))){
  nit_sub <- subset(sum_nit_realm, sum_nit_realm$realm == unique(sum_nit_realm$realm)[i])
  # Run individual model
  m <- lm(nit_sub$mean ~ as.numeric(nit_sub$year), na.action="na.omit")
  
  # Extract R2 value and save to dataframe
  lin_reg_nit[i,1] <- unique(sum_nit_realm$realm)[i]
  lin_reg_nit[i,2] <- summary(m)$r.squared
  lin_reg_nit[i,3] <- summary(m)$coefficients[2]
  lin_reg_nit[i,4] <- summary(m)$coefficients[8]
}; rm(nit_sub, m)

#Define names of dataframe
colnames(lin_reg_nit) <- c("Realm", "Adj. R2", "Intercept", "p-value")

#Produce table
pander::pandoc.table(lin_reg_nit, round=4, split.tables=Inf)

## ---- figS5_6 ----

if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_05.png"))){
  an_nit <- readRDS(paste0(workdir, "data/an_nit_", resolution, "_realm.rds"))
  # Convert nitrate data into long df format for plotting
  colnames(an_nit) <- c("x", "y", 2002:2015)
  an_nit1 <- tidyr::gather(an_nit[,c(1:10)], "year", "nit", -c(x,y))
  an_nit2 <- tidyr::gather(an_nit[,c(1,2,11:ncol(an_nit))], "year", "nit", -c(x,y)); rm(an_nit)
  an_nit1 <- na.omit(an_nit1)
  an_nit2 <- na.omit(an_nit2)
  
  # Plot annual nitrate maps  
  p1 <- ggplot() + 
    geom_raster(data=an_nit1, aes(x,y,fill=nit)) + 
    scale_fill_gradientn(name="Nitrate", colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines")); rm(an_nit1)
  ggsave("figures/SupplementaryFigure_05.png", width=10, height=12, units="in", dpi=300); rm(p1)
  p2 <- ggplot() + 
    geom_raster(data=an_nit2, aes(x,y,fill=nit)) + 
    scale_fill_gradientn(name="Nitrate", colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines")); rm(an_nit2)
  ggsave("figures/SupplementaryFigure_06.png", width=10, height=9, units="in", dpi=300); rm(p2)
}

## ---- figS7_8 ----

if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_07.png"))){
  an_pho <- readRDS(paste0(workdir, "data/an_pho_", resolution, "_realm.rds"))
  # Convert chlorophyll data into long df format for plotting
  colnames(an_pho) <- c("x", "y", 2002:2015)
  an_pho1 <- tidyr::gather(an_pho[,c(1:10)], "year", "pho", -c(x,y))
  an_pho2 <- tidyr::gather(an_pho[,c(1,2,11:ncol(an_pho))], "year", "pho",-c(x,y)); rm(an_pho)
  
  # Plot annual chlorophyll maps
  p1 <- ggplot() + 
    geom_raster(data=an_pho1, aes(x,y,fill=pho)) + 
    scale_fill_gradientn(name="Phosphate", breaks=c(0, 0.1,1,10,40), labels=c("0", "0.1","1","10","40"),
                         colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines"))
  ggsave(paste0("figures/SupplementaryFigure_07.png"), width=10, height=12, units="in", dpi=300); rm(an_pho1)
  p1 <- ggplot() + 
    geom_raster(data=an_pho2, aes(x,y,fill=pho)) + 
    scale_fill_gradientn(name="Phosphate", breaks=c(0, 0.1,1,10,40), labels=c("0", "0.1","1","10","40"),
                         colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines"))
  ggsave(paste0("figures/SupplementaryFigure_08.png"), width=10, height=9, units="in", dpi=300); rm(an_pho2)
}

## ---- figS9_10 ----

if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_09.png"))){
  an_dO2 <- readRDS(paste0(workdir, "data/an_dO2_", resolution, "_realm.rds"))
  # Convert dO2 data into long df format for plotting
  colnames(an_dO2) <- c("x", "y", 2002:2015)
  an_dO21 <- tidyr::gather(an_dO2[,c(1:10)], "year", "dO2", -c(x,y))
  an_dO22 <- tidyr::gather(an_dO2[,c(1,2,11:ncol(an_dO2))], "year", "dO2",-c(x,y)); rm(an_dO2)
  
  # Plot annual dO2 maps
  p1 <- ggplot() + 
    geom_raster(data=an_dO21, aes(x,y,fill=dO2)) + 
    scale_fill_gradientn(name="aD%O2", breaks=c(0, 0.1,1,10,40), labels=c("0", "0.1","1","10","40"),
                         colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines")); rm(an_dO21)
  ggsave(paste0("figures/SupplementaryFigure_09.png"), width=10, height=12, units="in", dpi=300)
  p2 <- ggplot() + 
    geom_raster(data=an_dO22, aes(x,y,fill=dO2)) + 
    scale_fill_gradientn(name="aD%O2", breaks=c(0, 0.1,1,10,40), labels=c("0", "0.1","1","10","40"),
                         colours=colourtheme, na.value="transparent", trans="log10") + 
    scale_x_continuous(name=expression(paste("Longitude (",degree,")")), 
                       expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(name=expression(paste("Latitude (",degree,")")), 
                       expand=c(0.01,3), breaks=seq(-90,90,30)) + 
    geom_sf(data=countries10, fill="gray", color="black") + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines")); rm(an_dO22)
  ggsave(paste0("figures/SupplementaryFigure_10.png"), width=10, height=9, units="in", dpi=300)
}

## ---- an_trix_realm ----

if(!file.exists(paste0(workdir, "data/an_trix_realm.rds"))){
  # Define years
  years <- 2002:2015
  
  # Load data
  chlorophyll <- readRDS(paste0(workdir, "data/chl_", resolution, "_realm.rds"))
  nitrate <- readRDS(paste0(workdir, "data/an_nit_", resolution, "_realm.rds"))
  phosphate <- readRDS(paste0(workdir, "data/an_pho_", resolution, "_realm.rds"))
  absdO2 <- readRDS(paste0(workdir, "data/an_dO2_", resolution, "_realm.rds"))
  
  chlorophyll <- rasterFromXYZ(chlorophyll)
  nitrate <- rasterFromXYZ(nitrate)
  phosphate <- rasterFromXYZ(phosphate)
  absdO2 <- rasterFromXYZ(absdO2)
  chlorophyll <- extend(chlorophyll, nitrate)
  invisible(gc())
  
  # Calculate yearly TRIX
  for(i in 1:length(years)){
    # Load yearly env data
    chl1 <- chlorophyll[[i]]
    nit1 <- nitrate[[i]]
    pho1 <- phosphate[[i]]
    absdO21 <- absdO2[[i]]
    
    # Stack data
    env_data <- stack(chl1, nit1, pho1, absdO21)
    rm(chl1, nit1, pho1, absdO21)
    
    # Calculate TRIX
    an_trix_realm <- trix(env_data, sp_realm, output=FALSE); rm(env_data)
    an_trix_realm <- as.data.frame(rasterToPoints(an_trix_realm))
    
    # Write annual trix to file
    saveRDS(an_trix_realm, paste0("data/an_trix_", years[i], "_realm.rds"), 
            compress="xz"); rm(an_trix_realm)
    
    # Print status
    print(i)
  }
  files <- list.files(path=paste0(workdir, "data"), 
                      pattern="an_trix_20", full.names=TRUE)
  an_trix_realm <- lapply(files, readRDS)
  an_trix_realm <- Reduce(function(...) dplyr::full_join(..., by=c("x","y"), all.x=TRUE), 
                          an_trix_realm)
  colnames(an_trix_realm) <- c("x", "y", 2002:2015)
  saveRDS(an_trix_realm, "data/an_trix_realm.rds", compress="xz"); file.remove(files)
}

## ---- fig2     ----

if(!file.exists(paste0(workdir, "figures/Figure_02.png"))){
  if(!file.exists(paste0(workdir, "data/lin_reg_an_trix.rds"))){
    # Fortify trix
    an_trix_realm <- readRDS("data/an_trix_realm.rds")
    an_trix_df <- t(an_trix_realm); rm(an_trix_realm)
    date <- seq(from = as.Date(paste(2002, "06-15", sep="-")), 
                to = as.Date(paste(2015, "06-15", sep="-")), 
                by="year")
    lin_reg_an_trix <- data.frame(x = an_trix_df[1,], 
                                  y = an_trix_df[2,])
    for(j in 1:ncol(an_trix_df)){
      # Run individual model
      m <- lm(an_trix_df[c(-1,-2),j] ~ date)
      
      # Extract R2 value and save to dataframe
      lin_reg_an_trix$slope[j] <- m$coefficients[2]
      lin_reg_an_trix$adjrsq[j] <- summary(m)$adj.r.squared
      lin_reg_an_trix$pvalue[j] <- summary(m)$coefficients[4]
    }; rm(an_trix_df)
    
    # Group into decreasing, stable, increase
    lin_reg_an_trix$var <- NA
    lin_reg_an_trix$var[lin_reg_an_trix$pvalue < 0.05 & lin_reg_an_trix$slope > 0] <- 1
    lin_reg_an_trix$var[lin_reg_an_trix$pvalue < 0.05 & lin_reg_an_trix$slope < 0] <- -1
    lin_reg_an_trix$var[lin_reg_an_trix$pvalue >= 0.05] <- 0
    saveRDS(lin_reg_an_trix, "data/lin_reg_an_trix.rds", compress="xz")
  } else{
    lin_reg_an_trix <- readRDS("data/lin_reg_an_trix.rds")
  }
  
  # Create TRIX maps
  p1 <- ggplot() + geom_raster(data=lin_reg_an_trix, aes(x=x, y=y, fill=factor(var))) + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_manual("Change", breaks=c(-1,0,1), labels = c("Decrease", "Stable", "Increase"), 
                      values=c("blue", "yellow", "red"), na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") + 
    geom_sf(data=countries10, fill="gray", color="black") + coord_sf() + 
    geom_polygon(data=df_realm, aes(x=long,y=lat, group=group), 
                 fill="transparent", colour="black") + ggtitle("a)") +
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.10))
  p2 <- ggplot() + geom_raster(data=lin_reg_an_trix, aes(x=x, y=y, fill=adjrsq)) + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("Adjusted R2", 
                         colours=colourtheme, na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") + 
    geom_sf(data=countries10, fill="gray", color="black") + coord_sf() + 
    geom_polygon(data=df_realm, aes(x=long,y=lat, group=group), 
                 fill="transparent", colour="black") + ggtitle("b)") +
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.10))
  
  # Issues with Legends require Grob workaround for Multiplot following: 
  # http://stackoverflow.com/questions/16255579/how-can-i-make-consistent-width-plots-in-ggplot-with-legends
  gA <- ggplotGrob(p1)
  gB <- ggplotGrob(p2)
    # The parts that differs in width
  leg1 <- grid::convertX(sum(with(gA$grobs[[15]], grobs[[1]]$widths)), "mm")
  leg2 <- grid::convertX(sum(with(gB$grobs[[15]], grobs[[1]]$widths)), "mm")
  # Add an empty column of "abs(diff(widths)) mm" width on the right of 
  # legend box for gA (the smaller legend box)
  gA$grobs[[15]] <- gtable::gtable_add_cols(gA$grobs[[15]], unit(abs(diff(c(leg1, leg2))), "mm"))
  # Combine the plots
  g = gridExtra::gtable_rbind(gA, gB, size = "max")
  grid.newpage()
  png(file=paste0(getwd(),"/figures/Figure_02.png"), 
      width=9, height=9, units="in", res=300)
  grid.draw(g)
  dev.off()
  rm(g,p1,p2,gA,gB,leg1,leg2)
}

## ---- an_mn_trix_realm ----

if(!file.exists(paste0(workdir, "data/summary_an_mn_trix.rds"))){
  # Load annual TRIX data
  an_trix_realm <- rasterFromXYZ(readRDS("data/an_trix_realm.rds"))
  
  # Calculate mean and sd of trix for each year
  an_mn_trix <- data.frame(mean=cellStats(an_trix_realm, stat=mean, na.rm=TRUE), 
                           sd=cellStats(an_trix_realm, stat=sd, na.rm=TRUE),
                           year=2002:2015); rm(an_trix_realm)
  # Write data to file
  saveRDS(an_mn_trix, file="data/summary_an_mn_trix.rds", compress="xz"); rm(an_mn_trix)
}

if(!file.exists("data/summary_an_trix_realm.rds")){
  # Load annual TRIX data
  an_trix_realm <- rasterFromXYZ(readRDS("data/an_trix_realm.rds"))
  
  # Calculate mean and sd of trix for each realm and each year
  mn_trix <- extract(an_trix_realm, sp_realm, fun=mean, df=TRUE, na.rm=TRUE)
  colnames(mn_trix) <- c("ID", 2002:2015)
  mn_trix <- tidyr::gather(mn_trix, "year", "mn_chl", -c(ID))
  sd_trix <- extract(an_trix_realm, sp_realm, fun=sd, df=TRUE, na.rm=TRUE)
  colnames(sd_trix) <- c("ID", 2002:2015)
  sd_trix <- tidyr::gather(sd_trix, "year", "sd_chl", -c(ID))
  sum_trix_realm <- merge(mn_trix, sd_trix, by=c("ID", "year")); rm(mn_trix, sd_trix)
  sum_trix_realm <- plyr::join(sum_trix_realm, data.frame(ID=c(1:length(names(sp_realm))), 
                                                          realm = names(sp_realm)), by="ID")
  colnames(sum_trix_realm) <- c("ID", "year", "mean", "sd", "realm")
  
  # Write Time Series Data to File
  saveRDS(sum_trix_realm, file="data/summary_an_trix_realm.rds", compress="xz")
  rm(sum_trix_realm, an_trix_realm)
}

## ---- an_mn_trix_lin_reg ----

# Load annual mean TRIX data
an_mn_trix <- readRDS("data/summary_an_mn_trix.rds")

# Run linear regression
m_trix <- lm(mean ~ year, data=an_mn_trix)

# Load summary of annual TRIX data per marine realm
sum_trix_realm <- readRDS("data/summary_an_trix_realm.rds")

# Create data frame for results of linear model
lin_reg_trix <- data.frame()

# Run linear regression between each explanatory variable and the different environmental impacts
for (i in 1:length(unique(sum_trix_realm$realm))){
  trix_sub <- subset(sum_trix_realm, sum_trix_realm$realm == unique(sum_trix_realm$realm)[i])
  trix_sub$year <- as.numeric(trix_sub$year)
  # Run individual model
  m <- lm(trix_sub$mean ~ trix_sub$year)
  
  # Extract R2 value and save to dataframe
  lin_reg_trix[i,1] <- unique(sum_trix_realm$realm)[i]
  lin_reg_trix[i,2] <- summary(m)$adj.r.squared
  lin_reg_trix[i,3] <- summary(m)$coefficients[2]
  lin_reg_trix[i,4] <- summary(m)$coefficients[8]
}
rm(trix_sub, m)

#Define names of dataframe
colnames(lin_reg_trix) <- c("Realm", "Adj. R2", "Intercept", "p-value")

#Produce table
pander::pandoc.table(lin_reg_trix, round=4, split.tables=Inf)

## ---- figX ----
if(!file.exists("figures/Figure_05.png")){
  
  # Load annual mean chlorophyll data
  an_mn_chl <- readRDS("data/summary_an_mn_chl.rds")
  
  # Load annual mean WOD data
  an_mn_nit <- readRDS("data/summary_an_mn_nit.rds")
  an_mn_pho <- readRDS("data/summary_an_mn_pho.rds")
  an_mn_dO2 <- readRDS("data/summary_an_mn_dO2.rds")
  
  # Load annual mean TRIX data
  an_mn_trix <- readRDS("data/summary_an_mn_trix.rds")
  
  # Create plot
  an_mn_chl$var <- "chl"
  an_mn_trix$var <- "trix"
  an_mn_nit$var <- "nit"
  an_mn_pho$var <- "pho"
  an_mn_dO2$var <- "dO2"
  an_mn_data <- plyr::rbind.fill(an_mn_chl, an_mn_nit, an_mn_pho, an_mn_dO2, an_mn_trix)
  an_mn_data$var <- factor(an_mn_data$var, 
                           labels = c("Chl-a", "Nitrate", "Phosphate", "aD%O", "TRIX")) 
  
  ggplot(data = an_mn_data, aes(x = year, y = mean)) + 
    geom_point() + scale_shape(solid=TRUE) + 
    # geom_linerange(aes(ymin=mean+sd, ymax=mean-sd)) + # Lines are way to long!
    geom_smooth(method = "lm", se=TRUE, formula = y ~ x, colour="black") + 
    labs(x="Year", y="") + 
    scale_y_continuous() +
    scale_x_continuous(breaks=seq(2002,2015,2)) + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "bottom", colour="black", parse=TRUE) + 
    facet_wrap(~ var, scale="free_y", strip.position="left") + 
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.y = element_text(size=16))
  ggsave(paste0("figures/Figure_05.png"), width=12, height=8, units="in", dpi=300)
}

## ---- fig3 ----
if(!file.exists("figures/Figure_03.png")){
  
  ## Fig. 3. Time series of Mean trophic index (TRIX) for each realm from 2002-2015.
  
  # Load annual TRIX data
  # an_trix_realm <- stack("data/an_trix_realm.tif")
  
  # Load summary of annual TRIX data
  sum_trix_realm <- readRDS("data/summary_an_trix_realm.rds")
  head(sum_trix_realm)
  
  # Time series graph with line of TRIX for each realm
  ggplot(sum_trix_realm, aes(x = as.numeric(year), y = mean)) + 
    labs(list(x = "Year", y = "Trophic index")) + 
    geom_point(shape=1) + 
    geom_smooth(method = "lm", se=TRUE, formula = y ~ x, colour="black") + 
    #scale_y_continuous(limits=c(4,8), breaks=seq(4,8,1)) + 
    scale_y_continuous(name="TRIX") + 
    scale_x_continuous(name="Year", breaks=seq(2002,2015,2)) + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.x.npc = "left", label.y.npc = "top", parse=TRUE, colour="black") + 
    facet_wrap(~ realm, ncol=3) + theme(strip.text.x = element_text(size=12, face="bold"),
                                        strip.background = element_blank(), 
                                        axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0("figures/Figure_03.png"), width=9, height=12, units="in", dpi=300)
}

## ---- figS11 ----

if(!file.exists("figures/SupplementaryFigure_11.png")){
  # Read trix
  an_trix_df <- readRDS("data/an_trix_realm.rds")
  years <- 2002:2015
  an_trix_long <- tidyr::gather(an_trix_df[,c(1:10)], "year", "trix",-c(x,y))
  colnames(an_trix_long) <- c("x", "y", "year", "trix")
  
  # Create TRIX maps
  ggplot() + 
    geom_sf(data=countries10, fill = "gray", color="black") + 
    geom_raster(data=an_trix_long, aes(x=x, y=y, fill=trix)) + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("TRIX", limits=c(2,9), breaks=seq(2,9,1), 
                         colours=trix_colours(7), na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") +
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines"))
  ggsave("figures/SupplementaryFigure_11.png", width=10, height=12, units="in", dpi=300)
}

## ---- figS12 ----

if(!file.exists("figures/SupplementaryFigure_12.png")){
  # Read trix
  an_trix_df <- readRDS("data/an_trix_realm.rds")
  years <- 2002:2015
  an_trix_long <- tidyr::gather(an_trix_df[,c(1,2,11:ncol(an_trix_df))], "year", "trix", -c(x,y))
  colnames(an_trix_long) <- c("x", "y", "year", "trix")
  
  # Create TRIX maps
  ggplot() + 
    geom_sf(data=countries10, fill = "gray", color="black") + 
    geom_raster(data=an_trix_long, aes(x,y, fill=trix)) + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    scale_fill_gradientn("TRIX", limits=c(2,9), breaks=seq(2,9,1), 
                         colours=trix_colours(7), na.value="transparent") +
    labs(x = "Longitude", y = "Latitude") +
    geom_polygon(data=df_realm, aes(long,lat, group=group), fill="transparent", colour="black") + 
    facet_wrap(~year, ncol=2) + coord_sf() + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines"))
  ggsave("figures/SupplementaryFigure_12.png", width=10, height=9, units="in", dpi=300)
}

## ---- halpern_data ----

## Table 2. Linear regression of each environmental variable and TRIX with transformed 
## direct place-based impact of shipping traffic frequency, organic pollution quantity, 
## inorganic pollution quantity, artisanal fishing rates, the number of invasive species 
## and ocean based pollution extracted from Halpern et al. [-@Halpern2008].
if(!file.exists("data/halpern_realm.rds")){
  if(!file.exists("data/halpern_realm_mask.rds")){
    halpern_all <- stack(lapply(list.files(path = paste0(filedir, "/Halpern/2008/"), 
                                    pattern=".tif", full.names=TRUE), FUN=function(x) raster(x)))
    #halpern_2008 <- extent(-18040095, 18040134, -9020047, 9020067)
    
    #halpern_2013 <- lapply(list.files(path = paste0(filedir, "/Halpern/2013/"), 
    #                                  pattern=".tif", full.names=TRUE), FUN=function(x) raster(x))
    
    #halpern_all <- stack(halpern_2008, halpern_2013)
    
    # Set projection of Halpern data (Molweiler)
    projection(halpern_all) <- CRS("+init=epsg:3857")
    
    # Transform sp_realm to Halpern data projection
    sp_realm_halpern <- spTransform(marinerealms, crs(halpern_all))
    
    # Mask Halpern data by sp_realm
    halpern_realm <- mask(halpern_all, sp_realm_halpern); rm(sp_realm_halpern, halpern_all)
    halpern_realm <- as.data.frame(rasterToPoints(halpern_realm))
    saveRDS(halpern_realm, "data/halpern_realm_mask.rds", compress="xz")
    } else{
    halpern_realm <- rasterFromXYZ(readRDS("data/halpern_realm_mask.rds"))
  }
  
  # Change projection of Halpern data
  if(!file.exists("data/halpern_realm_proj.rds")){
    beginCluster()
    projection(halpern_realm) <- CRS("+init=epsg:3857")
    halpern_realm <- projectRaster(halpern_realm, crs=crs.wgs84)
    halpern_realm <- as.data.frame(rasterToPoints(halpern_realm))
    saveRDS(halpern_realm, "data/halpern_realm_proj.rds", compress="xz")
    endCluster()
  } else{
    halpern_realm <- readRDS("data/halpern_realm_proj.rds")
    halpern_realm <- as.data.frame(rasterToPoints(halpern_realm))
  }
  
  trix_realm <- rasterFromXYZ(readRDS(paste0(workdir, "data/trix_REALM.rds")))
  
  # Need to be of same spatial resolution
  beginCluster()
  halpern_realm <- resample(halpern_realm, trix_realm)
  endCluster()
  halpern_realm <- as.data.frame(rasterToPoints(halpern_realm))
  saveRDS(halpern_realm, file="data/halpern_realm.rds", compress="xz")
}

## ---- figS13 ----

if(!file.exists("data/halpern_trix_realm_df.rds")){
  # Read halpern data from file
  halpern_realm <- readRDS("data/halpern_realm.rds")
  halpern_realm <- rasterFromXYZ(halpern_realm)
  
  # Read env data
  ma_env_realm <- readRDS("data/ma_env_realm.rds")
  ma_env_realm <- rasterFromXYZ(ma_env_realm)
  # Change to environmental conditions of 2008!
  
  # Read trix data
  trix_2008 <- readRDS("data/an_trix_realm.rds") %>% select(x,y, `2008`)
  trix_2008 <- rasterFromXYZ(trix_2008)
  
  # Read realm data
  r_realm <- rasterize(sp_realm, ma_env_realm)

  # Stack data
  halpern_trix_realm <- stack(ma_env_realm, trix_2008, halpern_realm, r_realm)
  
  # Convert Halpern data to dataframe
  halpern_trix_realm_df <- data.frame(rasterToPoints(halpern_trix_realm))
  rm(ma_env_realm, trix_2008, halpern_realm, halpern_trix_realm, r_realm)
  colnames(halpern_trix_realm_df) <- c("x", "y", "chl", "nit", "pho", "dO2", "trix", 
                                       "fish", "inorg", "invas", "organ", "pollu", "shipp", "realm")
  
  # Write Halpern Data to File
  saveRDS(halpern_trix_realm_df, file="data/halpern_trix_realm_df.rds", compress="xz")
} else{
  halpern_trix_realm_df <- readRDS("data/halpern_trix_realm_df.rds")
}

if(!file.exists("figures/SupplementaryFigure_13.png")){
  # Adapt data for plotting
  halpern_trix_realm_df <- halpern_trix_realm_df[,c("x", "y", "organ", "pollu", "shipp")]
  halpern_trix_realm_df <- halpern_trix_realm_df[apply(halpern_trix_realm_df[,c("organ", "pollu", "shipp")],1,function(x)any(!is.na(x))),]
  halpern_trix_long <- tidyr::gather(halpern_trix_realm_df, "var", "impact", -c(x,y)); rm(halpern_trix_realm_df)
  halpern_trix_long$var <- factor(halpern_trix_long$var, 
                                  labels=c("Organic pollution", "Ocean based pollution", "Shipping traffic"))
  
  # Create plot of Halpern data, MEOW and country outline
  p1 <- ggplot() + 
    geom_sf(data=countries10, fill = "gray", color="black") + 
    geom_tile(data=halpern_trix_long, aes(x=x,y=y, fill=impact)) + coord_sf() + 
    scale_fill_gradientn(name="Impact", colours=colourtheme, na.value="transparent", 
                         limits=c(0,1), breaks=seq(0,1,0.1)) + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    labs(x = "Longitude", y = "Latitude") +
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    facet_wrap(~var, ncol=1) + 
    theme(strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank(),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines"))
  ggsave("figures/SupplementaryFigure_13.png", width=8, height=12, units="in", dpi=300)
}

## ---- table2 ----

# Open Halpern-TRIX data
halpern_trix_realm_df <- readRDS("data/halpern_trix_realm_df.rds")

# Create data frame for results of linear model
lin_reg_halpern <- data.frame()

# Run linear regression between
for (i in 1:5){
  for(j in 1:6){
    # Run individual model
    m <- lm(halpern_trix_realm_df[,3+i] ~ halpern_trix_realm_df[,8+j])
    
    # Extract R2 value and save to dataframe
    lin_reg_halpern[i,j] <- summary(m)$coefficients[8]
  }
}

# Create row headers
lin_reg_halpern$Variable <- c("Chl-a", "Nitrate", "Phosphate", "aD%O2", "TRIX")

#Re-order dataframe
# lin_reg_halpern <- lin_reg_halpern[,c(7,1,2,3,4,5,6)]
lin_reg_halpern <- lin_reg_halpern[,c(7,4,5,6)]

#Define names of dataframe
#colnames(lin_reg_halpern) <- c("Variable", "Fishing", "Inorganic", "Invasive", 
#                               "Organic", "Pollution", "Shipping")
colnames(lin_reg_halpern) <- c("Variable", "Organic", "Pollution", "Shipping")

#Produce table of reference conditions 
pander::pandoc.table(lin_reg_halpern, round=4, split.tables=Inf)

## ---- fig4 ----
if(!file.exists("figures/Figure_04.png")){
  library(ecodist)
  #library(ade4)
  library(vegan)
  
  halpern_trix_realm_df <- readRDS("data/halpern_trix_realm_df.rds")
  halpern_trix_realm_df <- halpern_trix_realm_df[,c("trix", "fish", "inorg", "invas", 
                                                    "organ", "pollu", "shipp")]
  halpern_trix_realm_df <- na.omit(halpern_trix_realm_df) 
  
  # Classify TRIX into trophic state categories
  halpern_trix_realm_df$health <- NA 
  halpern_trix_realm_df$health[halpern_trix_realm_df$trix >= 8] <- NA 
  halpern_trix_realm_df$health[halpern_trix_realm_df$trix < 2] <- NA 
  halpern_trix_realm_df$health[halpern_trix_realm_df$trix >= 2 & halpern_trix_realm_df$trix < 4] <- 1 
  halpern_trix_realm_df$health[halpern_trix_realm_df$trix >= 4 & halpern_trix_realm_df$trix < 5] <- 2 
  halpern_trix_realm_df$health[halpern_trix_realm_df$trix >= 5 & halpern_trix_realm_df$trix < 6] <- 3 
  halpern_trix_realm_df$health[halpern_trix_realm_df$trix >= 6 & halpern_trix_realm_df$trix < 8] <- 4
  
  halpernPCA <- vegan::rda(halpern_trix_realm_df[,c("fish", "inorg", "invas", 
                                                    "organ", "pollu", "shipp")], scale=FALSE)
  halpernPCA
  smry <- summary(halpernPCA)
  df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2
  df2  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
  
  # Plot of PCA
  ggplot(data=df1, aes(x=PC1, y=PC2)) + coord_fixed() + 
    labs(x="Comp1, Axis1", y="Comp2, Axis2") +
    geom_hline(yintercept=0, col="darkgrey") + 
    geom_vline(xintercept=0, col="darkgrey") + 
    geom_point(aes(col=factor(halpern_trix_realm_df$health))) +
    scale_colour_manual("Status", breaks=c(4,3,2,1), 
                        labels = c("Bad", "Poor","Moderate","Good and above"), 
                        values=c("blue", "cyan", "yellow", "red"), 
                        na.value="transparent")
  ggsave(paste0("figures/PCA_Figure_01.png"), width=12, height=9, units="in", dpi=300)
  
  ggplot(data=df1, aes(x=PC1, y=PC2)) + coord_fixed() + 
    labs(x="Comp1, Axis1", y="Comp2, Axis2") +
    geom_hline(yintercept=0, col="darkgrey") + 
    geom_vline(xintercept=0, col="darkgrey") + 
    geom_point(aes(col=factor(halpern_trix_realm_df$health))) +
    scale_colour_manual("Status", breaks=c(4,3,2,1), 
                        labels = c("Bad", "Poor","Moderate","Good and above"), 
                        values=c("blue", "cyan", "yellow", "red"), 
                        na.value="transparent") +
    geom_segment(data=df2, aes(x=0, xend=PC1, y=0, yend=PC2), 
                 color="red", arrow=arrow(length=unit(0.01,"npc"))) +
    geom_text(data=df2, 
              aes(x=PC1,y=PC2,label=rownames(df2),
                  hjust=0.5*(1-sign(PC1)),vjust=0.5*(1-sign(PC2))), 
              color="red", size=4)
  ggsave(paste0("figures/PCA_Figure_02.png"), width=12, height=9, units="in", dpi=300)
  
  ggplot(data=df1, aes(x=PC1, y=PC2)) + coord_fixed() + 
    labs(x="Comp1, Axis1", y="Comp2, Axis2") +
    geom_hline(yintercept=0, col="darkgrey") + 
    geom_vline(xintercept=0, col="darkgrey") + 
    geom_point(aes(col=halpern_trix_realm_df$trix)) +
    scale_colour_gradientn("TRIX", colours=colourtheme, 
                           na.value="transparent")
  ggsave(paste0("figures/Figure_04.png"), 
         width=12, height=9, units="in", dpi=300)
}

## ---- figX2 ----  
if(!file.exists("figures/Figure_X2.png")){
  
  # Fortify Halpern and TRIX data
  halpern_trix_realm_df <- readRDS("data/halpern_trix_realm_df.rds")
  halpern_trix_realm_df <- halpern_trix_realm_df[,c("trix", "organ", "pollu", "shipp", "realm")]
  halpern_trix_realm_df <- tidyr::gather(halpern_trix_realm_df, "var", "value", -c(trix, realm))
  halpern_trix_realm_df$var <- factor(halpern_trix_realm_df$var, 
                                      labels = c("Organic pollution", "Ocean based pollution", "Shipping")) 
  halpern_trix_realm_df$realm <-  factor(halpern_trix_realm_df$realm, 
                                         labels=c("Arctic", "Central Indo-Pacific", "Eastern Indo-Pacific", 
                                                  "Southern Ocean", "Temperate Australasia", 
                                                  "Temperate Northern Atlantic", "Temperate Northern Pacific", 
                                                  "Temperate South America", "Temperate Southern Africa",
                                                  "Tropical Atlantic", "Tropical Eastern Pacific", 
                                                  "Western Indo-Pacific"))
  # Summarise data  
  halpern_trix_realm_sum <- summarise(group_by(halpern_trix_realm_df, var, realm), 
                                      mn_trix=mean(trix, na.rm=TRUE), mn_value=mean(value, na.rm=TRUE),
                                      sd_trix=sd(trix, na.rm=TRUE), sd_value=sd(value, na.rm=TRUE))
  
  # Linear regression of Halpern data and TRIX
  # ggplot(data = halpern_trix_realm_df, aes(x = value, y = trix)) + 
  #    geom_point(shape=1) # Original plot with just points!
  p1 <- ggplot(data = halpern_trix_realm_sum, aes(x = mn_value, y = mn_trix)) + geom_point() + scale_shape(solid=TRUE) + 
    #geom_errorbar(aes(ymin=mn_trix+sd_trix, ymax=mn_trix-sd_trix), width=0.00001) + 
    #geom_errorbarh(aes(xmin=mn_value-sd_value, xmax=mn_value+sd_value), height=0.01) + 
    labs(x="", y="TRIX (2008)") + 
    geom_smooth(method="lm", formula = y ~ x, se=FALSE) + 
    scale_x_continuous(expand=c(0,0.02)) + 
    scale_y_continuous() + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.x.npc = "right", label.y.npc = "bottom", 
                 rr.digits = 3, coef.digits = 2, parse=TRUE) + 
    facet_wrap(~ var, ncol=1, scales="free_x", strip.position="bottom") + 
    theme(strip.text.x = element_text(size=16), 
          strip.background=element_blank(),
          strip.placement= "outside",
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines"))
  ggsave(paste0("figures/Figure_X2.png"), width=6, height=10, units="in", dpi=300)
}

## ---- figSX ----

if(!file.exists("figures/SupplementaryFigure_X.png")){
  
  # Read data and subset accordingly
  halpern_trix_realm_df <- readRDS("data/halpern_trix_realm_df.rds")
  halpern_trix_realm_df <- halpern_trix_realm_df[,c("chl", "pho", "nit", "dO2",
                                                    "organ", "pollu", "shipp")]
  
  # Only remove rows where all columns have NAs
  halpern_trix_realm_df <- halpern_trix_realm_df[apply(
    halpern_trix_realm_df[,c("organ", "pollu", "shipp")],1,function(x)any(!is.na(x))),] 
  
  # Format data for facet_grid
  halpern_trix_realm_df <- tidyr::gather(halpern_trix_realm_df, "var", "value",
                                         -c(chl, pho, nit, dO2))
  halpern_trix_realm_df$var <- factor(halpern_trix_realm_df$var, 
                                      labels=c("Organic pollution", 
                                               "Ocean based pollution", "Shipping"))
  halpern_trix_realm_df <- tidyr::gather(halpern_trix_realm_df, "expl", "value2",
                                         -c(var, value))
  halpern_trix_realm_df$expl <- factor(halpern_trix_realm_df$expl, 
                                       labels=c("Chlorophyll a", "Nitrate", "Phosphate", "aD%O"))
  
  # Remove unnecessary data points
  halpern_trix_realm_df <- na.omit(halpern_trix_realm_df)
  
  # Create plot
  p1 <- ggplot(data = halpern_trix_realm_df, aes(x = value, y = value2)) + 
    geom_point(shape=1) + labs(x="", y="") + 
    geom_smooth(method="lm", formula = y ~ x) + 
    scale_x_continuous() + 
    #scale_x_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.1)) + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 label.x.npc="right", label.y.npc="top", parse = TRUE) +
    facet_grid(expl ~ var, scales="free", switch="both") + 
    theme(strip.background=element_blank(), strip.placement = "outside", 
          strip.text.y = element_text(size=16),
          strip.text.x = element_text(size=16),
          panel.spacing.x=unit(1.5, "lines"),
          panel.spacing.y=unit(1, "lines"))
  ggsave(paste0("figures/SupplementaryFigure_X.png"), width=9, height=12, units="in", dpi=300)
}

## ---- sum_trix_mpa ----

if(!file.exists(paste0(workdir, "data/trix_realm_mpa.rds"))){
  if(!file.exists(paste0(workdir, "data/trix_mpa.rds")) | !file.exists(paste0(workdir, "data/trix_non_mpa.rds"))){
    if(!file.exists("data/WDPA_Mar2020_IUCN_Marine.rds")){
      if(!dir.exists(paste0(filedir, "WDPA/WDPA_Mar2020-shapefile"))){
        # Download geodatabase of PAs
        download.file("https://pp-import-production.s3.amazonaws.com/WDPA_Mar2020_Public.zip", 
                      destfile=paste0(filedir, "WDPA/WDPA_Mar2020_Public.zip"))
        
        # Unzip geodatabase
        unzip("WDPA_Mar2020_Public.zip")
      }
      
      # Read WDPA Shapefile
      wdpa <- readOGR(dsn=paste0(filedir, "WDPA/WDPA_Mar2020-shapefile"), 
                      layer="WDPA_Mar2020-shapefile-polygons", verbose=FALSE)
      
      # Only take Protected Areas that are within the IUCN categories
      wdpa_iucn <- wdpa[wdpa$IUCN_CAT %in% c("II", "III", "IV", "Ia", "Ib", "V", "VI"),]; rm(wdpa)
      
      # Select protected areas that have a size bigger than 0
      wdpa_spatial <- wdpa_iucn[wdpa_iucn$REP_AREA > 0,]; rm(wdpa_iucn)
      
      # Divide protected areas in marine and terrestrial protected areas
      mpa <- wdpa_spatial[wdpa_spatial$MARINE != 0,]
      pa <- wdpa_spatial[wdpa_spatial$MARINE == 1,]
      rm(wdpa_spatial)
      
      # Save shapefile as RDS
      saveRDS(mpa, "data/WDPA_Mar2020_IUCN_Marine.rds", compress="xz")
    }
    
    # Read mpa
    mpa <- readRDS("data/WDPA_Mar2020_IUCN_Marine.rds")
    
    # Read trix
    trix_realm <-readRDS("data/trix_9km_realm.rds")
    trix_realm <- rasterFromXYZ(trix_realm)
    
    # Mask trix by protected and non-protected areas
    trix_mpa <- mask(trix_realm, mpa)
    trix_non_mpa <- mask(trix_realm, mpa, inverse=TRUE); rm(trix_realm)
    trix_mpa <- as.data.frame(rasterToPoints(trix_mpa))
    trix_non_mpa <- as.data.frame(rasterToPoints(trix_non_mpa))
    saveRDS(trix_mpa, "data/trix_mpa.rds", compress="xz")
    saveRDS(trix_non_mpa, "data/trix_non_mpa.rds", compress="xz")
  } else{
    trix_mpa <- readRDS("data/trix_mpa.rds")
    trix_mpa <- rasterFromXYZ(trix_mpa)
    trix_non_mpa <- readRDS("data/trix_non_mpa.rds")
    trix_non_mpa <- rasterFromXYZ(trix_non_mpa)
  }
  
  # Read realm data
  r_realm <- rasterize(sp_realm, trix_mpa)
  
  trix_realm_mpa <- stack(trix_mpa, trix_non_mpa, r_realm); rm(trix_mpa, trix_non_mpa, r_realm)
  trix_realm_mpa <- as.data.frame(trix_realm_mpa)
  colnames(trix_realm_mpa) <- c("pro", "unp", "realm")
  
  # Only remove rows where all columns have NAs
  trix_realm_mpa <- trix_realm_mpa[apply(trix_realm_mpa,1,function(x)any(!is.na(x))),] 
  trix_realm_mpa$realm <-  factor(trix_realm_mpa$realm, 
                                  labels=c("Arctic", "Central Indo-Pacific", "Eastern Indo-Pacific", 
                                           "Southern Ocean", "Temperate Australasia", 
                                           "Temperate Northern Atlantic", "Temperate Northern Pacific", 
                                           "Temperate South America", "Temperate Southern Africa",
                                           "Tropical Atlantic", "Tropical Eastern Pacific", 
                                           "Western Indo-Pacific"))
  # labels=c("Arctic", "C. Indo-Pac.", "E. Indo-Pac.", "S. Ocean", "Temp. Austr.", 
  # "Temp. N. Atl.", "Temp. N. Pac.", "Temp. S. Amer.", "Temp. S. Afr.",
  # "Trop. Atl.", "Trop. E.-Pac.", "W. Indo-Pac.")
  
  # Convert into long format
  trix_realm_mpa <- tidyr::gather(trix_realm_mpa, status, trix, -realm)
  trix_realm_mpa$status <- factor(trix_realm_mpa$status, labels=c("Protected", "Unprotected"))
  trix_realm_mpa <- na.omit(trix_realm_mpa)
  
  # Save df_trix_protected to file
  saveRDS(trix_realm_mpa, file="data/trix_realm_mpa.rds", compress="xz")
} else{
  # Read file
  trix_realm_mpa <- readRDS("data/trix_realm_mpa.rds")
  trix_realm_mpa$realm <-  factor(trix_realm_mpa$realm, 
                                  labels=c("Arctic", "Central Indo-Pacific", "Eastern Indo-Pacific", 
                                           "Southern Ocean", "Temperate Australasia", 
                                           "Temperate Northern Atlantic", "Temperate Northern Pacific", 
                                           "Temperate South America", "Temperate Southern Africa",
                                           "Tropical Atlantic", "Tropical Eastern Pacific", 
                                           "Western Indo-Pacific"))
  library(dplyr)
  sum_trix_mpa <- summarise(group_by(trix_realm_mpa, realm, status), mean=mean(trix))
  
  # Read MPA files for calculations
  trix_mpa <- readRDS("data/trix_mpa.rds")
  trix_mpa <- rasterFromXYZ(trix_mpa)
  trix_non_mpa <- readRDS("data/trix_non_mpa.rds")
  trix_non_mpa <- rasterFromXYZ(trix_non_mpa)
  
  # Calculate percentage of area protected
  perc_mpa <- length(na.omit(getValues(trix_mpa)))/(length(na.omit(getValues(trix_mpa)))+length(na.omit(getValues(trix_non_mpa))))*100
  
  # Percentage of area protected by marine realm
  perc_mpa_realm <- list()
  perc_realm <- sapply(1:length(sp_realm), FUN=function(x) 
    length(na.omit(getValues(crop(trix_mpa, sp_realm[x]))))/
      (length(na.omit(getValues(crop(trix_mpa, sp_realm[x]))))+
         length(na.omit(getValues(crop(trix_non_mpa, sp_realm[x])))))*100)
  perc_mpa_realm$perc <- as.vector(perc_realm)
  perc_mpa_realm$zone <- names(sp_realm)
  perc_mpa_realm <- as.data.frame(perc_mpa_realm)
}

## ---- fig5 ----
if(!file.exists(paste0(workdir, "figures/Figure_08.png"))){
  
  trix_realm_mpa <- readRDS("data/trix_realm_mpa.rds")
  
  # Plot bargraph 
  p1 <- ggplot(data=trix_realm_mpa, aes(x=factor(status), y=trix, fill=factor(status)), na.rm=TRUE) + 
    labs(list(x = "Status", y = "TRIX")) + 
    geom_boxplot() + 
    scale_y_continuous(expand=c(0,0.5), limits=c(0,10), breaks=seq(0,10,2)) + 
    scale_fill_manual(values = c("grey40", "grey90")) +
    geom_signif(comparisons = list(c("Protected", "Unprotected")), 
                map_signif_level=TRUE, test="t.test") + 
    facet_wrap(~ realm, ncol=3) + 
    theme(legend.position = "none", 
          strip.text.x = element_text(size=12, face="bold"),
          strip.background = element_blank())
  ggsave(paste0("figures/Figure_08.png"), width=9, height=12, units="in", dpi=300)
}

## ---- obis_trix ----

if(!file.exists(paste0(workdir, "data/obis_1deg_trix.rds"))){
  # Download OBIS summary layers with 0.1 degree resolution as shapefile
  # from http://www.iobis.org/geoserver/OBIS/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=OBIS:summaries&viewparams=table:map1deg_with_geom&outputformat=shape-zip"
  
  # Load obis summary data
  obis_summaries <- readOGR(paste0(filedir, "/OBIS/summaries_1deg"), layer="summaries", verbose=FALSE)
  
  # Get trix realm data
  trix_realm <- raster("data/trix_REALM.tif")
  ma_env <- stack("data/ma_env_REALM.tif")
  
  # Calculate mean TRIX for each polygon
  obis_trix <- extract(trix_realm, obis_summaries, fun=mean, na.rm = TRUE); rm(trix_realm)
  obis_env <- extract(ma_env, obis_summaries, fun=mean, na.rm = TRUE); rm(trix_realm)
  
  # Add TRIX data to SPDF
  obis_summaries@data <- cbind(obis_summaries@data, obis_trix)
  obis_summaries@data <- cbind(obis_summaries@data, obis_env)
  
  # Save new SPDF to file
  saveRDS(obis_summaries, file="data/obis_1deg_trix.rds")
}

## ---- obis_proper_format ----

if(!file.exists(paste0(workdir, "data/obis_1deg_trix.rds"))){
  
  # Read OBIS Summaries file
  obis_summaries <- readRDS("data/obis_1deg_trix.rds")
  
  # Prepare data for plotting
  obis_data <- fortify(obis_summaries, region="cscode")
  summary(obis_summaries@data$obis_trix)
  obis_summaries[obis_summaries@data$obis_trix == "NaN"] <- NA
  obis_summaries[is.na(obis_summaries@data$obis_trix)] <- NA
  obis_data <- merge(obis_data, obis_summaries@data, by.x="id", by.y="cscode")
  
  # Set variables to NA, where TRIX is NA
  obis_data <- obis_data[,c(1:13, ncol(obis_data))]
  colnames(obis_data) <- c("id", "long", "lat", "order", "hole", "piece", "group", "n", "s", 
                           "number_of_", "shannon", "simpson", "es", "obis_trix") 
  
  # Save obis data as .rds
  saveRDS(obis_data, "data/obis_1deg_trix.rds", compress="xz")
} else{
  obis_data <- readRDS("data/obis_1deg_trix.rds")
}

## ---- fig6 ----
if(!file.exists(paste0(workdir, "figures/Figure_06.png"))){
  
  # Read OBIS Summaries file
  obis_data <- readRDS("data/obis_1deg_trix.rds")
  
  # Subset obis data
  obis_data <- obis_data[,c("s", "number_of_", "shannon", "obis_trix")]
  
  # Create linear model
  var <- c("n", "s", "number_of_")
  m_obis <- lapply(1:length(var), FUN=function(x) lm(obis_trix ~ obis_data[,x], data=obis_data))
  
  # Change data to long format
  obis_data <- tidyr::gather(obis_data, "var", "value", -obis_trix)
  obis_data$var <- factor(obis_data$var, labels = c("Species richness", 
                                                    "No. of phyla", "Shannon index"))
  
  #Run linear regression
  ggplot(data = obis_data, aes(x = obis_trix, y = value)) + 
    labs(x="Mean TRIX", y="") + geom_point(shape=1) + 
    geom_smooth(method = "lm", formula = y ~ x) + 
    scale_x_continuous(expand=c(0.01,0.01), limits=c(2,9), breaks=seq(2,9,1)) + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 colour= "black", parse = TRUE) + 
    facet_wrap(~ var, ncol=1, scales = "free_y", strip.position="left") + 
    theme(strip.background = element_blank(), 
          strip.text.y = element_text(size=16), 
          strip.placement= "outside")
  ggsave(paste0("figures/Figure_06.png"), width=6, height=10, units="in", dpi=300)
}

## ---- figS14----
if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_14.png"))){
  
  obis_data <- readRDS("data/obis_1deg_trix.rds")
  
  #Plot species richness, number of phyla, shannon index
  p1 <- ggplot() + 
    geom_polygon(data=obis_data, aes(x=long, y=lat, group=group, fill=s)) + 
    scale_fill_gradientn(name="Species richness", colours=colourtheme, 
                         trans="log10", na.value="transparent") + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    labs(x = "Longitude", y = "Latitude") + ggtitle("a)") + 
    geom_sf(data=countries10, fill="gray", color="black") + coord_sf() + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.12))
  p2 <- ggplot() + 
    geom_polygon(data=obis_data, aes(x=long, y=lat, group=group, fill=number_of_)) + 
    scale_fill_gradientn(name="No. of phyla", colours=colourtheme, 
                         trans="log10", na.value="transparent") + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    labs(x = "Longitude", y = "Latitude") + ggtitle("b)") + 
    geom_sf(data=countries10, fill="gray", color="black") + coord_sf() + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.12))
  p3 <- ggplot() + 
    geom_polygon(data=obis_data, aes(x=long, y=lat, group=group, fill=shannon)) + 
    scale_fill_gradientn(name="Shannon index", colours=colourtheme, 
                         trans="log10", na.value="transparent") + 
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
    labs(x = "Longitude", y = "Latitude") + ggtitle("c)") + 
    geom_sf(data=countries10, fill="gray", color="black") + coord_sf() + 
    geom_polygon(data=df_realm, aes(long,lat, group=group), 
                 fill="transparent", colour="black") + 
    theme(legend.key = element_blank(),
          plot.title=element_text(size=18, hjust=-0.12))
  
  #Convert to plottable grobs for GridExtra
  gA <- ggplotGrob(p1)
  gB <- ggplotGrob(p2)
  gC <- ggplotGrob(p3)
  
  # Issues with Legends require Grob workaround for Multiplot following: 
  # http://stackoverflow.com/questions/17462504/align-edges-of-ggplot-chloropleth-legend-title-varies?rq=1
  grid.newpage()
  png(file=paste(getwd(),"/figures/SupplementaryFigure_14.png",sep=""), 
      width=8, height=12, units="in", res=300)
  g <- rbind(gA,gB,gC, size = "first")
  for (i in which(g$layout$name == "guide-box")) {
    g$grobs[[i]] <- g$grobs[[i]]$grobs[[1]]
  }
  grid.draw(g)
  dev.off()
  rm(g,p1,p2,p3,gA,gB,gC)
}

## ---- figSX2 ----

#if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_X2.png"))){

# Read data and subset accordingly
#  obis_data <- readRDS("data/obis_1deg_trix.rds")
#  obis_data <- obis_data[,c("chl", "pho", "nit", "dO2", "s", "number_of_", "shannon")]

# Only remove rows where all columns have NAs
#  obis_data <- obis_data[apply(obis_data[,c("chl", "pho", "nit", "dO2")],1,function(x)any(!is.na(x))),] 

# Format data for facet_grid
#  obis_data <- tidyr::gather(obis_data, "var", "value",
#                                  -c(chl, pho, nit, dO2))
#  obis_data$var <- factor(obis_data$var,
#                          labels=c("Species richness", "Number of Phyla", "Shannon index"))
#  obis_data <- tidyr::gather(obis_data, "expl", "value2", -c(var, value))
#  obis_data$expl <- factor(obis_data$expl, 
#                           labels=c("Chlorophyll a", "Nitrate", "Phosphate", "aD%O"))

# Remove unnecessary data points
#  obis_data <- na.omit(obis_data)

# Create plot
#  p1 <- ggplot(data = obis_data, aes(x = value, y = value2)) + 
#    geom_point(shape=1) + labs(x="", y="") + 
#    geom_smooth(method="lm", formula = y ~ x) + 
#    scale_x_continuous() + 
#    #scale_x_continuous(expand=c(0,0.1), limits=c(0,1), breaks=seq(0,1,0.1)) + 
#    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                label.x.npc="right", label.y.npc="top", 
#                 rr.digits = 3, coef.digits = 2, parse = TRUE) +
##    facet_grid(expl ~ var, scales="free", switch="both") + 
#    theme(strip.background=element_blank(), strip.placement = "outside", 
#          strip.text.y = element_text(size=16),
#          strip.text.x = element_text(size=16),
#          panel.spacing.x=unit(1.5, "lines"),
#          panel.spacing.y=unit(1, "lines"))
#  ggsave(paste0("figures/SupplementaryFigure_X2.png"), width=9, height=12, units="in", dpi=300)
#}

## ---- fao_data ----

if(!file.exists(paste0(workdir, "data/rfisheries_fao_all_2002_2015.rds"))){
  # Capture data downloaded from http://www.fao.org/fishery/statistics/global-capture-production/en
  fao_capture <- read.csv("Capture_2017.1.1/TS_FI_CAPTURE.csv")
  # Production data downloaded from http://www.fao.org/fishery/statistics/global-production/en
  #fao_production <- read.csv("GlobalProduction_2017.1.1/TS_FI_PRODUCTION.csv")
  
  # Get country codes
  fao_countries <- read.csv("Capture_2017.1.1/CL_FI_COUNTRY_GROUPS.csv")
  fao_countries <- fao_countries[,c("UN_CODE", "Name_en")]
  
  fao_sum <- summarise(group_by(fao_capture, Country, Year), total=sum(Quantity)); rm(fao_capture, fao_production)
  fao_year <- fao_sum[fao_sum$Year %in% c(2002:2015),]
  
  fao_all <- merge(fao_year, fao_countries, by.x = "Country", by.y="UN_CODE")
  colnames(fao_all) <- c("UN_CODE", "year", "catch", "country")
  fao_all <- fao_all[,c("year", "catch", "country")]
  
  saveRDS(fao_all, "data/rfisheries_fao_all_2002_2015.rds", compress="xz")
}

## ---- trix_iso ----

if(!file.exists(paste0(workdir, "data/an_mean_trix_iso3.rds"))){
  # Open shapefile of EEZ boundaries as Spatial Polygons Data Frame
  data(eez, package="geodat")
  # EEZ boundary shapefile (World EEZ v8) has been downloaded manually 
  # from Marineregions.org on the 19/09/2015
  
  # Create vector with names of sovereign states
  names_iso3 <- levels(eez$ISO_3digit); rm(eez)
  
  if(!file.exists(paste0(workdir, "data/sp_eez_iso3.rds"))){
    
    # Join EEZs into Sovereign shapes!
    sp_iso3 <- gUnaryUnion(sp_eez, sp_eez$ISO_3digit)
    #Alternatively:
    # sp_iso3 <- maptools::unionSpatialPolygons(sp_eez, sp_eez$ISO_3digit)
    
    # Write sovereign shapefile to file
    saveRDS(sp_iso3, "data/sp_eez_iso3.rds")
  }
  sp_iso3 <- readRDS("data/sp_eez_iso3.rds")
  
  # Load annual TRIX data
  an_trix_realm <- rasterFromXYZ(readRDS("data/an_trix_realm.rds"))
  
  # Extract TRIX data by country and year
  an_mean_trix_iso3 <- extract(an_trix_realm, sp_iso3, fun=mean, df=TRUE, na.rm=TRUE);   rm(sp_iso3)
  an_mean_trix_iso3 <- cbind(names_iso3, an_mean_trix_iso3)
  colnames(an_mean_trix_iso3) <- c("country", "ID", 2002:2015)
  an_mean_trix_iso3 <- tidyr::gather(an_mean_trix_iso3, year, an_mn_trix, -c(country, ID))
  
  # Write Annual Mean Trix to File
  saveRDS(an_mean_trix_iso3, "data/an_mean_trix_iso3.rds", compress="xz")
}

## ---- fao_trix_lm ----

# Read FAO data from file
fao_all <- readRDS("data/rfisheries_fao_all_2002_2015.rds")

# Load data from file
an_mean_trix_iso3 <- readRDS("data/an_mean_trix_iso3.rds")

# Merge TRIX and Fisheries data
trix_fao_year <- merge(fao_all, an_mean_trix_iso3, by = c("country", "year"))
rm(an_mean_trix_iso3)

trix_fao_year$log_catch <- log10(trix_fao_year$catch) 
trix_fao_year$log_catch[trix_fao_year$log_catch == -Inf] <- 0

# Calculate model of total catch per year!
fao_year <- summarise(group_by(trix_fao_year, year), total=sum(catch))
m_landing <- lm(total ~ year, data=fao_year)
#summary(m_landing)
#plot(total ~ year, data=fao_year)

# Calculate model for trix and catch
m_fao <- lm(an_mn_trix ~ log_catch,  data=trix_fao_year)

#Produce table´
summary(m_fao)

## ---- fig7 ----
if(!file.exists(paste0(workdir, "figures/Figure_07.png"))){
  
  # Plot scatterplot with regression 
  ggplot(data = trix_fao_year, aes(x = log_catch, y = an_mn_trix, col=year)) + 
    scale_x_continuous(name="log Annual landings (t)", limits=c(0,8), breaks=seq(0,8,1)) + 
    scale_y_continuous(name="Annual mean TRIX", limits=c(4,9), breaks=seq(4,9,1)) + 
    scale_colour_gradientn(name="Year", breaks=seq(2002,2015,2), colours=colourtheme) + 
    geom_point(shape=1) + geom_smooth(method = "lm", se=TRUE, formula = y ~ x) + 
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                 colour="black", parse=TRUE, label.x.npc = "left", label.y.npc="top")
  ggsave(paste0("figures/Figure_07.png"), width=8, height=6, units="in", dpi=300)
}
## ---- figS16 ----
#if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_16.png"))){

# Read FAO data from file
#fao_all <- readRDS("data/rfisheries_fao_all_2002_2015.rds")
#fao_all$year <- as.factor(fao_all$year)
#fao_all <- reshape2::dcast(fao_all, country ~ year, value.var="catch", sum)
#colnames(fao_all) <- c("country", years)

# Read trix per ISO
#an_mean_trix_iso3 <- readRDS("data/an_mean_trix_iso3.rds")
#an_mean_trix_iso3$year <- as.factor(an_mean_trix_iso3$year)
#an_mean_trix_iso3 <- tidyr::spread(an_mean_trix_iso3, year, an_mn_trix)
#colnames(an_mean_trix_iso3) <- c("country", "trix2002", "trix2003", "trix2004", "trix2005", "trix2006", 
#                                 "trix2007", "trix2008", "trix2009", "trix2010", "trix2011", 
#                                 "trix2012", "trix2013", "trix2014", "trix2015")

# Add data to EEZ shapefile
#if(!file.exists("data/sp_eez_df.csv")){
#  load(paste0(workdir, "/data/sp_eez.rda"))
# convert polygon into a dataframe
#  sp_eez_df <- fortify(sp_eez)
# merge the "fortified" data with the data from our spatial object
#  sp_eez_df <- merge(sp_eez_df, sp_eez@data, by.x = "id", by.y="ID")
#  write.csv(sp_eez_df, "data/sp_eez_df.csv")
#} else{
#  sp_eez_df <- read.csv("data/sp_eez_df.csv")
#}

# Add fishery data for plotting
#sp_eez_df <- merge(sp_eez_df, fao_all, by.x="ISO_3digit", by.y="country"); rm(fao_all)

# Add trix data for plotting
# sp_eez_df <- merge(sp_eez_df, an_mean_trix_iso3, by.x="country", by.y="country")

# Convert into long format
#sp_eez_df <- sp_eez_df[,c("id", "long", "lat", "order", "hole", "piece", "group", "2002", "2003",
#                         "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012",
#                          "2013", "2014", "2015")]
#sp_eez1 <- tidyr::gather(sp_eez_df[,c(1:15)], "year", "fish", 
#                         -c(id, long, lat, order, hole, piece, group)); rm(sp_eez_df)

#Plot fishery landing per year
#p1 <- ggplot() + 
#  geom_polygon(data=sp_eez1, aes(x=long, y=lat, group=group, fill=fish), colour="black") + 
#  scale_fill_gradientn(name="Fishery landings (t)", colours=colourtheme, 
#                       na.value="transparent") + facet_wrap(~ year, ncol=2) + 
#  scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
#  scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
#  labs(x = "Longitude", y = "Latitude") + 
#  geom_sf(data=countries10, fill="gray", color="black") + coord_sf() +
#  theme(strip.text.x = element_text(size=12, face="bold"),
#      strip.background = element_blank(),
#      panel.spacing.x=unit(1.5, "lines"),
#      panel.spacing.y=unit(1, "lines")); rm(sp_eez1)
#ggsave(paste0("figures/SupplementaryFigure_16.png"), width=9, height=8, units="in", dpi=300); rm(p1)
#}

## ---- figS17 ----
#if(!file.exists(paste0(workdir, "figures/SupplementaryFigure_17.png"))){

# Read FAO data from file
#  fao_all <- readRDS("data/rfisheries_fao_all_2002_2015.rds")
#  fao_all$year <- as.factor(fao_all$year)
#  fao_all <- reshape2::dcast(fao_all, country ~ year, value.var="catch", sum)
#  colnames(fao_all) <- c("country", years)

# Read trix per ISO
#an_mean_trix_iso3 <- readRDS("data/an_mean_trix_iso3.rds")
#an_mean_trix_iso3$year <- as.factor(an_mean_trix_iso3$year)
#an_mean_trix_iso3 <- tidyr::spread(an_mean_trix_iso3, year, an_mn_trix)
#colnames(an_mean_trix_iso3) <- c("country", "trix2002", "trix2003", "trix2004", "trix2005", "trix2006", 
#                                 "trix2007", "trix2008", "trix2009", "trix2010", "trix2011", 
#                                 "trix2012", "trix2013", "trix2014", "trix2015")

# Add data to EEZ shapefile
#  if(!file.exists("data/sp_eez_df.csv")){
#    load(paste0(workdir, "/data/sp_eez.rda"))
#    # convert polygon into a dataframe
#    sp_eez_df <- fortify(sp_eez)
# merge the "fortified" data with the data from our spatial object
#    sp_eez_df <- merge(sp_eez_df, sp_eez@data, by.x = "id", by.y="ID")
#    write.csv(sp_eez_df, "data/sp_eez_df.csv")
#  } else{
#    sp_eez_df <- read.csv("data/sp_eez_df.csv")
#  }

# Add fishery data for plotting
#  sp_eez_df <- merge(sp_eez_df, fao_all, by.x="ISO_3digit", by.y="country"); rm(fao_all)

# Add trix data for plotting
# sp_eez_df <- merge(sp_eez_df, an_mean_trix_iso3, by.x="country", by.y="country")

# Convert into long format
#  sp_eez_df <- sp_eez_df[,c("id", "long", "lat", "order", "hole", "piece", "group", "2002", "2003",
#                            "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012",
#                            "2013", "2014", "2015")]
#  sp_eez2 <- tidyr::gather(sp_eez_df[,c(1:7,16:ncol(sp_eez_df))], "year", "fish", 
#                           -c(id, long, lat, order, hole, piece, group)); rm(sp_eez_df)

#Plot fishery landing per year
#  p2 <- ggplot() + 
#    geom_polygon(data=sp_eez2, aes(x=long, y=lat, group=group, fill=fish), colour="black") + 
#    scale_fill_gradientn(name="Fishery landings (t)", colours=colourtheme, 
#                         na.value="transparent") + facet_wrap(~ year, ncol=2) + 
#    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(-180,180,40)) + 
#    scale_y_continuous(expand=c(0.01,3), breaks=seq(-90,90,30)) +
#    labs(x = "Longitude", y = "Latitude") + 
#    geom_sf(data=countries10, fill="gray", color="black") + coord_sf() + 
#    theme(strip.text.x = element_text(size=12, face="bold"),
#          strip.background = element_blank(),
#         panel.spacing.x=unit(1.5, "lines"),
#          panel.spacing.y=unit(1, "lines")); rm(sp_eez2)
#  ggsave(paste0("figures/SupplementaryFigure_17.png"), width=9, height=8, units="in", dpi=300); rm(p2)
#}
