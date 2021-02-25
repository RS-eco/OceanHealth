#Function to calculate TRIX components
trix <- function(x, y, output=TRUE){
  r <- raster::rasterize(x=y, y=x)
  
  #Zonal stats for getting mean and stdev per polygon 
  mean_r <- raster::zonal(x, r, fun='mean', na.rm=T)
  std_r <- raster::zonal(x, r, fun='sd', na.rm=T)
  
  #Calculate upper limit (Mean + 2.5 * STD) per polygon
  ul_r <- mean_r[,2:5] + (2.5*std_r[,2:5])
  ul_r <- cbind(1:12, ul_r)
  colnames(ul_r)[1] <- c("id")
  ul_r <- subs(r, data.frame(ul_r), which=c(2,3,4,5))
  
  #Calculate lower limit (Mean - 2.5 * STD) per polygon
  ll_r <- mean_r[,2:5] - (2.5*std_r[,2:5])
  ll_r <- cbind(1:12, ll_r)
  colnames(ll_r)[1] <- c("id")
  ll_r <- subs(r, data.frame(ll_r), which=c(2,3,4,5))
  
  #Value above mean + 2.5 * STD
  ext_r <- x
  ext_r[ext_r > ul_r] <- NA
  ext_r <- raster::mask(x, ext_r)
  
  #Remove values below Mean - 2.5 * STD
  ext_r[ext_r < ll_r] <- NA
  ext_r <- raster::mask(x, ext_r)
  plot(ext_r)
  
  #Log-Transfom data
  log_r <- log(ext_r)
  
  #Calculate minimum & maximum
  max_r <- data.frame(raster::zonal(log_r, r, fun='max'))
  min_r <- data.frame(raster::zonal(log_r, r, fun='min'))
  ref_cond <- cbind(min_r, max_r)
  if(output==TRUE){
    saveRDS(ref_cond, "Results/ref_cond_REALM.rds", compress="xz")
  }
  
  #Turn min & max into raster
  min_r <- subs(r, min_r, which=c(2,3,4,5))
  max_r <- subs(r, max_r, which=c(2,3,4,5))
  
  #Calculate component values (Measured value - Lower limit/Upper limit - lower limit)
  comp_r <- ((log_r - min_r)/ (max_r - min_r))
  
  # Calculate sum of individual components
  comp_r <- calc(comp_r, fun=sum)
  
  #Calculation of TRIX
  trix <- 2.5*comp_r
  return(trix)
}
