############################################################################
#  2.4  create RR list data                                                
############################################################################
dat <- list()
dat$fyear <- fyear
dat$lyear <- lyear 
dat$fage <- fage
dat$lage <- lage
dat$tcatch <- tcatch
dat$h_low <- tlcatch
dat$h_up <- tucatch
dat$minesc <- minesc
dat$minrun <- minrun
# Inriver Run
dat$inriv <- kusko.data$In.river
dat$inriv_sd <- kusko.data$In.river.sd
# Radio telemetry data 
dat$mr_tr <- kusko.data$tmr
dat$mrsd_tr <- kusko.data$tmr.sd
# Sonar  data 
dat$sonar <- kusko.data$Sonar
dat$sonar_sd <- kusko.data$Sonar.sd
# MR data
dat$mr <- kusko.data$mr
dat$mrsd <- kusko.data$mr.sd
# MR above Aniak
dat$amr <- kusko.data$amr
dat$amrsd <- kusko.data$amr.sd
#MR Holitna
dat$hmr <- kusko.data$hmr
dat$hmrsd <- kusko.data$hmr.sd
dat$w_esc <- t(w.esc)
dat$a_esc <- t(a.esc)
dat$cpue <- t(c.cpue)
dat$testp <- t(c.prop)
# These are assumed CV for weir observation 
# 1:Kwethluk, 2:Tuluksuk, 3: Aniak Salmon, 4:George, 5:Koglukluk, 6:Tatluiksuk, 7:Takotna, 8:Pitka Salmon 
dat$cvw <- c(0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15)
# These are assumed CV for aerial for likelihood 
# 1:Kwethluk, 2:Tulukuk, 3:Kisaralik, 4:Aniak Salmon, 5:Kipchuk, 6:Anik
# 7:Holokuk, 8:Oskawalik, 9:Holitna, 10:Cheeneetnuk, 11:Gagaryah, 12:Pitka
# 13:Bear, 14: Pitka Salmon
dat$cva <- c(0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50)
# These are assumed CV for fishery cpue: Large, mix, restrict 
dat$cvf <- c(0.20,0.2,0.2)
# These are assumed cv for harvest 
dat$cvh <- c(0.05,0.05)
# This is a conrol vector for weir likelihood 
# 1:Kwethluk, 2:Tuluksuk, 3: Aniak Salmon, 4:George, 5:Koglukluk, 6:Tatluiksuk, 7:Takotna, 8:Pitka Salmon 
dat$wlike <- c(1,1,1,1,1,1,1,1)
# This is the conrol vector aerial for likelihood 
# 1:Kwethluk, 2:Tulukuk, 3:Kisaralik, 4:Aniak Salmon, 5:Kipchuk, 6:Anik
# 7:Holokuk, 8:Oskawalik, 9:Holitna, 10:Cheeneetnuk, 11:Gagaryah, 12:Pitka
# 13:Bear, 14: Pitka Salmon
dat$alike <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1)
# Fishery likelihood: Large, mix, restrict  
dat$flike <- c(1,1,1)
# 1, Lower MR, 2. Bethel Sonar, 3. Upper MR, 4. Above Aniak, 5. Holitna
dat$rlike <- c(1,1,1,0,0)
dat <- rapply( dat, f=function(x) ifelse(is.na(x),0,x), how="replace" )
