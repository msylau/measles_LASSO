


rm(list=ls())

require(geosphere)
require(glmnet)
require(dplyr)
require(tidyr)
require(foreach)
require(parallel)
require(zoo)
require(TSdist)
require(Metrics)
require(hydroGOF)

set.seed(1)


log_scale = "TRUE" # IF model It on log scale
off_set = 1

mc.cores = 8

country = "UK"

###

if (country=="US"){

	data_dir = "xx"
	result_dir = "xx"

	USUK_data = read.csv(file=paste(data_dir, "measlesUKUS.csv",sep=""))
	US_data = subset(USUK_data, USUK_data$country=="US")

	cases_urban = US_data %>% select(c("decimalYear","loc","cases")) %>% tidyr::spread(key=loc,value=cases) %>% rename(X=decimalYear)
	pop_urban = US_data %>% select(c("decimalYear","loc","pop")) %>% tidyr::spread(key=loc,value=pop) %>% rename(X=decimalYear)  
	births_urban = US_data %>% select(c("decimalYear","loc","rec")) %>% tidyr::spread(key=loc,value=rec) %>% rename(X=decimalYear) 

	
	larger_places = c(colnames(cases_urban)[-1][which(apply(pop_urban[,-1],2,FUN=median)>800000)]) 

	cases_urban = na.contiguous(cases_urban[,c("X",larger_places)]) 
	pop_urban = na.contiguous(pop_urban[,c("X",larger_places)])

	median_pop = apply(as.matrix(pop_urban[,-1]),2,FUN=median)
	time_in_year = cases_urban$X
	names_cities <- colnames(cases_urban)[-1] 

	### parameters ####
	###################

	Ta= 26 
	Tb=5*26 
	Tc=5*26 

	max_k_ahead = 26
	t_initial = max(Ta,Tb) + 1 

	model.versions.to.loop = c("V1_fullreduced","V1_fullreduced_with_births") 

}
###

if (country=="UK"){

	era = "prevac" #

	if(era=="prevac") data_dir = "xx"
	if(era=="vac") data_dir = "xx"

	result_dir = "xx"

	cases_urban <- read.csv(paste(data_dir,"cases_urban.csv",sep=""))

	cases_urban = cases_urban %>% dplyr::select(-c("City.of.London"))# drop some to avoid duplicated coordinates


	pop_urban <- read.csv(paste(data_dir,"pop_urban.csv",sep=""))[,colnames(cases_urban)]
	coordinates_urban <- read.csv(paste(data_dir,"coordinates_urban.csv",sep=""))[,colnames(cases_urban)]
	susceptibles_urban <- read.csv(paste(data_dir,"susceptibles_urban.csv",sep=""))[,colnames(cases_urban)]
	births_urban_by_intyear <- read.csv(paste(data_dir,"births_urban.csv",sep=""))[,colnames(cases_urban)]
	births_urban_by_intyear[,2:ncol(births_urban_by_intyear)] = births_urban_by_intyear[,2:ncol(births_urban_by_intyear)]/26 # biweek births
	births_urban_by_intyear = births_urban_by_intyear %>% rename(year=X)
	births_urban = data.frame(X=cases_urban$X, year=as.numeric(trunc(cases_urban$X)))
	births_urban = merge(births_urban, births_urban_by_intyear,by="year")
	births_urban = births_urban %>% select(-c("year"))


	if(era=="vac"){ 

		cases_urban = subset(cases_urban,cases_urban$X>70)
		pop_urban = subset(pop_urban,pop_urban$X>70)
		susceptibles_urban = subset(susceptibles_urban,susceptibles_urban$X>70)
		births_urban = subset(births_urban,births_urban$X>70)

	}


	time_in_year = cases_urban$X


	median_pop = apply(pop_urban[,-1],2,FUN=median)

	names_cities <- colnames(cases_urban)[-1] # the first column is "X": time or lat/long

	### parameters ####
	###################

	Ta= 5*26 
	Tb=5*26 
	Tc=5*26 

	max_k_ahead = 2*26
	t_initial = max(Ta,Tb) + 1 


	model.versions.to.loop = c("V1_fullreduced","V1_fullreduced_with_births") 


}

####



for (model.version in model.versions.to.loop){


	size.train = which(cases_urban[,1]==56.01923)- t_initial 



	places_to_fit = c(names_cities[median_pop>300000])



	sort(median_pop[places_to_fit])

	###################




	list.all.places = mclapply( places_to_fit, 
		function(test_place) { 

			pop_test_place = median_pop[which(names_cities==test_place)]



			if(country=="UK") dist_to_others = sapply(names_cities[-c(which(names_cities==test_place))], 
								FUN=function(x){distCosine(coordinates_urban[,test_place], coordinates_urban[,x])/1000}) # dist between test_place and other places



			if(model.version=="V1_fullreduced") X.v1.fullreduced = matrix(NA, nrow=nrow(cases_urban)-t_initial-max_k_ahead+1, ncol=1*Tb)
			if(model.version=="V1_fullreduced_with_births") X.v1.fullreduced.with.births = matrix(NA, nrow=nrow(cases_urban)-t_initial-max_k_ahead+1, ncol=1+Tb)

			if(log_scale==FALSE){
				Y_oneahead = cases_urban[,test_place][t_initial:(nrow(cases_urban)-max_k_ahead)]
				Y_twoahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+1]
				Y_threeahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+2]
				Y_fourahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+3]
				Y_fiveahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+4]
				Y_sixahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+5]
				Y_sevenahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+6]
				Y_eightahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+7]

				Y_9ahead = cases_urban[,test_place][t_initial:(nrow(cases_urban)-max_k_ahead)+8]
				Y_10ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+9]
				Y_11ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+10]
				Y_12ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+11]
				Y_13ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+12]
				Y_14ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+13]
				Y_15ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+14]
				Y_16ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+15]
				Y_17ahead = cases_urban[,test_place][t_initial:(nrow(cases_urban)-max_k_ahead)+16]
				Y_18ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+17]
				Y_19ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+18]
				Y_20ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+19]
				Y_21ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+20]
				Y_22ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+21]
				Y_23ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+22]
				Y_24ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+23]
				Y_25ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+24]
				Y_26ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+25]


				Y_27ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+26]
				Y_28ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+27]
				Y_29ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+28]
				Y_30ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+29]
				Y_31ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+30]
				Y_32ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+31]
				Y_33ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+32]
				Y_34ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+33]
				Y_35ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+34]
				Y_36ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+35]
				Y_37ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+36]
				Y_38ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+37]
				Y_39ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+38]
				Y_40ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+39]
				Y_41ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+40]
				Y_42ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+41]
				Y_43ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+42]
				Y_44ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+43]
				Y_45ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+44]
				Y_46ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+45]
				Y_47ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+46]
				Y_48ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+47]
				Y_49ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+48]
				Y_50ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+49]
				Y_51ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+50]
				Y_52ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+51]



			}

			

			if(log_scale==TRUE){


				Y_oneahead = cases_urban[,test_place][t_initial:(nrow(cases_urban)-max_k_ahead)]
				Y_twoahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+1]
				Y_threeahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+2]
				Y_fourahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+3]
				Y_fiveahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+4]
				Y_sixahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+5]
				Y_sevenahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+6]
				Y_eightahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+7]

				Y_9ahead = cases_urban[,test_place][t_initial:(nrow(cases_urban)-max_k_ahead)+8]
				Y_10ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+9]
				Y_11ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+10]
				Y_12ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+11]
				Y_13ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+12]
				Y_14ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+13]
				Y_15ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+14]
				Y_16ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+15]

				Y_17ahead = cases_urban[,test_place][t_initial:(nrow(cases_urban)-max_k_ahead)+16]
				Y_18ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+17]
				Y_19ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+18]
				Y_20ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+19]
				Y_21ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+20]
				Y_22ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+21]
				Y_23ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+22]
				Y_24ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+23]
				Y_25ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+24]
				Y_26ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+25]


				Y_27ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+26]
				Y_28ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+27]
				Y_29ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+28]
				Y_30ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+29]
				Y_31ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+30]
				Y_32ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+31]
				Y_33ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+32]
				Y_34ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+33]
				Y_35ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+34]
				Y_36ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+35]
				Y_37ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+36]
				Y_38ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+37]
				Y_39ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+38]
				Y_40ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+39]
				Y_41ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+40]
				Y_42ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+41]
				Y_43ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+42]
				Y_44ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+43]
				Y_45ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+44]
				Y_46ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+45]
				Y_47ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+46]
				Y_48ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+47]
				Y_49ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+48]
				Y_50ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+49]
				Y_51ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+50]
				Y_52ahead = cases_urban[,test_place][(t_initial:(nrow(cases_urban)-max_k_ahead))+51]


			}

		


	
	 		if(model.version=="V1_fullreduced"){

				X.v1.fullreduced = t(sapply(t_initial:(nrow(cases_urban)-max_k_ahead), 

					function(t){

						if(log_scale==FALSE) It_at_lags = cases_urban[,test_place][(t-Tb):(t-1)]

						if(log_scale==TRUE) It_at_lags = log(cases_urban[,test_place][(t-Tb):(t-1)]+off_set)


						local_I_lagged_predictors = It_at_lags # local predictors, length should be = Ta

						c(local_I_lagged_predictors)
					}
				))

			}


	 		if(model.version=="V1_fullreduced_with_births"){

				X.v1.fullreduced.with.births = t(sapply(t_initial:(nrow(cases_urban)-max_k_ahead), 

					function(t){

						if(log_scale==FALSE){
							It_at_lags = cases_urban[,test_place][(t-Tb):(t-1)]
							# Bt_at_lags =  (pop_urban[,test_place][(t-Tc):(t-1)]>350000)*births_urban[,test_place][(t-Tc):(t-1)]
							Bt_at_lags =  (median_pop[test_place]>300000)*births_urban[,test_place][(t-Tb):(t-1)]
						}


						if(log_scale==TRUE){

							It_at_lags = log(cases_urban[,test_place][(t-Tb):(t-1)]+off_set)

							Bt_at_lags =   log(mean(births_urban[,test_place][(t-Tc):(t-1)]+1))

						}



						local_I_lagged_predictors = It_at_lags # local predictors, length should be = Ta
						local_B_lagged_predictors = Bt_at_lags # local predictors, length should be = Ta

						c(local_B_lagged_predictors, local_I_lagged_predictors)
					}
				))

			}


			##



			if(model.version=="V1_fullreduced") colnames(X.v1.fullreduced) = c(paste("local_I_lagged_",Tb:1,sep=""))

			if(model.version=="V1_fullreduced_with_births") colnames(X.v1.fullreduced.with.births) = c(paste("births_lagged_",1:1,sep=""),paste("local_I_lagged_",Tb:1,sep=""))


			if(model.version=="V1_fullreduced") X = X.v1.fullreduced # chose a version of model to fit
			if(model.version=="V1_fullreduced_with_births") X = X.v1.fullreduced.with.births # chose a version of model to fit


		

			## store the covariate, observed and predicted incidence of places ##
			#####################################################################
			if (country=="UK") list.test.place = list(X=X, Y_oneahead=Y_oneahead, Y_twoahead=Y_twoahead,  Y_threeahead=Y_threeahead, Y_fourahead=Y_fourahead, Y_fiveahead=Y_fiveahead,Y_sixahead=Y_sixahead,Y_sevenahead=Y_sevenahead,Y_eightahead=Y_eightahead , 
				Y_9ahead=Y_9ahead ,Y_10ahead=Y_10ahead ,Y_11ahead=Y_11ahead ,Y_12ahead=Y_12ahead ,Y_13ahead=Y_13ahead ,Y_14ahead=Y_14ahead ,Y_15ahead=Y_15ahead ,Y_16ahead=Y_16ahead ,Y_17ahead=Y_17ahead ,Y_18ahead=Y_18ahead ,Y_19ahead=Y_19ahead ,Y_20ahead=Y_20ahead ,Y_21ahead=Y_21ahead ,Y_22ahead=Y_22ahead ,Y_23ahead=Y_23ahead ,Y_24ahead=Y_24ahead ,Y_25ahead=Y_25ahead ,Y_26ahead=Y_26ahead ,
				Y_27ahead=Y_27ahead ,
				Y_28ahead=Y_28ahead ,
				Y_29ahead=Y_29ahead ,
				Y_30ahead=Y_30ahead ,
				Y_31ahead=Y_31ahead ,
				Y_32ahead=Y_32ahead ,
				Y_33ahead=Y_33ahead ,
				Y_34ahead=Y_34ahead ,
				Y_35ahead=Y_35ahead ,
				Y_36ahead=Y_36ahead ,
				Y_37ahead=Y_37ahead ,
				Y_38ahead=Y_38ahead ,
				Y_39ahead=Y_39ahead ,
				Y_40ahead=Y_40ahead ,
				Y_41ahead=Y_41ahead ,
				Y_42ahead=Y_42ahead ,
				Y_43ahead=Y_43ahead ,
				Y_44ahead=Y_44ahead ,
				Y_45ahead=Y_45ahead ,
				Y_46ahead=Y_46ahead ,
				Y_47ahead=Y_47ahead ,
				Y_48ahead=Y_48ahead ,
				Y_49ahead=Y_49ahead ,
				Y_50ahead=Y_50ahead ,
				Y_51ahead=Y_51ahead ,
				Y_52ahead=Y_52ahead ,
				pop_test_place=pop_test_place, dist_to_others=dist_to_others, Ta=Ta, Tb=Tb, max_k_ahead=max_k_ahead, t_initial=t_initial,time_in_year=time_in_year, size.train=size.train)
			# list.all.places$test_place = list.test.place


			if (country=="US")  list.test.place = list(X=X, Y_oneahead=Y_oneahead, Y_twoahead=Y_twoahead,  Y_threeahead=Y_threeahead, Y_fourahead=Y_fourahead, Y_fiveahead=Y_fiveahead,Y_sixahead=Y_sixahead,Y_sevenahead=Y_sevenahead,Y_eightahead=Y_eightahead , 
				Y_9ahead=Y_9ahead ,Y_10ahead=Y_10ahead ,Y_11ahead=Y_11ahead ,Y_12ahead=Y_12ahead ,Y_13ahead=Y_13ahead ,Y_14ahead=Y_14ahead ,Y_15ahead=Y_15ahead ,Y_16ahead=Y_16ahead ,Y_17ahead=Y_17ahead ,Y_18ahead=Y_18ahead ,Y_19ahead=Y_19ahead ,Y_20ahead=Y_20ahead ,Y_21ahead=Y_21ahead ,Y_22ahead=Y_22ahead ,Y_23ahead=Y_23ahead ,Y_24ahead=Y_24ahead ,Y_25ahead=Y_25ahead ,Y_26ahead=Y_26ahead ,
				pop_test_place=pop_test_place, Ta=Ta, Tb=Tb, max_k_ahead=max_k_ahead, t_initial=t_initial,time_in_year=time_in_year, size.train=size.train)
			# list.all.places$test_place = list.test.place

			list.test.place

			}, mc.cores=mc.cores)


	names(list.all.places) = places_to_fit

	## training data ##
	##################



	X_all = list.all.places[[places_to_fit[1]]]$X[1:size.train,]

	Y_oneahead_all = list.all.places[[places_to_fit[1]]]$Y_oneahead[1:size.train]
	Y_twoahead_all = list.all.places[[places_to_fit[1]]]$Y_twoahead[1:size.train]
	Y_threeahead_all = list.all.places[[places_to_fit[1]]]$Y_threeahead[1:size.train]
	Y_fourahead_all = list.all.places[[places_to_fit[1]]]$Y_fourahead[1:size.train]
	Y_fiveahead_all = list.all.places[[places_to_fit[1]]]$Y_fiveahead[1:size.train]
	Y_sixahead_all = list.all.places[[places_to_fit[1]]]$Y_sixahead[1:size.train]
	Y_sevenahead_all = list.all.places[[places_to_fit[1]]]$Y_sevenahead[1:size.train]
	Y_eightahead_all = list.all.places[[places_to_fit[1]]]$Y_eightahead[1:size.train]

	Y_9ahead_all = list.all.places[[places_to_fit[1]]]$Y_9ahead[1:size.train]
	Y_10ahead_all = list.all.places[[places_to_fit[1]]]$Y_10ahead[1:size.train]
	Y_11ahead_all = list.all.places[[places_to_fit[1]]]$Y_11ahead[1:size.train]
	Y_12ahead_all = list.all.places[[places_to_fit[1]]]$Y_12ahead[1:size.train]
	Y_13ahead_all = list.all.places[[places_to_fit[1]]]$Y_13ahead[1:size.train]
	Y_14ahead_all = list.all.places[[places_to_fit[1]]]$Y_14ahead[1:size.train]
	Y_15ahead_all = list.all.places[[places_to_fit[1]]]$Y_15ahead[1:size.train]
	Y_16ahead_all = list.all.places[[places_to_fit[1]]]$Y_16ahead[1:size.train]


	Y_17ahead_all = list.all.places[[places_to_fit[1]]]$Y_17ahead[1:size.train]
	Y_18ahead_all = list.all.places[[places_to_fit[1]]]$Y_18ahead[1:size.train]
	Y_19ahead_all = list.all.places[[places_to_fit[1]]]$Y_19ahead[1:size.train]
	Y_20ahead_all = list.all.places[[places_to_fit[1]]]$Y_20ahead[1:size.train]
	Y_21ahead_all = list.all.places[[places_to_fit[1]]]$Y_21ahead[1:size.train]
	Y_22ahead_all = list.all.places[[places_to_fit[1]]]$Y_22ahead[1:size.train]
	Y_23ahead_all = list.all.places[[places_to_fit[1]]]$Y_23ahead[1:size.train]
	Y_24ahead_all = list.all.places[[places_to_fit[1]]]$Y_24ahead[1:size.train]
	Y_25ahead_all = list.all.places[[places_to_fit[1]]]$Y_25ahead[1:size.train]
	Y_26ahead_all = list.all.places[[places_to_fit[1]]]$Y_26ahead[1:size.train]

	Y_27ahead_all = list.all.places[[places_to_fit[1]]]$Y_27ahead[1:size.train]
	Y_28ahead_all = list.all.places[[places_to_fit[1]]]$Y_28ahead[1:size.train]
	Y_29ahead_all = list.all.places[[places_to_fit[1]]]$Y_29ahead[1:size.train]
	Y_30ahead_all = list.all.places[[places_to_fit[1]]]$Y_30ahead[1:size.train]
	Y_31ahead_all = list.all.places[[places_to_fit[1]]]$Y_31ahead[1:size.train]
	Y_32ahead_all = list.all.places[[places_to_fit[1]]]$Y_32ahead[1:size.train]
	Y_33ahead_all = list.all.places[[places_to_fit[1]]]$Y_33ahead[1:size.train]
	Y_34ahead_all = list.all.places[[places_to_fit[1]]]$Y_34ahead[1:size.train]
	Y_35ahead_all = list.all.places[[places_to_fit[1]]]$Y_35ahead[1:size.train]
	Y_36ahead_all = list.all.places[[places_to_fit[1]]]$Y_36ahead[1:size.train]
	Y_37ahead_all = list.all.places[[places_to_fit[1]]]$Y_37ahead[1:size.train]
	Y_38ahead_all = list.all.places[[places_to_fit[1]]]$Y_38ahead[1:size.train]
	Y_39ahead_all = list.all.places[[places_to_fit[1]]]$Y_39ahead[1:size.train]
	Y_40ahead_all = list.all.places[[places_to_fit[1]]]$Y_40ahead[1:size.train]
	Y_41ahead_all = list.all.places[[places_to_fit[1]]]$Y_41ahead[1:size.train]
	Y_42ahead_all = list.all.places[[places_to_fit[1]]]$Y_42ahead[1:size.train]
	Y_43ahead_all = list.all.places[[places_to_fit[1]]]$Y_43ahead[1:size.train]
	Y_44ahead_all = list.all.places[[places_to_fit[1]]]$Y_44ahead[1:size.train]
	Y_45ahead_all = list.all.places[[places_to_fit[1]]]$Y_45ahead[1:size.train]
	Y_46ahead_all = list.all.places[[places_to_fit[1]]]$Y_46ahead[1:size.train]
	Y_47ahead_all = list.all.places[[places_to_fit[1]]]$Y_47ahead[1:size.train]
	Y_48ahead_all = list.all.places[[places_to_fit[1]]]$Y_48ahead[1:size.train]
	Y_49ahead_all = list.all.places[[places_to_fit[1]]]$Y_49ahead[1:size.train]
	Y_50ahead_all = list.all.places[[places_to_fit[1]]]$Y_50ahead[1:size.train]
	Y_51ahead_all = list.all.places[[places_to_fit[1]]]$Y_51ahead[1:size.train]
	Y_52ahead_all = list.all.places[[places_to_fit[1]]]$Y_52ahead[1:size.train]

	for (test_place in places_to_fit[-1]){
		X_all = rbind(X_all,list.all.places[[test_place]]$X[1:size.train,])

		Y_oneahead_all = c(Y_oneahead_all, list.all.places[[test_place]]$Y_oneahead[1:size.train])
		Y_twoahead_all = c(Y_twoahead_all, list.all.places[[test_place]]$Y_twoahead[1:size.train])
		Y_threeahead_all = c(Y_threeahead_all, list.all.places[[test_place]]$Y_threeahead[1:size.train])
		Y_fourahead_all = c(Y_fourahead_all, list.all.places[[test_place]]$Y_fourahead[1:size.train])
		Y_fiveahead_all = c(Y_fiveahead_all, list.all.places[[test_place]]$Y_fiveahead[1:size.train])
		Y_sixahead_all = c(Y_sixahead_all, list.all.places[[test_place]]$Y_sixahead[1:size.train])
		Y_sevenahead_all = c(Y_sevenahead_all, list.all.places[[test_place]]$Y_sevenahead[1:size.train])
		Y_eightahead_all = c(Y_eightahead_all, list.all.places[[test_place]]$Y_eightahead[1:size.train])

		Y_9ahead_all = c(Y_9ahead_all, list.all.places[[test_place]]$Y_9ahead[1:size.train])
		Y_10ahead_all = c(Y_10ahead_all, list.all.places[[test_place]]$Y_10ahead[1:size.train])
		Y_11ahead_all = c(Y_11ahead_all, list.all.places[[test_place]]$Y_11ahead[1:size.train])
		Y_12ahead_all = c(Y_12ahead_all, list.all.places[[test_place]]$Y_12ahead[1:size.train])
		Y_13ahead_all = c(Y_13ahead_all, list.all.places[[test_place]]$Y_13ahead[1:size.train])
		Y_14ahead_all = c(Y_14ahead_all, list.all.places[[test_place]]$Y_14ahead[1:size.train])
		Y_15ahead_all = c(Y_15ahead_all, list.all.places[[test_place]]$Y_15ahead[1:size.train])
		Y_16ahead_all = c(Y_16ahead_all, list.all.places[[test_place]]$Y_16ahead[1:size.train])


		Y_17ahead_all = c(Y_17ahead_all, list.all.places[[test_place]]$Y_17ahead[1:size.train])
		Y_18ahead_all = c(Y_18ahead_all, list.all.places[[test_place]]$Y_18ahead[1:size.train])
		Y_19ahead_all = c(Y_19ahead_all, list.all.places[[test_place]]$Y_19ahead[1:size.train])
		Y_20ahead_all = c(Y_20ahead_all, list.all.places[[test_place]]$Y_20ahead[1:size.train])
		Y_21ahead_all = c(Y_21ahead_all, list.all.places[[test_place]]$Y_21ahead[1:size.train])
		Y_22ahead_all = c(Y_22ahead_all, list.all.places[[test_place]]$Y_22ahead[1:size.train])
		Y_23ahead_all = c(Y_23ahead_all, list.all.places[[test_place]]$Y_23ahead[1:size.train])
		Y_24ahead_all = c(Y_24ahead_all, list.all.places[[test_place]]$Y_24ahead[1:size.train])
		Y_25ahead_all = c(Y_25ahead_all, list.all.places[[test_place]]$Y_25ahead[1:size.train])
		Y_26ahead_all = c(Y_26ahead_all, list.all.places[[test_place]]$Y_26ahead[1:size.train])

		Y_27ahead_all = c(Y_27ahead_all, list.all.places[[test_place]]$Y_27ahead[1:size.train])
		Y_28ahead_all = c(Y_28ahead_all, list.all.places[[test_place]]$Y_28ahead[1:size.train])
		Y_29ahead_all = c(Y_29ahead_all, list.all.places[[test_place]]$Y_29ahead[1:size.train])
		Y_30ahead_all = c(Y_30ahead_all, list.all.places[[test_place]]$Y_30ahead[1:size.train])
		Y_31ahead_all = c(Y_31ahead_all, list.all.places[[test_place]]$Y_31ahead[1:size.train])
		Y_32ahead_all = c(Y_32ahead_all, list.all.places[[test_place]]$Y_32ahead[1:size.train])
		Y_33ahead_all = c(Y_33ahead_all, list.all.places[[test_place]]$Y_33ahead[1:size.train])
		Y_34ahead_all = c(Y_34ahead_all, list.all.places[[test_place]]$Y_34ahead[1:size.train])
		Y_35ahead_all = c(Y_35ahead_all, list.all.places[[test_place]]$Y_35ahead[1:size.train])
		Y_36ahead_all = c(Y_36ahead_all, list.all.places[[test_place]]$Y_36ahead[1:size.train])
		Y_37ahead_all = c(Y_37ahead_all, list.all.places[[test_place]]$Y_37ahead[1:size.train])
		Y_38ahead_all = c(Y_38ahead_all, list.all.places[[test_place]]$Y_38ahead[1:size.train])
		Y_39ahead_all = c(Y_39ahead_all, list.all.places[[test_place]]$Y_39ahead[1:size.train])
		Y_40ahead_all = c(Y_40ahead_all, list.all.places[[test_place]]$Y_40ahead[1:size.train])
		Y_41ahead_all = c(Y_41ahead_all, list.all.places[[test_place]]$Y_41ahead[1:size.train])
		Y_42ahead_all = c(Y_42ahead_all, list.all.places[[test_place]]$Y_42ahead[1:size.train])
		Y_43ahead_all = c(Y_43ahead_all, list.all.places[[test_place]]$Y_43ahead[1:size.train])
		Y_44ahead_all = c(Y_44ahead_all, list.all.places[[test_place]]$Y_44ahead[1:size.train])
		Y_45ahead_all = c(Y_45ahead_all, list.all.places[[test_place]]$Y_45ahead[1:size.train])
		Y_46ahead_all = c(Y_46ahead_all, list.all.places[[test_place]]$Y_46ahead[1:size.train])
		Y_47ahead_all = c(Y_47ahead_all, list.all.places[[test_place]]$Y_47ahead[1:size.train])
		Y_48ahead_all = c(Y_48ahead_all, list.all.places[[test_place]]$Y_48ahead[1:size.train])
		Y_49ahead_all = c(Y_49ahead_all, list.all.places[[test_place]]$Y_49ahead[1:size.train])
		Y_50ahead_all = c(Y_50ahead_all, list.all.places[[test_place]]$Y_50ahead[1:size.train])
		Y_51ahead_all = c(Y_51ahead_all, list.all.places[[test_place]]$Y_51ahead[1:size.train])
		Y_52ahead_all = c(Y_52ahead_all, list.all.places[[test_place]]$Y_52ahead[1:size.train])

	}



	## Fit Lasso ##
	###############

	## one-step prediction ##
	#########################

	lower_limits = -Inf
	intercept= TRUE

	cv.lasso.oneahead <- cv.glmnet(X_all, Y_oneahead_all, alpha = 1, lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.oneahead)

	model.oneahead <- glmnet(X_all, Y_oneahead_all, alpha = 1, lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.oneahead$lambda.1se, intercept=intercept)

	coef.oneahead = coef(model.oneahead)[,1]
	non.zero.coef.oneahead = coef.oneahead[which(coef.oneahead!=0)]


	## two-step prediction ##
	#########################

	cv.lasso.twoahead <- cv.glmnet(X_all, Y_twoahead_all, alpha = 1, lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.twoahead)

	model.twoahead <- glmnet(X_all, Y_twoahead_all, alpha = 1, lower.limits = lower_limits, family = "poisson",
	                lambda = cv.lasso.twoahead$lambda.1se, intercept=intercept)

	coef.twoahead = coef(model.twoahead)[,1]
	non.zero.coef.twoahead = coef.twoahead[which(coef.twoahead!=0)]



	## three-step prediction ##
	#########################

	cv.lasso.threeahead <- cv.glmnet(X_all, Y_threeahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",, intercept=intercept)
	# plot(cv.lasso.threeahead)

	model.threeahead <- glmnet(X_all, Y_threeahead_all, alpha = 1, lower.limits = lower_limits, family = "poisson",
	                lambda = cv.lasso.threeahead$lambda.1se, intercept=intercept)

	coef.threeahead = coef(model.threeahead)[,1]
	non.zero.coef.threeahead = coef.threeahead[which(coef.threeahead!=0)]


	## four-step prediction ##
	#########################

	cv.lasso.fourahead <- cv.glmnet(X_all, Y_fourahead_all, alpha = 1, lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.fourahead)

	model.fourahead <- glmnet(X_all, Y_fourahead_all, alpha = 1, lower.limits = lower_limits, family = "poisson",
	                lambda = cv.lasso.fourahead$lambda.1se, intercept=intercept)

	coef.fourahead = coef(model.fourahead)[,1]
	non.zero.coef.fourahead = coef.fourahead[which(coef.fourahead!=0)]





	## five-step prediction ##
	#########################

	cv.lasso.fiveahead <- cv.glmnet(X_all, Y_fiveahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.fiveahead)

	model.fiveahead <- glmnet(X_all, Y_fiveahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.fiveahead$lambda.1se, intercept=intercept)

	coef.fiveahead = coef(model.fiveahead)[,1]
	non.zero.coef.fiveahead = coef.fiveahead[which(coef.fiveahead!=0)]




	## six-step prediction ##
	#########################
	cv.lasso.sixahead <- cv.glmnet(X_all, Y_sixahead_all, alpha = 1, lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.sixahead)

	model.sixahead <- glmnet(X_all, Y_sixahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.sixahead$lambda.1se, intercept=intercept)

	coef.sixahead = coef(model.sixahead)[,1]
	non.zero.coef.sixahead = coef.sixahead[which(coef.sixahead!=0)]



	## seven-step prediction ##
	#########################

	cv.lasso.sevenahead <- cv.glmnet(X_all, Y_sevenahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson", intercept=intercept)
	# plot(cv.lasso.sevenahead)

	model.sevenahead <- glmnet(X_all, Y_sevenahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.sevenahead$lambda.1se, intercept=intercept)

	coef.sevenahead = coef(model.sevenahead)[,1]
	non.zero.coef.sevenahead = coef.sevenahead[which(coef.sevenahead!=0)]


	## eight-step prediction ##
	#########################

	cv.lasso.eightahead <- cv.glmnet(X_all, Y_eightahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.eightahead)

	model.eightahead <- glmnet(X_all, Y_eightahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.eightahead$lambda.1se, intercept=intercept)

	coef.eightahead = coef(model.eightahead)[,1]
	non.zero.coef.eightahead = coef.eightahead[which(coef.eightahead!=0)]



	## 9-step prediction ##
	#########################

	cv.lasso.9ahead <- cv.glmnet(X_all, Y_9ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.9ahead)

	model.9ahead <- glmnet(X_all, Y_9ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.9ahead$lambda.1se, intercept=intercept)

	coef.9ahead = coef(model.9ahead)[,1]
	non.zero.coef.9ahead = coef.9ahead[which(coef.9ahead!=0)]


	## 10-step prediction ##
	#########################

	cv.lasso.10ahead <- cv.glmnet(X_all, Y_10ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.10ahead)

	model.10ahead <- glmnet(X_all, Y_10ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.10ahead$lambda.1se, intercept=intercept)

	coef.10ahead = coef(model.10ahead)[,1]
	non.zero.coef.10ahead = coef.10ahead[which(coef.10ahead!=0)]



	## 11-step prediction ##
	#########################

	cv.lasso.11ahead <- cv.glmnet(X_all, Y_11ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.11ahead)

	model.11ahead <- glmnet(X_all, Y_11ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.11ahead$lambda.1se, intercept=intercept)

	coef.11ahead = coef(model.11ahead)[,1]
	non.zero.coef.11ahead = coef.11ahead[which(coef.11ahead!=0)]


	## 12-step prediction ##
	#########################

	cv.lasso.12ahead <- cv.glmnet(X_all, Y_12ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.12ahead)

	model.12ahead <- glmnet(X_all, Y_12ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.12ahead$lambda.1se, intercept=intercept)

	coef.12ahead = coef(model.12ahead)[,1]
	non.zero.coef.12ahead = coef.12ahead[which(coef.12ahead!=0)]


	## 13-step prediction ##
	#########################

	cv.lasso.13ahead <- cv.glmnet(X_all, Y_13ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.13ahead)

	model.13ahead <- glmnet(X_all, Y_13ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.13ahead$lambda.1se, intercept=intercept)

	coef.13ahead = coef(model.13ahead)[,1]
	non.zero.coef.13ahead = coef.13ahead[which(coef.13ahead!=0)]



	## 14-step prediction ##
	#########################

	cv.lasso.14ahead <- cv.glmnet(X_all, Y_14ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.25ahead)

	model.14ahead <- glmnet(X_all, Y_14ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.14ahead$lambda.1se, intercept=intercept)

	coef.14ahead = coef(model.14ahead)[,1]
	non.zero.coef.14ahead = coef.14ahead[which(coef.14ahead!=0)]


	## 15-step prediction ##
	#########################

	cv.lasso.15ahead <- cv.glmnet(X_all, Y_15ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.15ahead)

	model.15ahead <- glmnet(X_all, Y_15ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.15ahead$lambda.1se, intercept=intercept)

	coef.15ahead = coef(model.15ahead)[,1]
	non.zero.coef.15ahead = coef.15ahead[which(coef.15ahead!=0)]


	## 16-step prediction ##
	#########################

	cv.lasso.16ahead <- cv.glmnet(X_all, Y_16ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.16ahead)

	model.16ahead <- glmnet(X_all, Y_16ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.16ahead$lambda.1se, intercept=intercept)

	coef.16ahead = coef(model.16ahead)[,1]
	non.zero.coef.16ahead = coef.16ahead[which(coef.16ahead!=0)]


	## 17-step prediction ##
	#########################

	cv.lasso.17ahead <- cv.glmnet(X_all, Y_17ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.17ahead)

	model.17ahead <- glmnet(X_all, Y_17ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.17ahead$lambda.1se, intercept=intercept)

	coef.17ahead = coef(model.17ahead)[,1]
	non.zero.coef.17ahead = coef.17ahead[which(coef.17ahead!=0)]


	## 18-step prediction ##
	#########################

	cv.lasso.18ahead <- cv.glmnet(X_all, Y_18ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.18ahead)

	model.18ahead <- glmnet(X_all, Y_18ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.18ahead$lambda.1se, intercept=intercept)

	coef.18ahead = coef(model.18ahead)[,1]
	non.zero.coef.18ahead = coef.18ahead[which(coef.18ahead!=0)]



	## 19-step prediction ##
	#########################

	cv.lasso.19ahead <- cv.glmnet(X_all, Y_19ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.19ahead)

	model.19ahead <- glmnet(X_all, Y_19ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.19ahead$lambda.1se, intercept=intercept)

	coef.19ahead = coef(model.19ahead)[,1]
	non.zero.coef.19ahead = coef.19ahead[which(coef.19ahead!=0)]



	## 20-step prediction ##
	#########################

	cv.lasso.20ahead <- cv.glmnet(X_all, Y_20ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.20ahead)

	model.20ahead <- glmnet(X_all, Y_20ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.20ahead$lambda.1se, intercept=intercept)

	coef.20ahead = coef(model.20ahead)[,1]
	non.zero.coef.20ahead = coef.20ahead[which(coef.20ahead!=0)]


	## 21-step prediction ##
	#########################

	cv.lasso.21ahead <- cv.glmnet(X_all, Y_21ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.21ahead)

	model.21ahead <- glmnet(X_all, Y_21ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.21ahead$lambda.1se, intercept=intercept)

	coef.21ahead = coef(model.21ahead)[,1]
	non.zero.coef.21ahead = coef.21ahead[which(coef.21ahead!=0)]


	## 22-step prediction ##
	#########################

	cv.lasso.22ahead <- cv.glmnet(X_all, Y_22ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.22ahead)

	model.22ahead <- glmnet(X_all, Y_22ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.22ahead$lambda.1se, intercept=intercept)

	coef.22ahead = coef(model.22ahead)[,1]
	non.zero.coef.22ahead = coef.22ahead[which(coef.22ahead!=0)]


	## 23-step prediction ##
	#########################

	cv.lasso.23ahead <- cv.glmnet(X_all, Y_23ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.23ahead)

	model.23ahead <- glmnet(X_all, Y_23ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.23ahead$lambda.1se, intercept=intercept)

	coef.23ahead = coef(model.23ahead)[,1]
	non.zero.coef.23ahead = coef.23ahead[which(coef.23ahead!=0)]



	## 24-step prediction ##
	#########################

	cv.lasso.24ahead <- cv.glmnet(X_all, Y_24ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.24ahead)

	model.24ahead <- glmnet(X_all, Y_24ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.24ahead$lambda.1se, intercept=intercept)

	coef.24ahead = coef(model.24ahead)[,1]
	non.zero.coef.24ahead = coef.24ahead[which(coef.24ahead!=0)]


	## 25-step prediction ##
	#########################

	cv.lasso.25ahead <- cv.glmnet(X_all, Y_25ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.25ahead)

	model.25ahead <- glmnet(X_all, Y_25ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.25ahead$lambda.1se, intercept=intercept)

	coef.25ahead = coef(model.25ahead)[,1]
	non.zero.coef.25ahead = coef.25ahead[which(coef.25ahead!=0)]



	## 26-step prediction ##
	#########################

	cv.lasso.26ahead <- cv.glmnet(X_all, Y_26ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.26ahead)

	model.26ahead <- glmnet(X_all, Y_26ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.26ahead$lambda.1se, intercept=intercept)

	coef.26ahead = coef(model.26ahead)[,1]
	non.zero.coef.26ahead = coef.26ahead[which(coef.26ahead!=0)]

	## 27-step prediction ##
	#########################

	cv.lasso.27ahead <- cv.glmnet(X_all, Y_27ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.27ahead)

	model.27ahead <- glmnet(X_all, Y_27ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.27ahead$lambda.1se, intercept=intercept)

	coef.27ahead = coef(model.27ahead)[,1]
	non.zero.coef.27ahead = coef.27ahead[which(coef.27ahead!=0)]


	## 28-step prediction ##
	#########################

	cv.lasso.28ahead <- cv.glmnet(X_all, Y_28ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.28ahead)

	model.28ahead <- glmnet(X_all, Y_28ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.28ahead$lambda.1se, intercept=intercept)

	coef.28ahead = coef(model.28ahead)[,1]
	non.zero.coef.28ahead = coef.28ahead[which(coef.28ahead!=0)]


	## 29-step prediction ##
	#########################

	cv.lasso.29ahead <- cv.glmnet(X_all, Y_29ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.29ahead)

	model.29ahead <- glmnet(X_all, Y_29ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.29ahead$lambda.1se, intercept=intercept)

	coef.29ahead = coef(model.29ahead)[,1]
	non.zero.coef.29ahead = coef.29ahead[which(coef.29ahead!=0)]


	## 30-step prediction ##
	#########################

	cv.lasso.30ahead <- cv.glmnet(X_all, Y_30ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.30ahead)

	model.30ahead <- glmnet(X_all, Y_30ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.30ahead$lambda.1se, intercept=intercept)

	coef.30ahead = coef(model.30ahead)[,1]
	non.zero.coef.30ahead = coef.30ahead[which(coef.30ahead!=0)]

	## 31-step prediction ##
	#########################

	cv.lasso.31ahead <- cv.glmnet(X_all, Y_31ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.31ahead)

	model.31ahead <- glmnet(X_all, Y_31ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.31ahead$lambda.1se, intercept=intercept)

	coef.31ahead = coef(model.31ahead)[,1]
	non.zero.coef.31ahead = coef.31ahead[which(coef.31ahead!=0)]

	## 32-step prediction ##
	#########################

	cv.lasso.32ahead <- cv.glmnet(X_all, Y_32ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.32ahead)

	model.32ahead <- glmnet(X_all, Y_32ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.32ahead$lambda.1se, intercept=intercept)

	coef.32ahead = coef(model.32ahead)[,1]
	non.zero.coef.32ahead = coef.32ahead[which(coef.32ahead!=0)]

	## 33-step prediction ##
	#########################

	cv.lasso.33ahead <- cv.glmnet(X_all, Y_33ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.33ahead)

	model.33ahead <- glmnet(X_all, Y_33ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.33ahead$lambda.1se, intercept=intercept)

	coef.33ahead = coef(model.33ahead)[,1]
	non.zero.coef.33ahead = coef.33ahead[which(coef.33ahead!=0)]

	## 34-step prediction ##
	#########################

	cv.lasso.34ahead <- cv.glmnet(X_all, Y_34ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.34ahead)

	model.34ahead <- glmnet(X_all, Y_34ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.34ahead$lambda.1se, intercept=intercept)

	coef.34ahead = coef(model.34ahead)[,1]
	non.zero.coef.34ahead = coef.34ahead[which(coef.34ahead!=0)]

	## 35-step prediction ##
	#########################

	cv.lasso.35ahead <- cv.glmnet(X_all, Y_35ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.35ahead)

	model.35ahead <- glmnet(X_all, Y_35ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.35ahead$lambda.1se, intercept=intercept)

	coef.35ahead = coef(model.35ahead)[,1]
	non.zero.coef.35ahead = coef.35ahead[which(coef.35ahead!=0)]


	## 36-step prediction ##
	#########################

	cv.lasso.36ahead <- cv.glmnet(X_all, Y_36ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.36ahead)

	model.36ahead <- glmnet(X_all, Y_36ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.36ahead$lambda.1se, intercept=intercept)

	coef.36ahead = coef(model.36ahead)[,1]
	non.zero.coef.36ahead = coef.36ahead[which(coef.36ahead!=0)]


	## 37-step prediction ##
	#########################

	cv.lasso.37ahead <- cv.glmnet(X_all, Y_37ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.37ahead)

	model.37ahead <- glmnet(X_all, Y_37ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.37ahead$lambda.1se, intercept=intercept)

	coef.37ahead = coef(model.37ahead)[,1]
	non.zero.coef.37ahead = coef.37ahead[which(coef.37ahead!=0)]

	## 38-step prediction ##
	#########################

	cv.lasso.38ahead <- cv.glmnet(X_all, Y_38ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.38ahead)

	model.38ahead <- glmnet(X_all, Y_38ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.38ahead$lambda.1se, intercept=intercept)

	coef.38ahead = coef(model.38ahead)[,1]
	non.zero.coef.38ahead = coef.38ahead[which(coef.38ahead!=0)]


	## 39-step prediction ##
	#########################

	cv.lasso.39ahead <- cv.glmnet(X_all, Y_39ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.39ahead)

	model.39ahead <- glmnet(X_all, Y_39ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.39ahead$lambda.1se, intercept=intercept)

	coef.39ahead = coef(model.39ahead)[,1]
	non.zero.coef.39ahead = coef.39ahead[which(coef.39ahead!=0)]



	## 40-step prediction ##
	#########################

	cv.lasso.40ahead <- cv.glmnet(X_all, Y_40ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.40ahead)

	model.40ahead <- glmnet(X_all, Y_40ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.40ahead$lambda.1se, intercept=intercept)

	coef.40ahead = coef(model.40ahead)[,1]
	non.zero.coef.40ahead = coef.40ahead[which(coef.40ahead!=0)]


	## 41-step prediction ##
	#########################

	cv.lasso.41ahead <- cv.glmnet(X_all, Y_41ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.41ahead)

	model.41ahead <- glmnet(X_all, Y_41ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.41ahead$lambda.1se, intercept=intercept)

	coef.41ahead = coef(model.41ahead)[,1]
	non.zero.coef.41ahead = coef.41ahead[which(coef.41ahead!=0)]


	## 26-step prediction ##
	#########################

	cv.lasso.42ahead <- cv.glmnet(X_all, Y_42ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.42ahead)

	model.42ahead <- glmnet(X_all, Y_42ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.42ahead$lambda.1se, intercept=intercept)

	coef.42ahead = coef(model.42ahead)[,1]
	non.zero.coef.42ahead = coef.42ahead[which(coef.42ahead!=0)]


	## 26-step prediction ##
	#########################

	cv.lasso.43ahead <- cv.glmnet(X_all, Y_43ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.43ahead)

	model.43ahead <- glmnet(X_all, Y_43ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.43ahead$lambda.1se, intercept=intercept)

	coef.43ahead = coef(model.43ahead)[,1]
	non.zero.coef.43ahead = coef.43ahead[which(coef.43ahead!=0)]


	## 44-step prediction ##
	#########################

	cv.lasso.44ahead <- cv.glmnet(X_all, Y_44ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.44ahead)

	model.44ahead <- glmnet(X_all, Y_44ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.44ahead$lambda.1se, intercept=intercept)

	coef.44ahead = coef(model.44ahead)[,1]
	non.zero.coef.44ahead = coef.44ahead[which(coef.26ahead!=0)]


	## 45-step prediction ##
	#########################

	cv.lasso.45ahead <- cv.glmnet(X_all, Y_45ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.45ahead)

	model.45ahead <- glmnet(X_all, Y_45ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.45ahead$lambda.1se, intercept=intercept)

	coef.45ahead = coef(model.45ahead)[,1]
	non.zero.coef.45ahead = coef.45ahead[which(coef.45ahead!=0)]


	## 46-step prediction ##
	#########################

	cv.lasso.46ahead <- cv.glmnet(X_all, Y_46ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.46ahead)

	model.46ahead <- glmnet(X_all, Y_46ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.46ahead$lambda.1se, intercept=intercept)

	coef.46ahead = coef(model.46ahead)[,1]
	non.zero.coef.46ahead = coef.46ahead[which(coef.46ahead!=0)]


	## 26-step prediction ##
	#########################

	cv.lasso.47ahead <- cv.glmnet(X_all, Y_47ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.47ahead)

	model.47ahead <- glmnet(X_all, Y_47ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.47ahead$lambda.1se, intercept=intercept)

	coef.47ahead = coef(model.47ahead)[,1]
	non.zero.coef.47ahead = coef.47ahead[which(coef.47ahead!=0)]


	## 48-step prediction ##
	#########################

	cv.lasso.48ahead <- cv.glmnet(X_all, Y_48ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.48ahead)

	model.48ahead <- glmnet(X_all, Y_48ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.48ahead$lambda.1se, intercept=intercept)

	coef.48ahead = coef(model.48ahead)[,1]
	non.zero.coef.48ahead = coef.48ahead[which(coef.48ahead!=0)]

	## 49-step prediction ##
	#########################

	cv.lasso.49ahead <- cv.glmnet(X_all, Y_49ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.49ahead)

	model.49ahead <- glmnet(X_all, Y_49ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.49ahead$lambda.1se, intercept=intercept)

	coef.49ahead = coef(model.49ahead)[,1]
	non.zero.coef.49ahead = coef.49ahead[which(coef.49ahead!=0)]

	## 50-step prediction ##
	#########################

	cv.lasso.50ahead <- cv.glmnet(X_all, Y_50ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.50ahead)

	model.50ahead <- glmnet(X_all, Y_50ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.50ahead$lambda.1se, intercept=intercept)

	coef.50ahead = coef(model.50ahead)[,1]
	non.zero.coef.50ahead = coef.50ahead[which(coef.50ahead!=0)]


	## 51-step prediction ##
	#########################

	cv.lasso.51ahead <- cv.glmnet(X_all, Y_51ahead_all, alpha = 1,lower.limits = lower_limits, family = "poisson", intercept=intercept)
	# plot(cv.lasso.51ahead)

	model.51ahead <- glmnet(X_all, Y_51ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.51ahead$lambda.1se, intercept=intercept)

	coef.51ahead = coef(model.51ahead)[,1]
	non.zero.coef.51ahead = coef.51ahead[which(coef.51ahead!=0)]

	## 52-step prediction ##
	#########################

	cv.lasso.52ahead <-  cv.glmnet(X_all, Y_52ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",, intercept=intercept)

	model.52ahead <- glmnet(X_all, Y_52ahead_all, alpha = 1,lower.limits = lower_limits,  family = "poisson",
	                lambda = cv.lasso.52ahead$lambda.1se, intercept=intercept)

	coef.52ahead = coef(model.52ahead)[,1]
	non.zero.coef.52ahead = coef.52ahead[which(coef.52ahead!=0)]

	for (test_place in places_to_fit){

		y.hat.oneahead = predict(model.oneahead,newx=list.all.places[[test_place]]$X)
		y.hat.twoahead = predict(model.twoahead,newx=list.all.places[[test_place]]$X)
		y.hat.threeahead = predict(model.threeahead,newx=list.all.places[[test_place]]$X)
		y.hat.fourahead = predict(model.fourahead,newx=list.all.places[[test_place]]$X)
		y.hat.fiveahead = predict(model.fiveahead,newx=list.all.places[[test_place]]$X)
		y.hat.sixahead = predict(model.sixahead,newx=list.all.places[[test_place]]$X)
		y.hat.sevenahead = predict(model.sevenahead,newx=list.all.places[[test_place]]$X)
		y.hat.eightahead = predict(model.eightahead,newx=list.all.places[[test_place]]$X)

		y.hat.9ahead = predict(model.9ahead,newx=list.all.places[[test_place]]$X)
		y.hat.10ahead = predict(model.10ahead,newx=list.all.places[[test_place]]$X)
		y.hat.11ahead = predict(model.11ahead,newx=list.all.places[[test_place]]$X)
		y.hat.12ahead = predict(model.12ahead,newx=list.all.places[[test_place]]$X)
		y.hat.13ahead = predict(model.13ahead,newx=list.all.places[[test_place]]$X)
		y.hat.14ahead = predict(model.14ahead,newx=list.all.places[[test_place]]$X)
		y.hat.15ahead = predict(model.15ahead,newx=list.all.places[[test_place]]$X)
		y.hat.16ahead = predict(model.16ahead,newx=list.all.places[[test_place]]$X)

		y.hat.17ahead = predict(model.17ahead,newx=list.all.places[[test_place]]$X)
		y.hat.18ahead = predict(model.18ahead,newx=list.all.places[[test_place]]$X)
		y.hat.19ahead = predict(model.19ahead,newx=list.all.places[[test_place]]$X)
		y.hat.20ahead = predict(model.20ahead,newx=list.all.places[[test_place]]$X)
		y.hat.21ahead = predict(model.21ahead,newx=list.all.places[[test_place]]$X)
		y.hat.22ahead = predict(model.22ahead,newx=list.all.places[[test_place]]$X)
		y.hat.23ahead = predict(model.23ahead,newx=list.all.places[[test_place]]$X)
		y.hat.24ahead = predict(model.24ahead,newx=list.all.places[[test_place]]$X)
		y.hat.25ahead = predict(model.25ahead,newx=list.all.places[[test_place]]$X)
		y.hat.26ahead = predict(model.26ahead,newx=list.all.places[[test_place]]$X)


		y.hat.27ahead = predict(model.27ahead,newx=list.all.places[[test_place]]$X)
		y.hat.28ahead = predict(model.28ahead,newx=list.all.places[[test_place]]$X)
		y.hat.29ahead = predict(model.29ahead,newx=list.all.places[[test_place]]$X)
		y.hat.30ahead = predict(model.30ahead,newx=list.all.places[[test_place]]$X)
		y.hat.31ahead = predict(model.31ahead,newx=list.all.places[[test_place]]$X)
		y.hat.32ahead = predict(model.32ahead,newx=list.all.places[[test_place]]$X)
		y.hat.33ahead = predict(model.33ahead,newx=list.all.places[[test_place]]$X)
		y.hat.34ahead = predict(model.34ahead,newx=list.all.places[[test_place]]$X)
		y.hat.35ahead = predict(model.35ahead,newx=list.all.places[[test_place]]$X)
		y.hat.36ahead = predict(model.36ahead,newx=list.all.places[[test_place]]$X)
		y.hat.37ahead = predict(model.37ahead,newx=list.all.places[[test_place]]$X)
		y.hat.38ahead = predict(model.38ahead,newx=list.all.places[[test_place]]$X)
		y.hat.39ahead = predict(model.39ahead,newx=list.all.places[[test_place]]$X)
		y.hat.40ahead = predict(model.40ahead,newx=list.all.places[[test_place]]$X)
		y.hat.41ahead = predict(model.41ahead,newx=list.all.places[[test_place]]$X)
		y.hat.42ahead = predict(model.42ahead,newx=list.all.places[[test_place]]$X)
		y.hat.43ahead = predict(model.43ahead,newx=list.all.places[[test_place]]$X)
		y.hat.44ahead = predict(model.44ahead,newx=list.all.places[[test_place]]$X)
		y.hat.45ahead = predict(model.45ahead,newx=list.all.places[[test_place]]$X)
		y.hat.46ahead = predict(model.46ahead,newx=list.all.places[[test_place]]$X)
		y.hat.47ahead = predict(model.47ahead,newx=list.all.places[[test_place]]$X)
		y.hat.48ahead = predict(model.48ahead,newx=list.all.places[[test_place]]$X)
		y.hat.49ahead = predict(model.49ahead,newx=list.all.places[[test_place]]$X)
		y.hat.50ahead = predict(model.50ahead,newx=list.all.places[[test_place]]$X)
		y.hat.51ahead = predict(model.51ahead,newx=list.all.places[[test_place]]$X)
		y.hat.52ahead = predict(model.52ahead,newx=list.all.places[[test_place]]$X)
	
	

	
		list.all.places[[test_place]]$y.hat.oneahead = y.hat.oneahead
		list.all.places[[test_place]]$y.hat.twoahead = y.hat.twoahead
		list.all.places[[test_place]]$y.hat.threeahead = y.hat.threeahead
		list.all.places[[test_place]]$y.hat.fourahead = y.hat.fourahead
		list.all.places[[test_place]]$y.hat.fiveahead = y.hat.fiveahead
		list.all.places[[test_place]]$y.hat.sixahead = y.hat.sixahead
		list.all.places[[test_place]]$y.hat.sevenahead = y.hat.sevenahead
		list.all.places[[test_place]]$y.hat.eightahead = y.hat.eightahead


		list.all.places[[test_place]]$y.hat.9ahead = y.hat.9ahead
		list.all.places[[test_place]]$y.hat.10ahead = y.hat.10ahead
		list.all.places[[test_place]]$y.hat.11ahead = y.hat.11ahead
		list.all.places[[test_place]]$y.hat.12ahead = y.hat.12ahead
		list.all.places[[test_place]]$y.hat.13ahead = y.hat.13ahead
		list.all.places[[test_place]]$y.hat.14ahead = y.hat.14ahead
		list.all.places[[test_place]]$y.hat.15ahead = y.hat.15ahead
		list.all.places[[test_place]]$y.hat.16ahead = y.hat.16ahead

		list.all.places[[test_place]]$y.hat.17ahead = y.hat.17ahead
		list.all.places[[test_place]]$y.hat.18ahead = y.hat.18ahead
		list.all.places[[test_place]]$y.hat.19ahead = y.hat.19ahead
		list.all.places[[test_place]]$y.hat.20ahead = y.hat.20ahead
		list.all.places[[test_place]]$y.hat.21ahead = y.hat.21ahead
		list.all.places[[test_place]]$y.hat.22ahead = y.hat.22ahead
		list.all.places[[test_place]]$y.hat.23ahead = y.hat.23ahead
		list.all.places[[test_place]]$y.hat.24ahead = y.hat.24ahead
		list.all.places[[test_place]]$y.hat.25ahead = y.hat.25ahead
		list.all.places[[test_place]]$y.hat.26ahead = y.hat.26ahead

		list.all.places[[test_place]]$y.hat.27ahead = y.hat.27ahead
		list.all.places[[test_place]]$y.hat.28ahead = y.hat.28ahead
		list.all.places[[test_place]]$y.hat.29ahead = y.hat.29ahead
		list.all.places[[test_place]]$y.hat.30ahead = y.hat.30ahead
		list.all.places[[test_place]]$y.hat.31ahead = y.hat.31ahead
		list.all.places[[test_place]]$y.hat.32ahead = y.hat.32ahead
		list.all.places[[test_place]]$y.hat.33ahead = y.hat.33ahead
		list.all.places[[test_place]]$y.hat.34ahead = y.hat.34ahead
		list.all.places[[test_place]]$y.hat.35ahead = y.hat.35ahead
		list.all.places[[test_place]]$y.hat.36ahead = y.hat.36ahead
		list.all.places[[test_place]]$y.hat.37ahead = y.hat.37ahead
		list.all.places[[test_place]]$y.hat.38ahead = y.hat.38ahead
		list.all.places[[test_place]]$y.hat.39ahead = y.hat.39ahead
		list.all.places[[test_place]]$y.hat.40ahead = y.hat.40ahead
		list.all.places[[test_place]]$y.hat.41ahead = y.hat.41ahead
		list.all.places[[test_place]]$y.hat.42ahead = y.hat.42ahead
		list.all.places[[test_place]]$y.hat.43ahead = y.hat.43ahead
		list.all.places[[test_place]]$y.hat.44ahead = y.hat.44ahead
		list.all.places[[test_place]]$y.hat.45ahead = y.hat.45ahead
		list.all.places[[test_place]]$y.hat.46ahead = y.hat.46ahead
		list.all.places[[test_place]]$y.hat.47ahead = y.hat.47ahead
		list.all.places[[test_place]]$y.hat.48ahead = y.hat.48ahead
		list.all.places[[test_place]]$y.hat.49ahead = y.hat.49ahead
		list.all.places[[test_place]]$y.hat.50ahead = y.hat.50ahead
		list.all.places[[test_place]]$y.hat.51ahead = y.hat.51ahead
		list.all.places[[test_place]]$y.hat.52ahead = y.hat.52ahead


	}
		



} # end of looping model versions

