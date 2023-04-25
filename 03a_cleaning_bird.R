# rm(list=ls())
gc()
options(scipen=999)
# setwd(dir = "~/Desktop/analysis/")

#libraries
library(data.table)
library(dplyr)
sf::sf_use_s2(FALSE)

# Load utility functions
source("utils/utils.R")
source("utils/utils_plot.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # GENERATE JOINT SHAPEFILE  # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# define necessary shapefile - join all the shapefiles from each country on interest
paths_list = c("egypt_shapefiles/EGY_adm2.shp",
               "tun_adm1/TUN_adm1.shp",
               "JOR_adm/JOR_adm1.shp",
               "dza_adm_unhcr2020_shp/dza_admbnda_adm1_unhcr2020.shp",
               "kwt_adm1/KWT_adm1.shp",
               "lbn_adm_cdr_20200810/lbn_admbnda_adm1_cdr_20200810.shp",
               "lby_adm_unosat_lbsc_20180507_SHP/lby_admbnda_adm2_unosat_lbsc_20180507.shp",
               "mar_adm_unhcr2020_shp/mar_admbnda_adm1_unhcr_20201203.shp",
               "yem_adm_govyem_cso_ochayemen_20191002_shp/yem_admbnda_adm1_govyem_cso_20191002.shp",
               "sau_adm1/SAU_adm1.shp",
               "gadm36_ISR_shp/gadm36_ISR_1.shp"
)
adm_colnames_list = c("ADM1_EN",
                      "NAME_1",
                      "NAME_1",
                      "ADM1_EN",
                      "NAME_1",
                      "admin1Name",
                      "ADM2_EN",
                      "ADM1_EN",
                      "ADM1_EN",
                      "NAME_1",
                      "NAME_1"
)
plotnames_list = c("egypt_shape",
                   "tunisia_shape",
                   "jordan_shape",
                   "algeria_shape",
                   "kuwait_shape",
                   "lebanon_shape",
                   "libya_shape",
                   "morocco_shape",
                   "yemen_shape",
                   "saudi arabia_shape",
                   "israel_shape"
)
shape = collect_and_connect_shapes(shape_paths_list = sapply(paths_list,function(x){paste("shapefiles/",x,sep="")}),
                                    adm_colnames_list = adm_colnames_list,
                                    plotnames_list = plotnames_list ,
                                    plot_path = "Plots/Bird/")
# plot full shapefile 
pdf(file = 'Plots/Bird/middle_east.pdf',width =10, height = 10)
plot(shape)
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # CLEAN INDIVIDUAL-LEVEL DATA # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# load individual-level training data 
individual.data_path <- 'data/matchdata_isisall.csv'
survey = fread(file = individual.data_path)

# # # clean Governorate Names in Data
survey $governorate[which(survey $governorate=="Oum el-Bouaghi")] = "Oum El Bouaghi"
survey $governorate[which(survey $governorate=="Algiers")] ="Alger"
survey $governorate[which(survey $governorate=="Médéa")] = "Medea"
survey $governorate[which(survey $governorate=="Bordj Bou Arréridj")] ="Bordj Bou Arrer"
survey $governorate[which(survey $governorate=="El Taref")] ="El-Tarf"
survey $governorate[which(survey $governorate=="Souk Ahras")] ="Souk-Ahras"
survey $governorate[which(survey $governorate=="Tipasa")] ="Tipaza"
survey $governorate[which(survey $governorate=="Aïn Defla")] ="Ain-Defla"
survey $governorate[which(survey $governorate=="Ain Tecmouchent")] ="Ain-Temouchent"
survey $governorate[which(survey $governorate=="Ghardaïa")] ="Ghardaia"
survey $governorate[which(survey $governorate=="Alexander")] = "Alexandria"
survey $governorate[which(survey $governorate=="El Sharqiya")] = "Sharkia"
survey $governorate[which(survey $governorate=="Qalyubia")] ="Kalyoubia"
survey $governorate[which(survey $governorate=="Kafir el-Sheikh")] = "Kafr El-Shikh"
survey $governorate[which(survey $governorate=="The West Bank")] = NA
survey $governorate[which(survey $governorate=="Monufia")] = "Menoufia" 
survey $governorate[which(survey $governorate=="Beheira")] ="Behera"
survey $governorate[which(survey $governorate=="Faiyum")] ="Fayoum"
survey $governorate[which(survey $governorate=="Minya")] ="Menia"
survey $governorate[which(survey $governorate=="Asyut")] ="Assiut"
survey $governorate[which(survey $governorate=="Sohag")] ="Suhag"
survey $governorate[which(survey $governorate=="The Red Sea")] ="Red Sea"
survey $governorate[which(survey $governorate=="The capital")] = "Amman" 
survey $governorate[which(survey $governorate=="Jerash")] ="Jarash"
survey $governorate[which(survey $governorate=="Ajloun")] ="Ajlun"
survey $governorate[which(survey $governorate=="al Karak")] ="Karak"
survey $governorate[which(survey $governorate=="Ma'an")] ="Ma`an"
survey $governorate[which(survey $governorate=="Ahmadi")] ="Al Ahmadi"
survey $governorate[which(survey $governorate=="al-Jahra")] ="Al Jahrah"
survey $governorate[which(survey $governorate=="Mubarak al-Sabah")] ="Mubarak Al-Kabeer"
survey $governorate[which(survey $governorate=="al-Farwaniyah")]  = "Al Farwaniyah"
survey $governorate[which(survey $governorate=="Nabtieh")] ="El Nabatieh"
survey $governorate[which(survey $governorate=="Beqaa")] ="Bekaa"
survey $governorate[which(survey $governorate=="Northern")] ="North"
survey $governorate[which(survey $governorate=="Southern")] ="South"
survey $governorate[which(survey $governorate=="Butnan")] ="Tobruk"
survey $governorate[which(survey $governorate=="Darnah")] ="Derna"
survey $governorate[which(survey $governorate=="Jabal Akhdar")] ="Al Jabal Al Akhdar"
survey $governorate[which(survey $governorate=="Marj")] = "Almarj"
survey $country[which(survey $governorate=="Bahariya")] = "Egypt"
survey $governorate[which(survey $governorate=="Bahariya")] = "Giza"
survey $governorate[which(survey $governorate=="Kufra")] = "Alkufra"
survey $governorate[which(survey $governorate=="Murqub")] ="Almargeb"
survey $governorate[which(survey $governorate=="Misurata")] ="Misrata"
survey $governorate[which(survey $governorate=="Ajafarh")] = "Aljfara"
survey $governorate[which(survey $governorate=="Zawiya")] ="Zwara"
survey $governorate[which(survey $governorate=="Nuqat al Khams")] ="Zwara"
survey $governorate[which(survey $governorate=="Jabal al Gharbi")] ="Al Jabal Al Gharbi"
survey $governorate[which(survey $governorate=="Sabha")] ="Sebha"
survey $governorate[which(survey $governorate=="Murzuk")] ="Murzuq"
survey $governorate[which(survey $governorate=="Wadi al Hayaa")] ="Ubari"
survey $governorate[which(survey $governorate=="Wadi al Shatii")] ="Wadi Ashshati"
survey $governorate[which(survey $governorate=="Jufra")] ="Aljufra"
survey $governorate[which(survey $governorate=="Sirte")] ="Sirt"
survey $governorate[which(survey $governorate=="Guélmim-Es Semara")] ="Guelmim Oued Noun" # new shapefile
survey $governorate[which(survey $governorate=="Laayoune-Boujdour-Sakia El Hamra")] ="Guelmim Oued Noun" #new shapefile
survey $governorate[which(survey $governorate=="Souss-Massa-Draâ")] = "Souss Massa"#new shapefile
survey $governorate[which(survey $governorate=="Gharb-Cherarda-Béni Hssen")] ="Rabat Sale Kenitra"#new shapefile
survey $governorate[which(survey $governorate=="Chaouia-Ouardigha")] ="Casablanca Settat"#new shapefile
survey $governorate[which(survey $governorate=="Marrakech-Tensift-Al Haouz")] ="Marrakech Safi"#new shapefile
survey $governorate[which(survey $governorate=="Grand-Casablanca")] ="Casablanca Settat"#new shapefile
survey $governorate[which(survey $governorate=="Rabat-Salé-Zemmour-Zaér")] ="Rabat Sale Kenitra"#new shapefile
survey $governorate[which(survey $governorate=="Tadla-Azilal")] ="Beni Mellal Khenifra"#new shapefile
survey $governorate[which(survey $governorate=="Méknès-Tafilalet")] ="Draa Tafilalet"
survey $governorate[which(survey $governorate=="Fès-Boulemane")] ="Fez Meknes"
survey $governorate[which(survey $governorate=="Taza-Al Hoceima-Taounate")] ="Tangier Tetouan Al Hoceima"
survey $governorate[which(survey $governorate=="Tanger-Tétouan")] ="Tangier Tetouan Al Hoceima"
survey $governorate[which(survey $governorate=="Oued Ed Dahab")] ="Guelmim Oued Noun"# this is in western sahara... assign to closest region
survey $governorate[which(survey $governorate=="Ben Arous")] = "Ben Arous (Tunis Sud)"
survey $governorate[which(survey $governorate=="Manouba")] = "Manubah"
survey $governorate[which(survey $governorate=="Zaghouane")] = "Zaghouan"
survey $governorate[which(survey $governorate=="Beja")] = "Béja"
survey $governorate[which(survey $governorate=="Kasserine")] = "Kassérine"
survey $governorate[which(survey $governorate=="Sidi Bouzid")] ="Sidi Bou Zid"
survey $governorate[which(survey $governorate=="Gabes")] ="Gabès"
survey $governorate[which(survey $governorate=="Mednine")] ="Médenine"
survey $governorate[which(survey $governorate=="al-Bayda'")] ="Al Bayda"
survey $governorate[which(survey $governorate=="Ta'izz")] ="Ta'iz"
survey $governorate[which(survey $governorate=="al-Jawf")] ="Al Jawf"
survey $governorate[which(survey $governorate=="al-Hudaydah")] ="Al Hodeidah"
survey $governorate[which(survey $governorate=="Hadhramaut")] ="Hadramawt"
survey $governorate[which(survey $governorate=="Saada")] ="Sa'dah"
survey $governorate[which(survey $governorate=="Lahij")] = "Lahj"
survey $governorate[which(survey $governorate=="al-Mahwit")] = "Al Mahwit"
survey $governorate[which(survey $governorate=="al-Mahrah")] = "Al Maharah"
survey $governorate[which(survey $governorate=="Ad Dali")] = "Ad Dali'"
survey $governorate[which(survey $governorate=="")] =NA
survey $country[which(survey $governorate=="Saudi Arabia")] = "Saudi Arabia"
survey $governorate[which(survey $governorate=="Saudi Arabia")] = NA#"Ar Riyad" 
#https://icsr.info/wp-content/uploads/2019/02/ICSR-Report-Saudi-Foreign-Fighters-Analysis-of-Leaked-Islamic-State-Entry-Documents.pdf
# most likely from 
survey $governorate[which(survey $governorate=="Syria")] = NA 
survey $governorate[which(survey $governorate=="Al Wahat")] = "Ejdabia"
survey $governorate[which(survey $governorate=="Laâyoune-Boujdour-Sakia El Hamra")] ="Guelmim Oued Noun"
survey $governorate[which(survey $governorate=="Béni Mellal-Khénifra")] = "Beni Mellal Khenifra"
survey $governorate[which(survey $governorate=="Tripoli" & survey $country=="Lebanon")] = "North"
# amman is just capital city in arabic
survey $governorate[which(survey $governorate=="Amman" & survey $country=="Kuwait")] = "Al Kuwayt"
# save clean individual-level data 

save(survey ,file = 'generated_data/Bird/clean_survey.RData')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # CLEAN OFFSET DENOMINATOR  # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# get denominator for adjustment 
denom = fread(file = 'data/mena_pops.csv')
# add saudi
denom = bind_rows(denom,data.frame(country = 'Saudi Arabia', 
                                   male_18plus = 15466126, 
                                   pct_sunni = 0.93*(1-0.125), # P(sunni, muslim) = P(muslim)*P(sunni|muslim) (1-0.125) = % sunni (amnogst muslims) - 0.93 = % muslim
                                   male_18_sunni= 15466126*0.125*0.93, 
                                   source = NA,
                                   year = NA,
                                   url = NA)) 
# saudi demo: https://www.stats.gov.sa/
# saudi rel : https://www.pewresearch.org/fact-tank/2018/04/12/5-facts-about-religion-in-saudi-arabia/
# add israel
denom = bind_rows(denom,data.frame(country = 'Israel', 
                                   male_18plus = 2607890,  
                                   pct_sunni = 0.82*0.21, # P(sunni, arab) = P(arab)*P(sunni|arab) 0.82  = % sunni amongst arabs - 0.21 = % arabs
                                   male_18_sunni= 0.82*0.21*2607890, 
                                   source = NA,
                                   year = NA,
                                   url = NA)) 
# https://www.cbs.gov.il/en/subjects/Pages/The-2008-Census-of-Population.aspx
# https://en.wikipedia.org/wiki/Demographics_of_Israel#:~:text=As%20of%202019%2C%20Arab%20citizens,Eastern%20Orthodox%20and%20Catholic%20denominations).
denom = denom [match(unique(shape$ADM0_EN),denom $country),]
denom $N1 = NA
denom $N1[denom $country=="Algeria"]=225
denom $N1[denom $country=="Egypt"]=1000
denom $N1[denom $country=="Jordan"]=2500
denom $N1[denom $country=="Kuwait"]=70
denom $N1[denom $country=="Lebanon"]=900
denom $N1[denom $country=="Libya"]=600
denom $N1[denom $country=="Morocco"]=1500
denom $N1[denom $country=="Tunisia"]=7000
denom $N1[denom $country=="Saudi Arabia"]=2500
denom $N1[denom $country=="Yemen"]= 110 # https://icsr.info/wp-content/uploads/2018/07/Women-in-ISIS-report_20180719_web.pdf
denom $N1[denom $country=="Israel"] = 60 # http://gppreview.com/2019/04/01/islamic-state-recruits-israel/
# # #
denom $N1_low = NA
denom $N1_low[denom $country=="Algeria"]=170
denom $N1_low[denom $country=="Egypt"]=600
denom $N1_low[denom $country=="Jordan"]=2000
denom $N1_low[denom $country=="Kuwait"]=70
denom $N1_low[denom $country=="Lebanon"]=900
denom $N1_low[denom $country=="Libya"]=600
denom $N1_low[denom $country=="Morocco"]=1200
denom $N1_low[denom $country=="Tunisia"]=6000
denom $N1_low[denom $country=="Yemen"]= 110 

save(denom ,file = 'generated_data/Bird/clean_denom.RData')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # GENERATE STAN DATA OBJECT # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# restrict analysis to complete cases
survey_complete = survey[complete.cases(survey),]   

# extract levels for country which actually exist in the dataset
survey_complete$governorate = as.factor(paste(survey_complete$country,survey_complete$governorate,sep=" - "))
levels(survey_complete$governorate) = c(levels(survey_complete$governorate),
                               paste(shape$ADM0_EN,shape$ADM1_EN,sep = " - ")[which(is.na(match(shape$ADM1_EN,levels(survey_complete$governorate))))])
survey_complete$country = as.factor(survey_complete$country)
levels(survey_complete$country) = c(levels(survey_complete$country),shape$ADM0_EN[which(is.na(match(shape$ADM0_EN,levels(survey_complete$country))))])
# disambiiguate shapefile names
shape$ADM1_EN = paste(shape$ADM0_EN,shape$ADM1_EN,sep = " - ")
# order shape same way as levels of survey
shape = shape[match(levels(survey_complete$governorate),shape$ADM1_EN),]
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# # # # # # # # GENERATE NEIGHBOURHOOD OBJECT # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# get distance matrix 
lat_lon = st_coordinates(st_centroid(shape))
colnames(lat_lon) = c("lat","lon")
D = geodist(lat_lon,measure = "geodesic") 

# # # district neighbours 
nb = poly2nb(as(shape,'Spatial'), row.names =shape$ADM1_EN, snap=sqrt(.Machine$double.eps),queen=TRUE, useC=TRUE, foundInBox=NULL)
nb = addnbs(sp.sample = as(shape,'Spatial'),ID = shape$ADM1_EN,D=D);

# borders between Algeria and Morocco
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tindouf",name2 = "Morocco - Guelmim Oued Noun",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tindouf",name2 = "Morocco - Souss Massa",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tindouf",name2 = "Morocco - Draa Tafilalet",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Bechar",name2 = "Morocco - Draa Tafilalet",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Bechar",name2 = "Morocco - Oriental",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Naama",name2 = "Morocco - Oriental",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tlemcen",name2 = "Morocco - Oriental",IDs = shape$ADM1_EN)
# borders between Algeria and Tunisia
nb = add_specific_nbs(nb = nb,name1 = "Algeria - El-Tarf",name2 = "Tunisia - Jendouba",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - El-Tarf",name2 = "Tunisia - Le Kef",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Souk-Ahras",name2 = "Tunisia - Jendouba",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Souk-Ahras",name2 = "Tunisia - Le Kef",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tebessa",name2 = "Tunisia - Le Kef",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tebessa",name2 = "Tunisia - Kassérine",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tebessa",name2 = "Tunisia - Gafsa",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Tebessa",name2 = "Tunisia - Tozeur",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - El Oued",name2 = "Tunisia - Kebili",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - El Oued",name2 = "Tunisia - Tataouine",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Ouargla",name2 = "Tunisia - Tataouine",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Illizi",name2 = "Tunisia - Tataouine",IDs = shape$ADM1_EN)
# borders between Tunisia and Libya
nb = add_specific_nbs(nb = nb,name1 = "Tunisia - Médenine",name2 = "Libya - Zwara",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Tunisia - Médenine",name2 = "Libya - Nalut",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Tunisia - Tataouine",name2 = "Libya - Zwara",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Tunisia - Tataouine",name2 = "Libya - Nalut",IDs = shape$ADM1_EN)
# borders between Libya and Egypt
nb = add_specific_nbs(nb = nb,name1 = "Libya - Tobruk",name2 = "Egypt - Matrouh",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Libya - Ejdabia",name2 = "Egypt - Matrouh",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Libya - Ejdabia",name2 = "Egypt - New Valley",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Libya - Alkufra",name2 = "Egypt - New Valley",IDs = shape$ADM1_EN)
# borders between Israel and Egypt
nb = add_specific_nbs(nb = nb,name1 = "Israel - HaDarom",name2 = "Egypt - North Sinai",IDs = shape$ADM1_EN)
# borders between Lebanon and Israel
nb = add_specific_nbs(nb = nb,name1 = "Lebanon - South",name2 = "Israel - HaZafon",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Lebanon - El Nabatieh",name2 = "Israel - HaZafon",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Lebanon - South",name2 = "Israel - Golan",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Lebanon - El Nabatieh",name2 = "Israel - Golan",IDs = shape$ADM1_EN)
# borders between Algeria and Libya
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Illizi",name2 = "Libya - Nalut",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Illizi",name2 = "Libya - Wadi Ashshati",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Illizi",name2 = "Libya - Ghat",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Algeria - Illizi",name2 = "Libya - Murzuq",IDs = shape$ADM1_EN)
# borders between Saudi and Yemen
nb = add_specific_nbs(nb = nb,name1 = "Saudi Arabia - Jizan",name2 = "Yemen - Hajjah",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Saudi Arabia - Jizan",name2 = "Yemen - Sa'dah",IDs = shape$ADM1_EN)

nb = add_specific_nbs(nb = nb,name1 = "Saudi Arabia - Najran",name2 = "Yemen - Sa'dah",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Saudi Arabia - Najran",name2 = "Yemen - Al Jawf",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Saudi Arabia - Najran",name2 = "Yemen - Hadramawt",IDs = shape$ADM1_EN)

nb = add_specific_nbs(nb = nb,name1 = "Saudi Arabia - Ash Sharqiyah",name2 = "Yemen - Hadramawt",IDs = shape$ADM1_EN)
nb = add_specific_nbs(nb = nb,name1 = "Saudi Arabia - Ash Sharqiyah",name2 = "Yemen - Al Maharah",IDs = shape$ADM1_EN)

# save this shape, distance matrix and neighbourhood matrix 
save(shape, file = 'generated_data/Bird/middle.east_shape.RData',compress = T)
save(nb,file = 'generated_data/Bird/middle.east_nb.RData',compress = T)
save(D,file = 'generated_data/Bird/middle.east_D.RData',compress = T)

# plot edges and nodes over map
pdf(file = 'Plots/Bird/connected.pdf',width = 5, height = 5)
plot_network(shape = shape, nb_object = nb)
dev.off()

# sample size
n = dim(survey_complete)[1]     

# vector of cases
y = survey_complete$case               


# known country-wide probability of recruitment 
# What Explains the Flow of Foreign Fighters to ISIS? Efraim Benmelech & Esteban F. Klor, 2020
pi_large_area = denom$N1[match(levels(survey_complete$country),denom$country)]/denom$male_18_sunni[match(levels(survey_complete$country),denom$country)]
names(pi_large_area ) = levels(survey_complete$country) 


# cases per area 
n1_large_area =  sapply(levels(survey_complete$country),function(x){
  sum(y[which(survey_complete$country==x)])
  })

# controls per area
n0_large_area = sapply(levels(survey_complete$country),function(x){
  sum(survey_complete$country==x) - sum(y[which(survey_complete$country==x)])
  })

# sample size per area
n_large_area = n1_large_area + n0_large_area
# calculate offset per area - no contamination case
offset_large_area = c(); 
for(i in 1:length(n1_large_area)){ 
  offset_large_area = c(offset_large_area, 
                        log( n1_large_area[i]/(pi_large_area[i]*n0_large_area[i]) + 1 ) )
}
# calculate offset per area - contamination
offset_large_area_no.contamination = c(); 
for(i in 1:length(n1_large_area)){ 
  offset_large_area_no.contamination = c(offset_large_area_no.contamination, 
                                         log(((1-pi_large_area[i])/pi_large_area[i])*( (n1_large_area[i]/n_large_area[i])  /(1-(n1_large_area[i]/n_large_area[i])))))
}


theta_0_large_area = rep(0,length(n1_large_area))
theta_1_large_area = n1_large_area/(n1_large_area+pi_large_area*n0_large_area)


# design-matrix 
X.original = survey_complete[,c("coledu","age","married","student","lowstat")]
X.original$age2 = X.original$age^2 
X.original$coledu_lowstat = X.original$coledu*X.original$lowstat
X = cbind(intercept = 1,scale(X.original))

# number of regression coefficients 
p = dim(X)[2]                  

# mean of covariates in original scale
mu_X = apply(X.original,2,mean)

# sd of covariates in original scale
sd_X  = apply(X.original,2,sd) 

# what kind of scaling is used
scale_X = "1sd"            

# area ids and names in the survey
small_area_id = as.integer(survey_complete$governorate)
small_area_names = levels(survey_complete$governorate)
large_area_id = as.integer(survey_complete$country)
large_area_names =  levels(survey_complete$country)

# number of large areas
N_large_area = nlevels(survey_complete$country)

# small-area neighborhood object 
# Have we achieved a fully connected graph ? 
if(!isDisconnected(nb)){print("Success! The Graph is Fully Connected")}else{print("Failure... Some parts of the graph are still disconnected...")}
nb_object = nb

# extract nb details
tmp = nb2graph(nb_object);
N_small_area =  tmp$N
node1_small_area = tmp$node1
node2_small_area = tmp$node2
N_small_area_edges = tmp$N_edges

# extract scaling factor
scaling_factor = scale_nb_components(nb_object)[1]     
# translate neighbourhood to matrix
nb.matrix = nb2mat(neighbours = nb_object,style = 'B',zero.policy = TRUE)  
# get small area shape
shape_small = shape       
# get large area shape
shape_large = shape %>% group_by(ADM0_EN) %>% summarise(geometry = sf::st_union(geometry)) %>% ungroup()
shape_large = shape_large[match(levels(survey_complete$country),shape_large$ADM0_EN),]

# generate data list for feeding to model 
data.list = 
         list(survey_complete = survey_complete,
              
              n = n,
              y = y,
              
              pi = NA,
              
              n1 = NA, 
              n0 = NA,
              offset = NA,
              theta_1 = NA,
              theta_0 = NA,
              
              pi_large_area = pi_large_area,
              
              n1_large_area = n1_large_area,
              n0_large_area = n0_large_area,
              
              offset_large_area = ifelse(is.na(offset_large_area),0,offset_large_area), #  the NAs are not going to be used by the model, but STAN won't accept a vector with missing values. 
              offset_large_area_no.contamination = ifelse(is.na(offset_large_area_no.contamination),0,offset_large_area_no.contamination), #  the NAs are not going to be used by the model, but STAN won't accept a vector with missing values. 
              
              theta_1_large_area = ifelse(is.na(theta_1_large_area),1,theta_1_large_area),
              theta_0_large_area = ifelse(is.na(theta_0_large_area),0,theta_0_large_area),
              
              
              
              X = X,
              p = p,
              mu_X = mu_X,
              sd_X = sd_X,
              scale_X = scale_X,
              
              small_area_id = small_area_id,
              small_area_names = small_area_names,
              large_area_id =  large_area_id,
              large_area_names = large_area_names,
              
              N_large_area = N_large_area,
              nb_object = nb_object,
              N_small_area = N_small_area,
              node1_small_area = node1_small_area,
              node2_small_area = node2_small_area ,
              N_small_area_edges =N_small_area_edges,
              scaling_factor = scaling_factor,
              nb.matrix  = nb.matrix ,
              
              shape_small = shape_small,
              shape_large = shape_large
)

save(data.list,file = 'generated_data/Bird/data.list.RData')






