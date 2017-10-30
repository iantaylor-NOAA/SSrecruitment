if (system("hostname", intern=TRUE) %in% c("NWCLW04223033") ){
  # location on Ian's computer
  mydir <- 'c:/github/SSrecruitment/'
}

# load r4ss package
require(r4ss)

# read example model into R
dir.ex1 <- file.path(mydir, 'models/rockfish_example')
ex1 <- SS_output(dir.ex1)
SS_plots(ex1)

# remake plot of data availability with larger scale
SSplotData(ex1)

dir.all_main <- file.path(mydir, 'models/rockfish_all_main')

# get function to copy input files
source(file.path(mydir, 'R/copy_SS_files.R'))


## create new model with all main recdevs 
if(FALSE){
  copy_SS_files(source=dir.ex1, target=dir.all_main)
  ctl.ex1 <- readLines(file.path(dir.ex1, 'rockfish_example_ctl.ss'))

  ctl.all_main <- ctl.ex1
  # change phases for early and late recdevs to be negative
  ctl.all_main[grep("early_phase", ctl.all_main)] <- " -5 #_recdev_early_phase"
  ctl.all_main[grep("forecast_recruitment", ctl.all_main)] <- " -5 #_forecast_recruitment phase"
  ctl.all_main[grep("first year of main recr_devs", ctl.all_main)] <- "1896 # first year of main recr_devs; early devs can preceed this era"
  ctl.all_main[grep("last year of main recr_devs", ctl.all_main)] <- "2016 # last year of main recr_devs; forecast devs start in following year"
  writeLines(ctl.all_main,
             file.path(dir.all_main, 'rockfish_example_ctl.ss'))
}

ex.all_main <- SS_output(dir.all_main)

dir.ex2 <- file.path(mydir, 'models/rockfish_example_bad_survey')
dir.ex3 <- file.path(mydir, 'models/rockfish_example_bad_survey_1955')
dir.all_main2 <- file.path(mydir, 'models/rockfish_all_main_bad_survey')

ex2 <- SS_output(dir.ex2)
ex3 <- SS_output(dir.ex3)
ex.all_main2 <- SS_output(dir.all_main2)

SSplotComparisons(SSsummarize(list(ex1, ex.all_main, ex2, ex.all_main2)))
SSplotComparisons(SSsummarize(list(ex2, ex3, ex.all_main2)))

rs1 <- SS_output(file.path(mydir, 'models/rockfish_simplified'))


### profile over sigmaR
dir.prof <- file.path(mydir, 'models/rockfish_example_sigmaR_profile')

# vector of values to profile over
sig.vec <- seq(1, 0, -.1)
Nprofile <- length(sig.vec)
# run SS_profile command
profile <- SS_profile(dir=dir.prof, # directory
                      model="ss",
                      masterctlfile="control.ss_new",
                      newctlfile="control_modified.ss",
                      string="SR_sigmaR",
                      extras = "-nox -nohess",
                      profilevec=sig.vec)

# read the output files (with names like Report1.sso, Report2.sso, etc.)
profilemodels <- SSgetoutput(dirvec=dir.prof, keyvec=1:Nprofile, getcovar=FALSE)
# replace final 2 models with manually run cases
profilemodels[[10]] <- SS_output(file.path(mydir,
                   'models/rockfish_example_sigmaR_0.1'))
profilemodels[[11]] <- SS_output(file.path(mydir,
                   'models/rockfish_example_sigmaR_0.01'))
profilesummary <- SSsummarize(profilemodels)

# plot sigmaR profile likelihood
SSplotProfile(profilesummary, profile.string="SR_sigmaR", legendloc='right',
              print=TRUE, plotdir=file.path(mydir, 'plots'),
              profile.label=expression(paste("Recruit deviation variability parameter: ",
                  sigma[italic(R)])))
recdev.sd.vec <- 0*sig.vec
for(i in 1:length(sig.vec)){
  recdev.sd.vec[i] <- profilemodels[[i]]$sigma_R_info$SD_of_devs[1]
}

# plot showing variability of recdevs as a function of sigmaR 
png(file.path(mydir, "plots/sigmaR_vs_recdev_variability.png"),
    width=6.5, height=5, pointsize=10, units='in', res=200, bg="transparent")
plot(sig.vec, recdev.sd.vec, lwd=4, col=rgb(0,0,1,.7), type='o', xaxs='i', yaxs='i',
     xlim=c(0,1.05), ylim=c(0,.6), las=1,
     xlab=expression(paste("Recruit deviation variability parameter: ", sigma[italic(R)])),
     ylab="S.D of 'main' recdevs (1960-2012)")
abline(0,1, lty=3)
dev.off()

### profile over h for B-H model
dir.prof <- file.path(mydir, 'models/rockfish_simplified_BH_steep_profile')

# vector of values to profile over
h.vec <- seq(0.2, 1.0, .1)
Nprofile <- length(h.vec)

# run SS_profile command
profile <- SS_profile(dir=dir.prof, # directory
                      model="ss",
                      masterctlfile="control.ss_new",
                      newctlfile="control_modified.ss",
                      string="steep",
                      extras = "-nox -nohess",
                      profilevec=h.vec)

# read the output files (with names like Report1.sso, Report2.sso, etc.)
profilemodels <- SSgetoutput(dirvec=dir.prof, keyvec=1:Nprofile, getcovar=FALSE)

png(file.path(mydir, "plots/Stock-recruit-curves_B-H.png"),
    width=5, height=5, pointsize=12, units='in', res=200, bg="transparent")
par(mar=c(4.1,4.1,1,1), las=1)
col.vec <- rich.colors.short(9)
for(i in seq(1,9,2)){
  SSplotSpawnrecruit(profilemodels[[i]], subplot=1, add=i>1,
                     colvec=c(NA, NA, col.vec[i], NA),
                     legend=FALSE, bias_adjusted=FALSE, estimated=FALSE,
                     relative=TRUE, init=FALSE, virg=FALSE)
  text(0.2, h.vec[i], paste0("h=",h.vec[i]), col=col.vec[i], adj=c(0,1))
}
dev.off()


### profile over Shepherd c with h = 0.6
dir.prof <- file.path(mydir, 'models/rockfish_simplified_Shep_h.6')
c.vec <- seq(.4,1.6,.2)
profile <- SS_profile(dir=dir.prof, # directory
                      model="ss",
                      masterctlfile="rockfish_simplified_ctl.ss",
                      newctlfile="rockfish_simplified_ctl.ss",
                      string="SR_Shepard_c",
                      extras = "-nox -nohess",
                      profilevec=c.vec)
profilemodels_Shep_h.6 <- SSgetoutput(dirvec=dir.prof,
                                      keyvec=1:Nprofile, getcovar=FALSE)
# plot of Shepherd stock-recruit curves
png(file.path(mydir, "plots/Stock-recruit-curves_Shep.png"),
    width=5, height=5, pointsize=12, units='in', res=200, bg="transparent")
par(mar=c(4.1,4.1,1,1), las=1)
col.vec <- rich.colors.short(6)
for(i in seq(6,1,-1)){
  SSplotSpawnrecruit(profilemodels_Shep_h.6[[i]], subplot=1, add=i!=6,
                     colvec=c(NA, NA, col.vec[i], NA),
                     legend=FALSE, bias_adjusted=FALSE, estimated=FALSE,
                     relative=TRUE, init=FALSE, virg=FALSE)
  text(0.2, h_adj(h=.6, c=c.vec[i]), paste0("c=",c.vec[i]),
       col=col.vec[i], adj=c(0,1))
}
dev.off()



### profile over h for Shepard with c=1
dir.prof <- file.path(mydir, 'models/rockfish_simplified_Shep_c1')
# run SS_profile command
profile <- SS_profile(dir=dir.prof, # directory
                      model="ss",
                      masterctlfile="rockfish_simplified_ctl.ss",
                      newctlfile="rockfish_simplified_ctl.ss",
                      string="steep",
                      extras = "-nox -nohess",
                      profilevec=h.vec)
# read the output files (with names like Report1.sso, Report2.sso, etc.)
profilemodels_c1 <- SSgetoutput(dirvec=dir.prof, keyvec=1:Nprofile, getcovar=FALSE)

### profile over h for Shepard with c=2
dir.prof <- file.path(mydir, 'models/rockfish_simplified_Shep_c2')
# run SS_profile command
profile <- SS_profile(dir=dir.prof, # directory
                      model="ss",
                      masterctlfile="rockfish_simplified_ctl.ss",
                      newctlfile="rockfish_simplified_ctl.ss",
                      string="steep",
                      extras = "-nox -nohess",
                      profilevec=h.vec)
# read the output files (with names like Report1.sso, Report2.sso, etc.)
profilemodels_c2 <- SSgetoutput(dirvec=dir.prof, keyvec=1:Nprofile, getcovar=FALSE)

### profile over h for Shepard with c=.5
dir.prof <- file.path(mydir, 'models/rockfish_simplified_Shep_c.5')
# run SS_profile command
profile <- SS_profile(dir=dir.prof, # directory
                      model="ss",
                      masterctlfile="rockfish_simplified_ctl.ss",
                      newctlfile="rockfish_simplified_ctl.ss",
                      string="steep",
                      extras = "-nox -nohess",
                      profilevec=h.vec)
# read the output files (with names like Report1.sso, Report2.sso, etc.)
profilemodels_c.5 <- SSgetoutput(dirvec=dir.prof, keyvec=1:Nprofile, getcovar=FALSE)


h_adj <- function(h, c){
  # calculation of adjusted h for Shepherd stock-recruit curve
  h_adj <- 0.2 + (h - 0.2)/0.8*(1/(5*0.2^c) - 0.2)
  return(h_adj)
}



