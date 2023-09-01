write_R_script_file=function(args_list,jobs_dir,logs_dir,jobname,scriptfile_path,job_time,mem,dryRun=FALSE){
log_dir_name=paste0(logs_dir,jobname)
script_filename=paste0(jobs_dir,jobname,".sh")
text_job=paste0("#!/bin/bash
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=",mem,"]
#BSUB -n 1
#BSUB -W ",job_time,":00
#BSUB -o ",log_dir_name,".stdout
#BSUB -e ",log_dir_name,".stderr
#BSUB -J ",jobname,"\n\n",

"ml R/4.3.0; Rscript --vanilla ",scriptfile_path," ",paste0(args_list,collapse=" "))
fileConn<-file(script_filename)
writeLines(text_job,fileConn)
print(script_filename)
if (dryRun==TRUE){
system(paste0("bsub<",script_filename))
}
}

mem=95000

args = commandArgs(trailingOnly=TRUE)

jobs_dir="/sc/arion/projects/psychAD/aging/sleep_patterns/analysis/dreamlet/jobs/"
logs_dir="/sc/arion/projects/psychAD/aging/sleep_patterns/analysis/dreamlet/logs/"

if (args[1] == "subtype") { 

file="/sc/arion/projects/CommonMind/aging/subtype_all_celltypes.txt"

} else if (args[1] == "subclass") {

file="/sc/arion/projects/CommonMind/aging/subclass_all_celltypes.txt"

} else {

""
}


celltypes_df=read.table(file)
celltypes=celltypes_df$V1

contrast_age_groups_list=c("early_adulthood","late_adulthood","old","adulthood","adulthood_old")

for (n_celltypes in celltypes) {

for (ii in c(1:length(contrast_age_groups_list))) {

scriptfile_path="/sc/arion/projects/psychAD/aging/sleep_patterns/scripts/day_night_age_groups_subid.R"
jobname=paste0(args[1],"_",n_celltypes,"_",contrast_age_groups_list[ii],"_subid")
args_list=c(args[1],n_celltypes,contrast_age_groups_list[ii])
job_time=48
write_R_script_file(args_list,jobs_dir,logs_dir,jobname,scriptfile_path,job_time,mem)

}

}


