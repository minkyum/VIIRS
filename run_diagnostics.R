for(ss in 1:6){
  for(tt in c(73,72)){
    system(paste('qsub -V -pe omp 8 -l h_rt=04:00:00 /usr3/graduate/mkmoon/GitHub/VIIRS/run_diagnostics.sh ',ss,tt,sep=''))  
  }  
}




