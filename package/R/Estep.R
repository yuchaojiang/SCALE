Estep = function(a, b, phi, epsilon,n,nA,nB){
  Zhat=matrix(ncol=4,nrow=n)
  colnames(Zhat)=c('Off','A','B','AB')
  for(j in 1:n){
    if(nA[j]<=1 & nB[j]<=1){
      prob1=1; prob2=prob3=prob4=0 # silent
    } else if (nA[j]>20 & nB[j]>20){
      prob4=1;prob1=prob2=prob3=0  # biallelic
    } else{
      prob1=epsilon^(nA[j]+nB[j])*phi[1]
      prob2=(1-epsilon)^(nA[j])*epsilon^(nB[j])*phi[2]
      prob3=epsilon^(nA[j])*(1-epsilon)^(nB[j])*phi[3]
      prob4=(1-epsilon)^(nA[j]+nB[j])*beta(nA[j]+a,nB[j]+b)/beta(a,b)*phi[4]
    }
    Zhat[j,]=c(prob1,prob2,prob3,prob4)/(prob1+prob2+prob3+prob4)
  }
  return(Zhat)
}