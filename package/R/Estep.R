Estep = function(a, b, phi, epsilon,n,nA,nB){
  
  prob1=epsilon^(nA+nB)*phi[1]
  prob2=(1-epsilon)^nA*epsilon^nB*phi[2]
  prob3=epsilon^nA*(1-epsilon)^nB*phi[3]
  prob4=(1-epsilon)^(nA+nB)*beta(nA+a,nB+b)/beta(a,b)*phi[4]
  
  hardcall1 = ((nA<=1) & (nB<=1))
  prob1[hardcall1]=1; prob2[hardcall1]=0
  prob3[hardcall1]=0; prob4[hardcall1]=0
  
  hardcall2 = ((nA>=20) & (nB>=20))
  prob1[hardcall2]=0; prob2[hardcall2]=0
  prob3[hardcall2]=0; prob4[hardcall2]=1
  Zhat=cbind(prob1,prob2,prob3,prob4)
  Zhat=Zhat/matrix(ncol=ncol(Zhat),nrow=nrow(Zhat),data=apply(Zhat,1,sum))
  colnames(Zhat)=c('Off','A','B','AB')
  return(Zhat)
}