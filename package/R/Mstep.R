Mstep = function (Zhat,cura,n,nA,nB){
  phihat=apply(Zhat,2,sum)/n
  phihat=pmax(phihat,0.01)
  phihat=phihat/sum(phihat)
  a.b=optim(par=cura,
            function(x){a=b=x;
  l=sum(Zhat[,4]*(lbeta(nA+a,nB+b)-lbeta(a,b)))
  return(-l)},
  method='Brent',lower=3,upper=100)$par
  par.list=list(phihat=phihat,a.b=a.b)
  return(par.list)
}