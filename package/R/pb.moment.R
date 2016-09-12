pb.moment=function(Y.temp,cellsize){
    m1=sum(Y.temp)/sum(cellsize)
    m2=sum(Y.temp*(Y.temp-1))/sum(cellsize^2)
    m3=sum(Y.temp*(Y.temp-1)*(Y.temp-2))/sum(cellsize^3)
    kon.hat = -2*(-m1*m2^2+m1^2*m3)/(-m1*m2^2+2*m1^2*m3-m2*m3)
    koff.hat = 2*(m1^2-m2)*(m1*m2-m3)*(m1*m3-m2^2)/(m1^2*m2-2*m2^2+m1*m3)/(2*m1^2*m3-m1*m2^2-m2*m3)
    s.hat = (-m1*m2^2+2*m1^2*m3-m2*m3)/(m1^2*m2-2*m2^2+m1*m3)
    kinetic.estimates=c(round(kon.hat,4),round(koff.hat,4),round(s.hat,2))
    return(kinetic.estimates)
}