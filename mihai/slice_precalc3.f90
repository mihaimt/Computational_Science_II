 rad     = ratio_f(k,i)-ratio(k,i)*cos_diff(l,j)
 
 den     = rad*sqrt(rad)
 
 comm    = surf(k,l,i)/den

 acc(1) = acc(1)+comm*(proj_1(k,i,l)-cos_corner(j))
 
 acc(2) = acc(2)+comm*(proj_2(k,i,l)-sin_corner(j))

 rad     = ratio_f(k,i)+ratio(k,i)*cos_diff(l,j)

 den     = rad*sqrt(rad)

 comm    = surf(k,l,i)/den

 acc(3) = acc(3)+comm*(proj_1(k,i,l)+cos_corner(j))

 acc(4) = acc(4)+comm*(proj_2(k,i,l)+sin_corner(j))

