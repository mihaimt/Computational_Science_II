 rad     = ratio_f(k,i) -ratio(k,i)*cos_diff(l,j)

 den     = rad * sqrt(rad)

 den_inv = 1./den

 comm    = mass(k,l)*radius_corn_2_inv(i)*den_inv
 
 acc(1) = acc(1)+comm*(ratio(k,i)*cos_center(l)-cos_corner(j))

 acc(2) = acc(2)+comm*(ratio(k,i)*sin_center(l)-sin_corner(j))
            
 rad     = ratio_f(k,i)+ratio(k,i)*cos_diff(l,j)

 den     = rad*sqrt(rad)

 den_inv = 1./den

 comm    = mass(k,l)*radius_corn_2_inv(i)*den_inv

 acc(3) = acc(3)+comm*(ratio(k,i)*cos_center(l)+cos_corner(j))

 acc(4) = acc(4)+comm*(ratio(k,i)*sin_center(l)+sin_corner(j))

