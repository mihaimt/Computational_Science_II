! calculate forcel______________________________________________________

level=2
do j=1,dim_theta
	do i=1,dim_r
		dforce_r=0.0
		dforce_theta=0.0

		if (Mod(i,2) == 1 .and. Mod(j,2)==1) then
			! level 0 
			do l=j-2,j+1
				if(l>0.and.l<=dim_theta) then
					var2=cos_dtheta0(j,l)
					var3=sin_dtheta0(j,l)
					do k=i-2,k+1
						if(k>0.and.k<=dim_r) then
							var4=r_sub0(i,k)
							var5=1+var4*var4-2*var4*var2
							var5=var5*sqrt(var5)
							var6=mass_term0(k,l)/(var5)
							dforce_r=dforce_r+(var4-var2)*var6
							dforce_theta=dforce_theta+var3*var6
						end if
					end do
				end if
			end do
			! level 1
			do l=j-3,j+3,2
				if(l>0.and.l<=dim_theta) then
					var2=cos_dthetal(j,l)
					var3=sin_dthetal(j,l)
					do k=i-3,i+3,2
						if(.not.(abs(l)==2.and.abs(k)==2)) then
							if(k>0.and.k<=dim_r) then
								var4=r_subl(i,k)
								var5=1+var4*var4-2*var4*var2
								var5=var5*sqrt(var5)
								var6=mass_terml(1,k,l)/(var5)
								dforce_r=dforce_r+(var4-var2)*var6
								dforce_theta=dforce_theta+var3*var6
							end if
						end if
					end do
				end if
			end do
			! level 2
			do l=j-6-14*4,j+6+14*4,4
				if(l>0.and.l<=dim_theta) then
					var2=cos_dthetal(j,l)
					var3=sin_dthetal(j,l)
					do k=i-6-14*4,i+6+14*4,4
						if(.not.(abs(l)==2.and.abs(k)==2))
							if(k>0.and.k<=dim_r.and.(l) then
								var4=r_subl(i,k)
								var5=1+var4*var4-2*var4*var2
								var5=var5*sqrt(var5)
								var6=mass_terml(2,k,l)/(var5)
								dforce_r=dforce_r+(var4-var2)*var6
								dforce_theta=dforce_theta+var3*var6
							end if
						end if
					end do
				end if
			end do
						
			
					
					

			do l=2,dim_theta,2
				var2=cos_dthetal(j,l)
				var3=sin_dthetal(j,l)
				do k=2,dim_r,2
					var4=r_subl(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_terml(level,k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			forcel(i,j,1)=dforce_r
			forcel(i,j,2)=dforce_theta
		end if

		if (Mod(i,2) == 0 .and. Mod(j,2)==1) then
			do l=2,dim_theta,2
				var2=cos_dthetal(j,l)
				var3=sin_dthetal(j,l)
				do k=1,dim_r,2
					var4=r_subl(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_terml(level,k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			forcel(i,j,1)=dforce_r
			forcel(i,j,2)=dforce_theta
		end if

		if (Mod(i,2) == 1 .and. Mod(j,2)==0) then
			do l=1,dim_theta,2
				var2=cos_dthetal(j,l)
				var3=sin_dthetal(j,l)
				do k=2,dim_r,2
					var4=r_subl(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_terml(level,k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			forcel(i,j,1)=dforce_r
			forcel(i,j,2)=dforce_theta
		end if

		if (Mod(i,2) == 0 .and. Mod(j,2)==0) then
			do l=1,dim_theta,2
				var2=cos_dthetal(j,l)
				var3=sin_dthetal(j,l)
				do k=1,dim_r,2
					var4=r_subl(i,k)
					var5=1+var4*var4-2*var4*var2
					var5=var5*sqrt(var5)
					var6=mass_terml(level,k,l)/(var5)
					dforce_r=dforce_r+(var4-var2)*var6
					dforce_theta=dforce_theta+var3*var6
				end do
			end do
			forcel(i,j,1)=dforce_r
			forcel(i,j,2)=dforce_theta
		end if
	end do
end do
