!投影法
!方腔100cm×100cm
!------------------------------------------------------
!U-x方向速度分量,V-y方向速度
!P-压力值
!F-流函数
!Vor-涡量
!av——迎风特性系数  w——松弛因子
program projetion_method
implicit none
real*8, allocatable::U(:,:),Ut(:,:),V(:,:),Vt(:,:),P(:,:),Pt(:,:),d(:,:),Sp(:,:),f(:,:),vor(:,:)
real*8 max,con,dt,dx,dy,at,bt,av,w,re,max_du,max_dv,error
integer i,j,a,b,time,k,n
!--------------------------------------------划分网格,再扩充两个半网格边界外的虚拟网格
!V计算点——(i+1/2,j)
!U计算点——(i,j+1/2)
!P计算点——(i+1/2,j+1/2)
!F,vor计算点——(i,j)
a=128	
b=128
w=0.9
av=0.8
re=10000.
re=1./re
!dt=0.00001
!dt=1./a/b
dx=1./a;dy=1./b
dt=0.05*dx
at=1./dx/dx
bt=1./dy/dy
allocate( U(-a-2:a+2,-b-2:b+2),Ut(-a-2:a+2,-b-2:b+2),V(-a-2:a+2,-b-2:b+2),Vt(-a-2:a+2,-b-2:b+2) )
allocate( P(-a-2:a+2,-b-2:b+2),Pt(-a-2:a+2,-b-2:b+2),D(-a-2:a+2,-b-2:b+2),Sp(-a-2:a+2,-b-2:b+2) )
allocate( F(-a-2:a+2,-b-2:b+2),Vor(-a-2:a+2,-b-2:b+2) )
!*************************************************************
U=0.
V=0.
P=0.
U(:,b)=1.	!初始值
u(:,b+1)=2.*U(:,b)-u(:,b-1)
!----------------------------------边界条件----------------
!AD&BC侧边:p(i+1,j)-p(i,j)=0
!			U=V=0,  U(i+1,j)=-U(i,j), V(i,j)=V(i,J+1)=0
!AB&CD上下边：P(i,j+1)-P(i,j)=0
!			U+V=0,  U(i,j)=U(i+1,j)=0, V(i,j)=-V(i,j+1)
!壁面压力梯度近似为0
!---------------------------------------
call vir_mesh	!虚拟网格数值处理
!==============================================================
time=0

do while( 1.>0.)

max=0.

call eval_U  	!预估速度V(*)

call vir_mesh

call cal_P		!由上面那部分速度计算P(n+1),迭代出解！！！

call cal_U		!计算速度新值

call vir_mesh
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(mod(time,200)==0)then
	k=time/200
call position_UV  !求节点速度
call cal_stream !计算流函数
call cal_vorticity !计算涡量
call position_P  !计算压力在网格结点上的值

	open(10,file='data_'//char(48+mod(k/1000,10))//char(48+mod(k/100,10))//char&
	&(48+mod(k/10,10))//char(48+mod(k,10))//'.dat',status='unknown')
	write(10,"('TITLE=""grid=',I3,',Re=',f6.1,'"" ')")a,1./re
	write(10,*)'VARIABLES= "X" , "Y" ,"U","v","Ut","P","f","vor"'
	write(10,*)'ZONE I=',a+1,', J=',b+1,' F=POINT'

	do i=-a*0.5,a*0.5
	do j=-b*0.5,b*0.5
		write(10,*)i+0.5*a,j+0.5*b,U(2.*i,2.*j),v(2.*i,2.*j),ut(2.*i,2.*j),P(2.*i,2.*j),F(2.*i,2.*j),vor(2.*i,2.*j)
	enddo
	enddo
	close(10)


endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

time=time+1

!-------------------------判断是否达到计算时刻或者稳态收敛----
error=0.
do i=-a+1,a-1,2
do j=-b+1,b-1,2
	con= (U(i+1,j)-U(i-1,j))/dx +(v(i,j+1)-v(i,j-1))/dy 
	error=dmax1(error,abs(con) )
	!if(max<abs(con)) max=abs(con)
enddo
enddo
if (time>500)  then
if(mod(time,200)==0)then
write(*,*) time,error !,max_du,max_dv
end if
else 
print*,time,error
endif
!if(max<0.00000001) exit

enddo

call position_UV  !求节点速度
call cal_stream !计算流函数
call cal_vorticity !计算涡量
call position_P  !计算压力在网格结点上的值
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	open(10,file='output_f.dat')
	write(10,*)'TITLE=""grid=',a,',re=',1./re,'""'
	write(10,*)'VARIABLES= "X" , "Y" ,"U","v","Ut","P","f","vor"'
	write(10,*)'ZONE I=',a+1,', J=',b+1,' F=POINT'

	do i=-a*0.5,a*0.5
	do j=-b*0.5,b*0.5
		write(10,*)i+0.5*a,j+0.5*b,U(2.*i,2.*j),v(2.*i,2.*j),ut(2.*i,2.*j),P(2.*i,2.*j),F(2.*i,2.*j),vor(2.*i,2.*j)
	enddo
	enddo
	close(10)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

deallocate(Ut,U,V,Vt,P,Pt)
10 format(1x,f10.5)
contains
!*************************************************************
subroutine vir_mesh

	!AD 
	U(-a-2,:)=-1*u(-a+2,:)
	V(-a-1,:)=-1*v(-a+1,:)
	!BC  
	v(a+1,:)=-v(a-1,:)
	U(a+2,:)=-U(a-2,:)
	!CD  
	U(:,-b-1)=-U(:,-b+1)
	v(:,-b-2)=-v(:,-b+2)
	!AB 
	v(:,b+2)=-v(:,b-2)
	U(:,b+1)=2.*U(:,b)-U(:,b-1)

end subroutine vir_mesh
!***************************************************************
subroutine eval_U 		!忽略压力项，预估速度

real*8 FUX,FUY,FVX,FVY
do i=-a+1,a-1,2
do j=-b+1,b-1,2

	!FUX=-0.25*dt/dx*( (U(i+3,j)+U(i+1,j))**2 -(U(i-1,j)+U(i+1,j))**2 +av*abs(U(i+3,j)+U(i+1,j))*(U(i+1,j)-U(i+3,j)) -av*abs(U(i-1,j)+U(i+1,j))*(U(i-1,j)-U(i+1,j)) )
	FUX=-0.25*dt/dx*( U(i+3,j)*(U(i+3,j)+2.*U(i+1,j)) -U(i-1,j)*(U(i-1,j)+2.*U(i+1,j)) +av*abs(U(i+3,j)+U(i+1,j))*(U(i+1,j)-U(i+3,j)) -av*abs(U(i-1,j)+U(i+1,j))*(U(i-1,j)-U(i+1,j)) )
	FUY=-0.25*dt/dy*( (U(i+1,j)+U(i+1,j+2))*(v(i,j+1)+v(i+2,j+1)) +av*abs(v(i,j+1)+v(i+2,j+1))*(U(i+1,j)-U(i+1,j+2)) -(U(i+1,j-2)+U(i+1,j))*(v(i,j-1)+v(i+2,j-1)) -av*abs(v(i,j-1)+v(i+2,j-1))*(U(i+1,j-2)-U(i+1,j)) )
	Ut(i+1,j)= U(i+1,j)+FUX+FUY+dt*at*( U(i+3,j)-2.*U(i+1,j)+U(i-1,j) )*re +dt*bt*( U(i+1,j+2)-2.*U(i+1,j)+U(i+1,j-2) )*re
	
	!FVX=-0.25*dt/dy*( (v(i,j+3)+v(i,j+1))**2 +av*abs(v(i,j+3)+v(i,j+1))*(v(i,j+1)-v(i,j+3)) -(v(i,j-1)+v(i,j+1))**2 -av*abs(v(i,j-1)+v(i,j+1))*(v(i,j-1)-v(i,j+1)) )
	FVX=-0.25*dt/dy*( v(i,j+3)*(v(i,j+3)+2*v(i,j+1)) -v(i,j-1)*(v(i,j-1)+2.*v(i,j+1)) +av*abs(v(i,j+3)+v(i,j+1))*(v(i,j+1)-v(i,j+3)) -av*abs(v(i,j-1)+v(i,j+1))*(v(i,j-1)-v(i,j+1)) )
	FVY=-0.25*dt/dx*( (U(i+1,j+2)+U(i+1,j))*(v(i+2,j+1)+v(i,j+1)) +av*abs(U(i+1,j+2)+U(i+1,j))*(v(i,j+1)-v(i+2,j+1)) -(u(i-1,j+2)+u(i-1,j))*(v(i,j+1)+v(i-2,j+1)) -av*abs(u(i-1,j+2)+u(i-1,j))*(v(i-2,j+1)-v(i,j+1)) )
	vt(i,j+1)= v(i,j+1)+FVX+FVY+dt*at*( v(i+2,j+1)-2.*v(i,j+1)+v(i-2,j+1) )*re +dt*bt*( v(i,j+3)-2.*v(i,j+1)+v(i,j-1) )*re

enddo
enddo
Ut(:,b)=1.; Ut(:,-b)=0.; Ut(-a,:)=0.; Ut(a,:)=0.
vt(:,b)=0.; vt(:,-b)=0.; vt(-a,:)=0.; vt(a,:)=0.
U=Ut
v=vt

end subroutine eval_U
!**************************************************************
subroutine cal_P
do n=1,300
max=0.
do i=-a+1,a-1,2
do j=-b+1,b-1,2
	con= (U(i+1,j)-U(i-1,j))/dx +(v(i,j+1)-v(i,j-1))/dy 
	Pt(i,j)=0.5*( at*(P(i+2,j)+P(i-2,j))+ bt*(P(i,j+2)+P(i,j-2))- con/dt )/(at+bt)
	Pt(i,j)=w*Pt(i,j)+(1-w)*P(i,j)
	max = dmax1(max,abs(Pt(i,j)-P(i,j)))
	if(max<abs(Pt(i,j)-P(i,j)))  max=abs(Pt(i,j)-P(i,j))
	P(i,j)=Pt(i,j)
	!P(i,j)=w*0.5*( at*(P(i+2,j)+P(i-2,j))+ bt*(P(i,j+2)+P(i,j-2))- con/dt )/(at+bt) +(1-w)*P(i,j)
enddo
enddo  
if(max<1.E-007) exit   
P(-a-1,:)=P(-a+1,:)
P(a+1,:)=P(a-1,:)
P(:,-b-1)=P(:,-b+1)
P(:,b+1)=P(:,b-1)

enddo

end subroutine cal_P
!********************************************************************************
subroutine cal_U

max_du=0;max_dv=0.
do i=-a+1,a-1,2
do j=-b+1,b-1,2
U(i+1,j)=U(i+1,j) -dt/dx*(P(i+2,j)-P(i,j))
v(i,j+1)=v(i,j+1) -dt/dy*(P(i,j+2)-P(i,j)) 
enddo
enddo
U(a,:)=0.;v(:,b)=0.
U(-a,:)=0.;v(:,-b)=0.
end subroutine cal_U
!**********************************************************************************
!*********************************************************
subroutine cal_stream

real fu(-a-2:a+2,-b-2:b+2),fv(-a-2:a+2,-b-2:b+2)
fu= 0.          
fv= 0.
f(-a,:)=0.;f(a,:)=0.;f(:,-b)=0.;f(:,b)=0.
do i=-a+2,a-2,2
do j=-b+2,b-2,2
	fu(i,j)= dy*u(i,j-1)+fu(i,j-2)
end do
end do    
       
do j=-a+2,a-2,2
do i=-b+2,b-2,2
	fv(i,j)= -1.0*dx*v(i-1,j)+fv(i-2,j)
end do
end do    
       
do i=-a+2,a-2,2
do j=-b+2,b-2,2
	!f(0.5*i,0.5*j)= 0.5*(fu(i,j)+fv(i,j))
	 f(i,j)= 0.5*(fu(i,j)+fv(i,j))
end do
end do    

end subroutine cal_stream 

!*****************************************************************
subroutine cal_vorticity

vor=0.
do i=-a+2,a-2,2        !求涡量函数
do j=-b+2,b-2,2   
	!vor(0.5*i,0.5*j)= (v(i+1,j)-v(i-1,j))/dx - (u(i,j+1)-u(i,j-1))/dy
	vor(i,j)= (v(i+1,j)-v(i-1,j))/dx - (u(i,j+1)-u(i,j-1))/dy
end do
end do    
    
                       !求边界上的涡量函数
     do j=-b+2,b-2,2
       vor(-a,j)= (v(-a+1,j)-v(-a-1,j))/dx   !left wall
     end do
     
     do j=-b+2,b-2,2
       vor(a,j)= (v(a+1,j)-v(a-1,j))/dx    !right wall
     end do
     
     do i=-a+2,a-2,2
       vor(i,b)= (u(i,b+1)-u(i,b-1))/dy    !up wall
     end do
     
     do i=-a+2,a-2,2
       vor(i,-b)= (u(i,-b+1)-u(i,-b-1))/dy    !down wall
     end do
     
                     !计算角点涡量函数
                     
      vor(-a,-b)= 0.5*(vor(-a,-b+1) +  vor(-a+1,-b))
      vor(-a,b)= 0.5*(vor(-a,b-1) +  vor(-a+1,b))
      vor(a,-b)= 0.5*(vor(a-1,-b) +  vor(a,-b+1))
      vor(a,b)= 0.5*(vor(a,b-1) +  vor(a-1,b))
      
end subroutine cal_vorticity

!************************************************************
subroutine position_UV

do j=-a+2,a-2,2          !计算u 速度在网格节点上的值
do i=-b+2,b-2,2
U(i,j)=(u(i,j+1)+u(i,j-1))*0.5  
end do
end do
 
do i=-a+2,a-2,2         !计算v 速度在网格节点上的值
do j=-b+2,b-2,2
v(i,j)=(v(i+1,j)+v(i-1,j))*0.5
end do
end do   
UT=0. 
do j=-a+2,a-2,2           !求节点速度绝对值
do i=-b+2,b-2,2
ut(i,j)= sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j))
end do
end do   

end subroutine position_UV

!*************************************************************
subroutine position_P
							!计算压力在网格结点上的值
do i=-a,a,2
do j=-a,b,2
P(i,j)=0.25*(p(i+1,j+1)+p(i-1,j+1)+p(i-1,j-1)+p(i+1,j-1))
enddo	
enddo

end subroutine position_P

!*************************************************************
end
