subroutine myprepot
!-----------------------------------------------------------------------------
!..
!.. The O(3P) + HBr --> OH + Br MRCI+Q/CBS(aug-cc-pVnZ(-PP); n = Q,5)+SO 3A"
!.. potential surface of A. G. S. de Oliveira-Filho, F. R. Ornellas
!.. and K. A. Peterson
!..
!.. Reference:
!.. Antonio G. S. de Oliveira-Filho, Fernando R. Ornellas and Kirk A. Peterson
!.. J. Chem. Phys. 136, 174316 (2012). 
!..
!.. The Reproducing Kernel Hilbert Space method of Ho and Rabitz is used
!.. to interpolate the surface consisting of 1110 geometries spanning
!.. O-H-Br angles of 60-180 deg.  This surface does not contain the
!.. H + BrO arrangement and no claims of accuracy are made for O-H-Br
!.. angles smaller than 60 deg. 
!-----------------------------------------------------------------------------
!  USAGE:
!
!  On Input:
!      r1 = R(O-H)
!      r2 = R(H-Br)
!      r3 = R(Br-O).
!
! On Output:
!      value = Energy, relative to the asymptotic reactants 
!              valley (O + HBr(r_e)), in hartree 
!      dr1 = Derivative of the potential with respect to r1 (R(O-H))
!            in hartree/bohr 
!      dr2 = Derivative of the potential with respect to r2 (R(H-Br))
!            in hartree/bohr 
!      dr3 = Derivative of the potential with respect to r3 (R(O-Br))
!            in hartree/bohr 
!
!  NOTE:
!    Before any actual potential energy calculations are made, a single
!    call to myprepot must be made:
!      call myprepot
!
!    Later, the potential energy is computed by calling mypot:
!      call mypot(r1,r2,r3,value,dr1,dr2,dr3)
!
!   The parameters are read from the o_hbr_a_pp.par file 
!
!-----------------------------------------------------------------------------
!
!                              Saddle point
!
!        r1                          r2                        r3
!   2.6289099999999999        2.8486750000000001        5.1064215904894388
!        energy
!   7.9866812204093529E-003
!        dr1                         dr2                       dr3
!   1.5320730578638475E-007   1.5109533511581397E-006  -1.5577335148869720E-007
!
!                             Reactants side vdW well
!
!        r1                          r2                        r3
!   4.3746320000000001        2.6855400000000000        7.0601719999996142     
!        energy
!  -2.5958471320488563E-003
!        dr1                         dr2                       dr3
!   1.4650519388234480E-003   1.4648796113941348E-003  -1.4650322069311608E-003
!
!                             Products  side vdW well
!
!        r1                          r2                        r3
!   1.8372599999999999        4.8222319999999996        4.4678939513580840     
!        energy
!  -3.3549004949264938E-002
!        dr1                         dr2                       dr3
!   6.5041529057150577E-007  -4.8876198921465885E-010   6.6615884752874166E-009
!
! ----------------------------------------------------------------------------
implicit none
integer, parameter:: dp=kind(0.d0)                   ! double precision
real(dp), dimension(14) :: phi_oh, phi_hbr 
real(dp), dimension(1110) :: x_1,x_2,y_1,amn
real(dp) :: r1,r2,r3,value,dr1,dr2,dr3
real(dp) :: r_oh,r_hbr,r_bro
real(dp) :: sum1, sum2, dzdr, dfdr 
real(dp) :: dbetadz, dedr_oh, dedr_hbr
real(dp) :: v_oh, v_hbr, v_bro, v3
real(dp) :: th_y,theta
real(dp) :: z,re,de,phi_inf,alpha_s,r_s,z_i,phi_z,f_s
real(dp) :: x_1_max,x_2_max,y_1_max
real(dp) :: x_1_min,x_2_min,y_1_min
real(dp) :: q_1,q_2,q_3,dq_1,dq_2,dq_3
integer :: i, ncall
data ncall/1/
save ncall
save phi_oh,phi_hbr,x_1,x_2,y_1,amn

if(ncall.eq.1) then 
    print *,"Loading O+HBr PES"
    print*, "Reference:"
    print*, "Antonio G. S. de Oliveira-Filho, Fernando R. Ornellas and Kirk A. Peterson"
    print*, "J. Chem. Phys. 136, 174316 (2012)." 
    open(unit=45,file="o_hbr_a_pp.par",status="old")
    do i=1,14
        read(45,*) phi_oh(i)
    end do
    do i=1,14
        read(45,*) phi_hbr(i)
    end do
    do i=1,1110
        read(45,*) x_1(i), x_2(i), y_1(i), amn(i)
    end do
    close(45)
    ncall=2
    return
end if
 
entry mypot(r1,r2,r3,value,dr1,dr2,dr3)

r_oh=r1
r_hbr=r2
r_bro=r3

re=phi_oh(1)
de=phi_oh(2)
phi_inf=phi_oh(3)
alpha_s=phi_oh(4)
r_s=phi_oh(5)
z=(r_oh-re)/(r_oh+re)
phi_z=-phi_inf
z_i=1.0_dp
sum2=0.0_dp
do i=6,14,1
    phi_z=phi_z+phi_oh(i)*z_i
    sum2=sum2+(i-6)*phi_oh(i)*z_i/z
    z_i=z_i*z
end do
sum1=phi_z
f_s=1.0_dp/(exp(alpha_s*(r_oh-r_s))+1.0_dp)
phi_z=f_s*phi_z+phi_inf
v_oh=de*((1.0_dp-(re/r_oh)**6*exp(-phi_z*z))**2-1.0_dp)

dzdr=2.0_dp*re/((r_oh+re)**2)
dfdr=-alpha_s*exp(alpha_s*(r_oh-r_s))*f_s*f_s

dbetadz=phi_z+((1/dzdr)*dfdr*sum1+f_s*sum2)*z
dedr_oh=2.0_dp*de*(1.0_dp-(re/r_oh)**6*exp(-phi_z*z))*((re/r_oh)**6)*exp(-phi_z*z)*(dzdr*dbetadz+6.0_dp/r_oh)

re=phi_hbr(1)
de=phi_hbr(2)
phi_inf=phi_hbr(3)
alpha_s=phi_hbr(4)
r_s=phi_hbr(5)
z=(r_hbr-re)/(r_hbr+re)
phi_z=-phi_inf
z_i=1.0_dp
sum2=0.0_dp
do i=6,14,1
    phi_z=phi_z+phi_hbr(i)*z_i
    sum2=sum2+(i-6)*phi_hbr(i)*z_i/z
    z_i=z_i*z
end do
sum1=phi_z
f_s=1.0_dp/(exp(alpha_s*(r_hbr-r_s))+1.0_dp)
phi_z=f_s*phi_z+phi_inf
v_hbr=de*((1.0_dp-(re/r_hbr)**6*exp(-phi_z*z))**2-1.0_dp)

dzdr=2.0_dp*re/((r_hbr+re)**2)
dfdr=-alpha_s*exp(alpha_s*(r_hbr-r_s))*f_s*f_s

dbetadz=phi_z+((1/dzdr)*dfdr*sum1+f_s*sum2)*z
dedr_hbr=2.0_dp*de*(1.0_dp-(re/r_hbr)**6*exp(-phi_z*z))*((re/r_hbr)**6)*exp(-phi_z*z)*(dzdr*dbetadz+6.0_dp/r_hbr)

v_bro=0.0868406_dp*exp(-1.91_dp*(r_bro-3.2509_dp))

th_y=(1.0_dp-((r_oh*r_oh+r_hbr*r_hbr-r_bro*r_bro)/(2.0_dp*r_oh*r_hbr)))/2.0_dp

v3=0.0_dp
dr1=0.0_dp
dr2=0.0_dp
dr3=0.0_dp
do i=1,1110,1

    if (r_oh.ge.x_1(i)) then
        q_1=(1.0_dp/14.0_dp)*(r_oh**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(x_1(i)/r_oh))
        dq_1=0.5_dp*(r_oh**(-8))*((8.0_dp/9.0_dp)*(x_1(i)/r_oh)-1.0_dp)
    else
        q_1=(1.0_dp/14.0_dp)*(x_1(i)**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(r_oh/x_1(i)))
        dq_1=-x_1(i)**(-8)/18.0_dp
    end if

    if (r_hbr.ge.x_2(i)) then
        q_2=(1.0_dp/14.0_dp)*(r_hbr**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(x_2(i)/r_hbr))
        dq_2=0.5_dp*(r_hbr**(-8))*((8.0_dp/9.0_dp)*(x_2(i)/r_hbr)-1.0_dp)
    else
        q_2=(1.0_dp/14.0_dp)*(x_2(i)**(-7))*(1.0_dp-(7.0_dp/9.0_dp)*(r_hbr/x_2(i)))
        dq_2=-x_2(i)**(-8)/18.0_dp
    end if

    if (th_y.ge.y_1(i)) then
        q_3=1.0_dp+th_y*y_1(i)+2.0_dp*y_1(i)**2*th_y*(1.0_dp-(1.0_dp/3.0_dp)*(y_1(i)/th_y))
        dq_3=y_1(i)+2.0_dp*y_1(i)**2
    else
        q_3=1.0_dp+th_y*y_1(i)+2.0_dp*th_y**2*y_1(i)*(1.0_dp-(1.0_dp/3.0_dp)*(th_y/y_1(i)))
        dq_3=y_1(i)+4.0_dp*y_1(i)*th_y-2.0_dp*th_y**2
    end if

    v3=v3+amn(i)*q_1*q_2*q_3
    dr1=dr1+amn(i)*dq_1*q_2*q_3
    dr2=dr2+amn(i)*q_1*dq_2*q_3
    dr3=dr3+amn(i)*q_1*q_2*dq_3
end do
value=0.97714_dp*v3+v_oh+v_hbr+v_bro+0.14242006_dp

dr3=0.97714_dp*dr3*r3/(2.0_dp*r1*r2)
dr1=0.97714_dp*dr1-dr3*((r1-r2*((r_oh**2+r_hbr**2-r_bro**2)/(2*r_oh*r_hbr)))/r3)
dr2=0.97714_dp*dr2-dr3*((r2-r1*((r_oh**2+r_hbr**2-r_bro**2)/(2*r_oh*r_hbr)))/r3)

dr1=dedr_oh+dr1
dr2=dedr_hbr+dr2
dr3=dr3-1.91_dp*0.0868406_dp*exp(-1.91_dp*(r_bro-3.2509_dp))
end subroutine
