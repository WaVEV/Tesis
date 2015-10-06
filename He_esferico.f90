! Codigo que calculos los autovalores del Helio esfercio.
! El hamiltoniano es:
! H = h(1) + h(2) + lambda* 1/r_>
! donde h(i) = -1/2 d^2 /dr^2 + 1/r
! 
! El codigo calcula los autovalores para el estado singlete (espin antiparalelo) por lo tanto 
! las funcion de onda espacial es simetrica. Una vez calculado el hamiltoniano en la base producto
! hace una simetrizacion, se reduce el tamaño de la matriz.
! 
! El codigo usa B-splines para calcular el hamiltoniano y luega calcula los 
! autovalores y autovectores que elija.
! 
! El intervalo de integracion es [0, Rmax].
! 
! El codigo va variando la variable eta, por lo tanto devuelve los autovalores como 
! funcion de eta. Eta representa la intesidad del potencial Vef
!!!!!!!!!!!!!!!!!!!!!! preprocessing variable !!!!!!!

#ifndef Rmax_ 
#define Rmax_ 10.d0
#endif

#ifndef lum_
#define lum_ 0
#endif

#ifndef c_
#define c_ 0.d0
#endif

#ifndef gamma_
#define gamma_ 0.d0
#endif

#ifndef kord_
#define kord_ 5
#endif

#ifndef l_
#define l_ 30
#endif

#ifndef intg_
#define intg_ 100
#endif

#ifndef nev_
#define nev_ 10
#endif

#ifndef etai_
#define etai_ 0.1d0
#endif

#ifndef etaf_
#define etaf_ 1.d0
#endif

#ifndef num_puntos_eta_
#define num_puntos_eta_ 10
#endif

#ifndef me_
#define me_ 1.d0
#endif

!!!!!!!!!!!!!!!!!!!!!! MODULES !!!!!!!!!!!!!!!!!!!!!!!
 
module carga
integer :: kord, lum, intg, nev, num_puntos_eta
real(8) :: Rmin, Rmax, eta, etai, etaf
real(8) :: me, delta, c, gamma, ll
integer :: l
character(1) :: tip, l1, l2
end module  carga

module matrices
integer :: nk, nb
real(8), allocatable, dimension(:) :: norma
real(8), allocatable, dimension(:,:) :: s, v01, ke, aux
real(8), allocatable, dimension(:,:,:,:) :: Vef
integer, allocatable, dimension(:,:) :: ist, irt, tri
end module matrices

module integracion
integer, allocatable, dimension(:) :: k
real(8), allocatable, dimension(:) :: t, sp, dsp
real(8), allocatable, dimension(:,:) :: x, w, pl
end module  integracion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!   MAIN !!!!!!!!!!!!!!!!!!!!!!!!!
program Bsplines
use carga
use matrices
use integracion
implicit none
integer :: i, j
real(8) :: time, t1, t2, t3, t4, t5, t6
real(8) :: omp_get_wtime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Rmin = 0.d0

Rmax = 10.0 ! valores de inicio y final del intervalo de integracion

tip = 'u' ! tipo de distribucion de knots, u, e, m

lum = 0 ! # de intervalors u en m
c = 0.d0 ! parametod en m dist u en [a,c] y e en (c,b]

gamma = 0.d0 ! parametro de la exponencial

kord = 5 ! orden de los B-splines

l = 20 ! # de intervalos

intg = 100 ! grado de la intregracion de cuadratura gaussiana, >=k

nev = 10 ! # de autovalores que queremos calculara

etai = 0.1d0; etaf = 1.d0 

num_puntos_eta = 10

me = 1.d0

if( intg<kord ) intg = kord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nk = l+2*kord-1   ! # de knots
! No incluimos primer   y ultimo spline debido a psi(0)=psi(R)=0
nb = kord+l-3      ! size base splines 
if( nev<0 ) nev = nb

l1 = char(modulo(l, 10) + 48);
l2 = char(modulo(l, 100)/10 + 48);

!###########################################################
!###########################################################
!###########################################################
!open(9,file='electron-autovalores-L='//l2//l1//'.dat')
open(9,file='autovalores-prueba-2.dat')
!###########################################################
!###########################################################

write(9,*) '# El codigo usa B-splines para calcular el hamiltoniano y'
write(9,*) '# calcula los autovalores y autovectores que elija.'

write(9,*) '# El intervalo de integracion es [0, Rmax].'

write(9,*) '# El codigo va variando la variable eta, devuelve los autovalores como'
write(9,*) '# funcion de eta. Eta representa la intesidad del potencial Vef'

write(9,*) '# tip: entrar el tipo de distribucion, u, e, m'
write(9,*) '# tip =', tip

write(9,*) '# lum: entrar el # de intervalos u en m'
write(9,*) '# lum =', lum

write(9,*) '# c: parametro en m dist u en [a,c] y e en (c,b]'
write(9,*) '# c =', c

write(9,*) '# gamma: parametro en exponencial'
write(9,*) '# gamma =', gamma

write(9,*) '# kord: orden de los b-splines'
write(9,*) '# kord =', kord

write(9,*) '# l: # de intervalos'
write(9,*) '# l =', l

write(9,*) '# nb: tamaño base'
write(9,*) '# nb =', nb

write(9,*) '# intg: grado de integracion de la cuadratura >= k'
write(9,*) '# intg =', intg

write(9,*) '# nev: # de autovalores a calcular <=nb'
write(9,*) '# nev =', nev

write(9,*) '# rango de valores de eta='
write(9,*) '# etai =', etai, 'etaf =', etaf

write(9,*) '# masa de electron me =', me

write(9,*) '# autovalores calculados'
!###########################################################
!###########################################################
!###########################################################
close(9)

allocate( Sp(kord), dsp(kord-1))

allocate( x(l,intg), w(l,intg), pl(l,intg))

allocate( t(nk), k(nk))

allocate( norma(nb), s(nb,nb), v01(nb,nb), ke(nb,nb))

allocate( Vef(nb, nb, nb, nb))

t1 = omp_get_wtime();
call  KNOTS_PESOS( kord, tip, gamma, Rmin, Rmax, c, l, lum, intg, t, k, x, w, pl)

t2 = omp_get_wtime();
!do i = 1, nk
!!  write(*,*)i, t(i), k(i)
!end do

!!write(*,*)"calcula las matrices de una particulas"
t3 = omp_get_wtime();
call matrix_elements( )

t4 = omp_get_wtime();
!!write(*,*)"calcula la interaccion"
call interaccion( )

t5 = omp_get_wtime();

do i = 1,nb
  do j = i+1,nb
    s(j,i) = s(i,j)
    v01(j,i) = v01(i,j)
    ke(j,i) = ke(i,j)
  end do
end do

!!write(*,*)"diagonaliza"
call init_e( )
t6 = omp_get_wtime();

!write(*,*) " tiempo de KNOTS_PESOS =", t2-t1
!write(*,*) " tiempo de calculo_matrices =", t3-t2
!write(*,*) " tiempo de calculo_interaccion =", t4-t3
!write(*,*) " tiempo de calculo_autovalores =", t5-t4
!write(*,*) " tiempo total =", t5-t1

deallocate(Sp, dsp, x, w, pl)
deallocate(t, k, norma, s, v01, ke, Vef)

call cpu_time(time)
!write(*,*)time/60.d0

end !termina el main, termina el programa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_elements( )
use carga
use matrices
use integracion
implicit none
integer :: i, j, m, n, im, in
real(8) :: bm, bn, rr

s = 0.d0 ; v01 = 0.d0 ; ke = 0.d0

do i = kord, kord+l-1
  do j = 1, intg
    rr = x(k(i),j)
    call bsplvb(t, kord, 1, rr, i, Sp)
    
    do m = 1,kord
      im = i-kord+m-1
      if( im>0.and.im<nb+1 ) then
        do n = 1,kord
          in = i-kord+n-1
          
          if( in>0.and.in<nb+1 )then
            s(im,in) = s(im,in) + sp(m)*sp(n)*w(k(i),j);
                
            v01(im, in) = v01(im, in) + w(k(i),j)*sp(m)*sp(n)/rr;
            
          endif
        end do
      endif
    end do
  end do
end do

do i = kord, kord+l-1
  do m = i-kord+1, i
    if(m>1.and.m<nb+2)then
      do n = m,i
        if(n<nb+2)then
          do j = 1, intg
          
            call bder(x(k(i),j), t, bm, bn, kord, nk, m, n, i)
            
            ke(m-1,n-1) = ke(m-1,n-1) + 0.5d0*w(k(i),j)*bm*bn/me
                    
          end do
        endif
      end do
    endif
  end do
end do


end subroutine matrix_elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interaccion( )
use carga
use matrices
use integracion
implicit none
  integer :: i1, i2, j1, j2, i, j,it
  integer :: m, n, im, in, mp, imp, np
  real(8) :: rr1, rr2, w1, w2
  real(8), allocatable, dimension(:,:) :: f, g

  allocate(f(nb, nb), g(nb, nb))

  Vef = 0.d0;

  do i1 = kord, kord+l-1
    do j1 = 1, intg
      rr1 = x(k(i1), j1); w1 = w(k(i1), j1);

      f = 0.d0; g = 0.d0;
      do i2 = kord, kord+l-1
        do j2 = 1, intg
          rr2 = x(k(i2), j2); w2 = w(k(i2),j2);

          call bsplvb(t, kord, 1, rr2, i2, Sp)
          

          do m = 1, kord
            im = i2-kord+m-1
            if(im>0 .and. im<nb+1)then
              do n = 1, kord
                in = i2-kord+n-1

                if(in>0 .and. in<nb+1)then

                  if(rr2.LE.rr1)then

                    f(im, in) = f(im, in) + Sp(m)*Sp(n)*w2/rr1;

                  else

                    g(im, in) = g(im, in) + Sp(m)*Sp(n)*w2/rr2;

                  end if

                end if
              end do
            end if
          end do
        end do
      end do

      Sp = 0.d0;
      call bsplvb(t, kord, 1, rr1, i1, Sp)
      do m = 1, kord
        im = i1-kord+m-1
        if(im>0 .and. im<nb+1) then
          do mp = 1, kord
            imp = i1-kord+mp-1;
            if(imp>0 .and. imp<nb+1) then
              do n = 1, nb
                do np = 1, nb

                  Vef(im, n, imp, np) = Vef(im, n, imp, np) + Sp(m)*Sp(mp)*w1*(f(n, np) + g(n, np))/dsqrt(s(n,n)*s(np,np));

                end do
              end do
            end if
          end do
        end if
      end do
    end do
  end do
  write(*,*) Vef
  stop

  deallocate(f, g)

end subroutine interaccion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_e( )
use matrices
use integracion
use carga
implicit none
integer(4) :: m, n, np, mp
real(8) :: nor, nor1, nor2

do n = 1, nb
  do np = 1, nb
    nor1 = dsqrt(s(n,n)*s(np,np))
    do m = 1, nb
      do mp = 1, nb
 
!        nor2 = dsqrt(s(m,m)*s(mp,mp));
        Vef(n,m,np,mp) = Vef(n,m,np,mp)/(nor1)!*nor2)
      end do
    end do
  end do
end do

do m = 1, nb
   do n = m+1, nb
      nor = dsqrt(s(m,m)*s(n,n))
      s(m,n) = s(m,n)/nor ; s(n,m) = s(m,n)
      ke(m,n) = ke(m,n)/nor ; ke(n,m) = ke(m,n)
      v01(m,n) = v01(m,n)/nor ; v01(n,m) = v01(m,n)
   end do

  ke(m,m) = ke(m,m)/s(m,m)
  v01(m,m) = v01(m,m)/s(m,m)
end do
 
do m = 1, nb
   norma(m) = dsqrt(s(m,m))
   s(m,m) = 1.d0
end do

call sener( )

return

end subroutine init_e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! sener.f90 basado en ener.f90 
! pero es una subroutine de exp_bs_gs.f90 para calcular el GS

! Calcula E_i de una matriz que se le pasa por mudule

subroutine sener( )
use carga
use matrices
use integracion
implicit none
integer(4) :: i, j, m, n, mp, np, dp
integer(4) :: info, ind, indp, NN
real(8) :: raiz
real(8), allocatable, dimension(:,:) :: auvec, val_exp
real(8), allocatable, dimension(:) :: v, e
real(8), allocatable, dimension(:,:) :: mh, ms, mv, hsim


raiz = 1.d0/sqrt(2.d0);

dp = nb
NN = nb*(nb+1)/2

allocate( e(nev), v(nev), auvec( NN,nev), mh(dp,dp), val_exp(nev, nev))
allocate( hsim( NN, NN), ms(NN,NN), mv(NN,NN))

!###########################################################
!###########################################################
!###########################################################
!open(31,file='electron-autovalores-L='//l2//l1//'.dat',position='append')
open(11,file='autovalores-prueba-2.dat',position='append')
open(12,file='expectacion_interaccion.dat')
!###########################################################
!###########################################################
!###########################################################

mh = 0.d0

!!! hamiltoniano de una particula
do m = 1, dp
  do n = 1, m
    mh(m,n) = ke(m,n) - v01(m,n) ; 
    mh(n,m) = mh(m,n)
  end do
end do

delta = (etaf-etai)/dble(num_puntos_eta)

do i = 0, num_puntos_eta-1
!i = 0
  eta = etai + dble(i)*delta;
!!  write(*,*) i, eta
  hsim = 0.d0; ms = 0.d0;
 
  !!! halmiltoniano de dos particulas ya simetrizado 
  ind = 1
  do n = 1, dp
    do m = n, dp
      indp = 1
      do np = 1, dp
        do mp = np, dp

          if(m.eq.n .and. mp.eq.np)then

            hsim( ind, indp) = 2.d0*s(n,np)*mh(n,np) + dble(eta)*Vef(n,n,np,np);
            ms( ind, indp) = s(n,np)*s(n,np);
            mv( ind, indp) = Vef(n,n,np,np);

          elseif(m.ne.n .and. mp.eq.np )then

            hsim( ind, indp) = raiz*( 2.d0*s(m,np)*mh(n,np) + 2.d0*s(n,np)*mh(m,np) &
                               & + dble(eta)*Vef(n,m,np,np) + dble(eta)*Vef(m,n,np,np) );
            ms( ind, indp) = 2.d0*raiz*s(n,np)*s(m,np);
            mv( ind, indp) = 2.d0*raiz*(Vef(n,m,np,np) + Vef(m,n,np,np));

          elseif(m.eq.n .and. mp.ne.np)then

            hsim( ind, indp) = raiz*( 2.d0*s(n,mp)*mh(n,np) + 2.d0*s(n,np)*mh(n,mp) &
                               & + dble(eta)*Vef(n,n,np,mp) + dble(eta)*Vef(n,n,mp,np) );
            ms( ind, indp) = 2.d0*raiz*s(n,np)*s(n,mp);
            mv( ind, indp) = 2.d0*raiz*(Vef(n,n,np,mp) + Vef(n,n,mp,np));

          else

            hsim( ind, indp) = s(n,np)*mh(m,mp) + s(n,mp)*mh(m,np) &
                            &+ s(m,mp)*mh(n,np) + s(m,np)*mh(n,mp) &
                            &+ dble(eta)*0.5d0*(Vef(n,m,np,mp) + Vef(n,m,mp,np) + Vef(m,n,np,mp) + Vef(m,n,mp,np));
            ms( ind, indp) = s(n,np)*s(m,mp) + s(n,mp)*s(m,np);
            mv( ind, indp) = 0.5d0*(Vef(n,m,np,mp) + Vef(n,m,mp,np) + Vef(m,n,np,mp) + Vef(m,n,mp,np));

          endif
          indp = indp + 1
        end do
      end do
      ind = ind + 1
    end do
  end do

  if( NN.ne.(ind-1) ) stop
  if( NN.ne.(indp-1) ) stop  
  !!write(*,*) NN, ind-1, indp-1
  !do n = 1, NN
  !do m = 1, NN
  !!  write(*,*) n,m,hsim(n,m)!-hsim(m,n)
  !end do
  !end do

  call eigen_value( NN, nev, info, hsim, ms, e, auvec)

  val_exp = MATMUL( TRANSPOSE(auvec), MATMUL(mv, auvec)) 
  do j = 1, nev
    v(j) = val_exp(j,j);
  end do
  
  !muestro los valores
  write(11,6) eta, (e(m), m = 1, 25)
  write(12,6) eta, (v(m), m = 1, 25)
  write(*,6) eta, (e(m), m = 1, 2)
  
  
  e = 0.d0
  info = 0

end do
close(11)
close(12)

deallocate(e, v, auvec, mh, ms, hsim, mv, val_exp)

6 format(e22.14,1x,1000(1x,e22.14))

return
end subroutine sener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eigen_value(dp,nev,info,nh,s,x,avec)

!
! calcula usando LAPACK N=nev autovalores

implicit none
integer*4 dp,nev
real(8),dimension (nev)::x
real(8),dimension (dp,nev)::avec
real(8),dimension (dp,dp)::nh,s
! definiciones para DSYGVX
integer(4)::itype,numau,lwork,info,il,iu
character(1)::jobz,range,uplo
real(8)::vl,vu,abstol
real(8),dimension(dp,dp)::z
real(8),dimension(9*dp)::work
integer(4),dimension(5*dp)::iwork
integer(4),dimension(dp)::ifail

lwork=9*dp ! lwork=dim(work)

!!!!!!!!!!!!!!!!!!!!!!!!!
! reemplaza  tred2 y tqli del NR por dsygvx de LAPACK
! parametros para dsygvx:
itype=1  ! especifica que A x=lambda Bx
jobz='V' ! V (N) (NO) calcula autovectores
!range='V'! autovalores en el intervalo [VL,VU]
!vl=-1.01d0 ! pongo  vl menor que el menor autovalor
!vu=0.d0   ! por ahora solo interezan los autovalores negativos
uplo='U' ! la matriz es storeada superior
! N es dp
!A es nh
! lda = dp
! B = s  ! NO hay que hacer cholesky !!!
! ldb = dp
! il, iu  autov entre il y iu SI range='I'
 range='I'
il=1
iu=nev 
abstol=1.d-12  ! por decir algo (?)aconsejan abstol=2*DLAMCH('S') ; ????
! M = numau# de autovalores calculados : output
! W = x , de menor a mayor
! se debe definir un array Z(LDZ,M) NO conocemos M a priori! poner N=dp
! Z devuelve eigenvectors
! ldz=dp : dimension de Z
! definir work(lwork)
lwork=9*dp
! definir integer array iwork(5*N)
! ifail=output:
!If JOBZ = 'V', then if INFO = 0, the first M elements of
!          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!          indices of the eigenvectors that failed to converge.
!          If JOBZ = 'N', then IFAIL is not referenced.
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  DPOTRF or DSYEVX returned an error code:
!             <= N:  if INFO = i, DSYEVX failed to converge;
!                    i eigenvectors failed to converge.  Their indices
!                    are stored in array IFAIL.
!             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
!                    minor of order i of B is not positive definite.
!                    The factorization of B could not be completed and
!                    no eigenvalues or eigenvectors were computed.


call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, dp, nh, dp, s, dp,            &
     &                   VL, VU, IL, IU, ABSTOL, numau, x, Z, dp, WORK,  &
     &                   LWORK, IWORK, IFAIL, INFO )


!!!!!!!!!!!!!!!!!!

avec(1:dp,1:nev)=z(1:dp,1:nev)



        return
        end subroutine eigen_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  bajada de http://www.atom.physto.se/~lindroth/comp08/bget.f


       subroutine bder(rr,t,dm,dn,kord,np,indexm,indexn,left)
!     returns the value of (d/dx) Bspline(kord,index) in rr
!     The first Bspline is called spline #1 (i.e. index=1)
!     the first knot point is in t(1)
!     np= number of knot points (distinct or multiple) including
!     the ghost points: N phyical points np=N +2*(kord-1)


      implicit none
      integer i,left,kord,indexm,indexn,np
      real*8 dm,dn,t(np),Sp(kord),rr

      dm=0.d0  ; dn=0.d0
!     if rr=t(np) then the routine assumes that
!     there is kord knotpoints in the last physical point and
!     returns Bder.ne.zero if index is np-kord, or np-kord-1
      if(rr.gt.t(np)) return
      if(rr.lt.t(1)) return
!      do it=1,np
!        if(rr.ge.t(it)) left=it
!      end do

      if(abs(rr-t(np)).lt.1.d-10) then
        !if(index.lt.np-kord-1) return

        if(indexm.eq.np-kord) then
          dm=dble(kord-1)/(t(np)-t(np-kord))
        else if(indexm.eq.np-kord-1) then
          dm=-dble(kord-1)/(t(np)-t(np-kord))
        end if

       if(indexn.eq.np-kord) then
          dn=dble(kord-1)/(t(np)-t(np-kord))
        else if(indexn.eq.np-kord-1) then
          dn=-dble(kord-1)/(t(np)-t(np-kord))
        end if

        return
      end if

      call bsplvb(t,kord-1,1,rr,left,Sp)

  if(indexm-left+kord.ge.1.or.indexm-left+kord.le.kord)then
      i=indexm-left+kord
      if(i.eq.1) then
        dm=dble(kord-1)*(-Sp(i)/(t(indexm+kord)-t(indexm+1)))
      else if(i.eq.kord) then
        dm=dble(kord-1)*(Sp(i-1)/(t(indexm+kord-1)-t(indexm)))
      else
        dm=dble(kord-1)*(Sp(i-1)/(t(indexm+kord-1)-t(indexm))- &
     &  Sp(i)/(t(indexm+kord)-t(indexm+1)))
      end if
  endif
  if(indexn-left+kord.ge.1.or.indexn-left+kord.le.kord)then
      i=indexn-left+kord

      if(i.eq.1) then
        dn=dble(kord-1)*(-Sp(i)/(t(indexn+kord)-t(indexn+1)))
      else if(i.eq.kord) then
        dn=dble(kord-1)*(Sp(i-1)/(t(indexn+kord-1)-t(indexn)))
      else
        dn=dble(kord-1)*(Sp(i-1)/(t(indexn+kord-1)-t(indexn))- &
     &  Sp(i)/(t(indexn+kord)-t(indexn+1)))
      end if

  endif
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE bsplvb(t,jhigh,index,x,left,biatx)
!      INCLUDE 'standa.typ'
!      INCLUDE 'spldim.def'
!      include 'ndim_inc'
       implicit none
      integer(4),PARAMETER::JMAX=100
      integer index,jhigh,left,i,j,jp1
      real*8 t,x,biatx,deltal,deltar,saved,term
      DIMENSION biatx(jhigh),t(left+jhigh),deltal(jmax),deltar(jmax)
      SAVE deltal,deltar
      DATA j/1/
!      write(6,*) ' jmax=',jmax
      GO TO (10,20),index
 10   j = 1
      biatx(1) = 1.d0
      IF (j .GE. jhigh) GO TO 99

 20   CONTINUE
         jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
!         write(6,'(1pd12.4,2(i5,1pd14.6))')
!     :   x,left+j,t(left+j),left+1-j,t(left+1-j)
!         write(6,'(i3,1p3d12.4)') j,deltal(j),deltar(j),
!     :   abs(deltal(j)-deltar(j))
         saved = 0.d0
         DO i = 1,j
!         write(6,'(2i3,1p3d12.4)') i,j,deltal(jp1-1),deltar(i),
!     :   abs(deltal(jp1-1)-deltar(i))

             term = biatx(i)/(deltar(i) + deltal(jp1-i))
             biatx(i) = saved + deltar(i)*term
             saved = deltal(jp1-i)*term
         END DO
         biatx(jp1) = saved
         j = jp1
         IF (j .LT. jhigh) GO TO 20
 99   RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE KNOTS_PESOS(kord,tip,gamma,a,b,c,l,lum,intg,t,k,x,w,pl)

!lo único que hace esta subrutina principal es dividir los casos segun se quiera knots uniformes, exponencial, mixto
!cada una de las subrutinas que llama ademas calcula los x y los w, abscisas y pesos para la cuadratura gaussiana, y los polinomios


!  INPUT:
! kord: orden de los b-splines
! tip; character(1): 'u' ; 'e' ; 'm': dist uniforme, exp o mixta de knots
! gamma : param dist. e (no usado si tip='u')
! a
! b  ; a<b todo es calculado en el intervalo [a,b]
! c solo usado si tip='m' => a<c<b ; dist u  en [a,c]; e en (c,b]; c es a u-knot
! NOTA I : dist optima: m; con c~r_0 donde v(r) es apreciablemente no nulo
! l = # total de sub_intervalos en [a,b]
! lum: solo se usa si tip=m : # subint con dist u => l-lem con dist e
! intg = # de puntos que usamos en la integral x cuadratura en c/subintervalo

 
!  OUTPUT:
!  t(l+2*kord-1)  ! los (l+2*kord-1) knots (contando multiplicidades en a y b)
!  k(l+2*kord-1)  ! da el # de intervalo j para el i-esimo nodo 1<=k<=l
!  x(l,intg),w(l,intg) : posiciones y pesos de los intg puntos en c/intervalo
! pl(l,intg),w(l,intg) : polinomios de Legendre en cada punto
! NOTA II : los p_l los calcula gratis la subr gauleg de int. x cuad. del NR 
! y a veces hacen falta=> version modif. de gauleg que los da como output

implicit none
integer(4)::kord,l,lum,intg,i1,j1
character(1)::tip
real(8)::gamma,a,b,c
integer(4),dimension(l+2*kord-1)::k
real(8),dimension(l+2*kord-1)::t
real(8),dimension(l,intg)::x,w,pl
real(8)::rr2
!!!!!!


if(tip=='u')then
call dist_unif(kord,a,b,l,intg,t,k,x,w,pl)
elseif(tip=='e')then 
call dist_exp(kord,gamma,a,b,l,intg,t,k,x,w,pl)
elseif(tip=='m')then 
call dist_mix(kord,gamma,a,b,c,l,lum,intg,t,k,x,w,pl)
else
!write(*,*)'error 1 en KNOTS_PESOS :',tip,' no corresponde a una distribucion'
stop
endif

return 
end    SUBROUTINE KNOTS_PESOS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DIST_UNIF(kord,a,b,l,intg,t,k,x,w,pl)

implicit none
integer(4)::kord,l,intg
real(8)::a,b
integer(4),dimension(l+2*kord-1)::k
real(8),dimension(l+2*kord-1)::t
real(8),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,nk
real(8)::ri,rf,dr
real(8),dimension(intg)::vx,vw,vpl

nk=l+2*kord-1   ! # de knots

dr=(b-a)/dfloat(l)

! calcula los puntos y pesos para la cuadratura
x=0.d0;w=0.d0

do i=1,l
   
   ri=a+dfloat(i-1)*dr
   rf=ri+dr

   call gauleg_pl(ri,rf,vx,vw,vpl,intg)

   do j=1,intg
   x(i,j)=vx(j)
   w(i,j)=vw(j)
   pl(i,j)=vpl(j)
   end do

end do

t(1)=a
k(1)=1

do i=2,kord
   t(i)=t(1)
   k(i)=1
end do

do i=kord+1,kord+l
   t(i)=t(i-1)+dr
   k(i)=k(i-1)+1
end do

do i=kord+l+1,nk
   t(i)=t(i-1)
   k(i)=k(i-1)
end do

return

end    SUBROUTINE DIST_UNIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIST_EXP(kord,gamma,a,b,l,intg,t,k,x,w,pl)

implicit none
integer(4)::kord,l,intg
real(8)::gamma,a,b
integer(4),dimension(l+2*kord-1)::k
real(8),dimension(l+2*kord-1)::t
real(8),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,nk
real(8)::ri,rf,dr,ye
real(8),dimension(intg)::vx,vw,vpl

nk=l+2*kord-1   ! # de knots

dr=(b-a)/(dexp(gamma)-1.d0)
ye=gamma/dfloat(l)

! calcula los puntos y pesos para la cuadratura
x=0.d0;w=0.d0

do i=1,l

   ri=a+dr*(dexp(ye*dfloat(i-1))-1.d0)
   rf=a+dr*(dexp(ye*dfloat(i))-1.d0)

   call gauleg_pl(ri,rf,vx,vw,vpl,intg)

   do j=1,intg
      x(i,j)=vx(j)
      w(i,j)=vw(j)
      pl(i,j)=vpl(j)
   end do

end do

t(1)=a
k(1)=1

do i=2,kord
   t(i)=t(1)
   k(i)=1
end do

do i=kord+1,kord+l
   t(i)=a+dr*(dexp(ye*dfloat(k(i-1)))-1.d0)
   k(i)=k(i-1)+1
end do

do i=kord+l+1,nk
t(i)=t(i-1)
k(i)=k(i-1)
end do

return
end    SUBROUTINE DIST_EXP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DIST_MIX(kord,gamma,a,b,c,l,lum,intg,t,k,x,w,pl)
implicit none
integer(4)::kord,l,lum,intg
real(8)::gamma,a,b,c
integer(4),dimension(l+2*kord-1)::k
real(8),dimension(l+2*kord-1)::t
real(8),dimension(l,intg)::x,w,pl
!!!!!!
integer(4)::i,j,le,nk
real(8),dimension(lum,intg)::xu,wu,plu
real(8),dimension(l-lum,intg)::xe,we,ple
integer(4),dimension(lum+2*kord-1)::ku
real(8),dimension(lum+2*kord-1)::tu
integer(4),dimension(l-lum+2*kord-1)::ke
real(8),dimension(l-lum+2*kord-1)::te

nk=l+2*kord-1   ! # de knots

call DIST_UNIF(kord,a,c,lum,intg,tu,ku,xu,wu,plu)

le=l-lum

call DIST_EXP(kord,gamma,c,b,le,intg,te,ke,xe,we,ple)

do i=1,lum
x(i,1:intg)=xu(i,1:intg)
w(i,1:intg)=wu(i,1:intg)
pl(i,1:intg)=plu(i,1:intg)
end do

do i=1,le
x(lum+i,1:intg)=xe(i,1:intg)
w(lum+i,1:intg)=we(i,1:intg)
pl(lum+i,1:intg)=ple(i,1:intg)
end do

do i=1,kord+lum
t(i)=tu(i)
k(i)=ku(i)
end do

j=kord+lum
do i=kord+1,le+kord+1
j=j+1
t(j)=te(i)
k(j)=ke(i)+ku(kord+lum)-1
end do

do i=j,nk
t(i)=t(i-1)
k(i)=k(i-1)
end do

return
end    SUBROUTINE DIST_MIX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! subrutina de NR modificada
! calcula coordenadas y pesos GL y los valores de P_l(x)

      SUBROUTINE gauleg_pl(x1,x2,x,w,pl,n)
      INTEGER n
      DOUBLE PRECISION x1,x2,x(n),w(n),pl(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=dcos(3.1415926535897932385d0*(dfloat(i)-.25d0)/(dfloat(n)+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/dfloat(j)
11        continue
          pl(i)=p1
          pp=dfloat(n)*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(dabs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

