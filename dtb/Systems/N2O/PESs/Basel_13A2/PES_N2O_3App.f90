  !   atom 1: N
  !   atom 2: O
  !   atom 3: N
  !
  ! The PES for the NO+N in the ^3A" state is call 
  !
  !     call  call initialize_PES()
  !     call USERE(EU,X,Y,Z,DX,DY,DZ,QECONT,ECONT,NATOMX)
  !
  !     EU   - User energy to be returned
  !     X,Y,Z - Coordinates for energy evaluation
  !     DX,DY,DZ - Forces. Must be modified. Append (dE/dX,...)
  !     QECONT - Flag for analysis (0-no analysis,>0=analysis)
  !     ECONT(NATOMX) - Analysis array to be filled if QECONT>0.
  !     NATOMX - Number of atoms = 3
  !
  !     
module chm_kinds
    implicit none
    INTEGER, PARAMETER :: chm_real = 8   
    
end module    
!***** For Charmm just remove until here and the comments with 3 exclamations (!!!) ***********
!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains the 1-dimensional kernel functions and their derivatives
!
!///////////////////////////////////////////////////////////////////////////////
module RKHS
    implicit none
    private
    !declare all methods/variables available from the outside
    public :: q_ker, dq_ker, k_ker, dk_ker               
    
!functions and subroutines start here
contains
    !the reproducing kernel q_ker for distance like variables
    real*8 function q_ker(x,r)
        implicit none
        real*8 :: x,r,x_l,x_s
        x_l=r
           x_s=x
        if(r < x) then
            x_s=r
            x_l=x
        endif
        q_ker = (1d0/(14d0*x_l**7))*(1d0-7*x_s/(9*x_l))
    end function q_ker
    
    !derivative of q_ker
    real*8 function dq_ker(x,r)
        implicit none
        real*8 :: x,r
        if(x < r) then
            dq_ker = -1d0/(18d0*r**8)
        else
            dq_ker = -(1d0/18d0)*(9*x-8*r)/(x**9)
        end if
    end function dq_ker


    !the reproducing kernel k_ker for angle like variables
    real*8 function k_ker(x,r)      !k=2 Anglelike
        implicit none
        real*8 :: x,r,x_l,x_s
        x_l=r !Large value of the coordinate
        x_s=x !Small value of the coordinate
        if(r < x) then !Here the smaller and larger value are chosen
            x_s=r
            x_l=x
        endif
        !Only two terms (up to n=1) were taken of the Gauss´ hypergeometric function expansion  
        k_ker=1+x_s*x_l+2*x_s**2*x_l*(1.d0-x_s/(3*x_l))
    end function k_ker

    !derivative of k_ker
    real*8 function dk_ker(x,r)
        implicit none
        real*8 :: x,r
        if(x < r) then
            dk_ker = r+4*x*r*(1.d0-x/(3*r))-(2*x**2)/3d0
        else
            dk_ker = r+2*r**2*(1.d0-r/(3*x))+(2*r**3)/(3*x)
        end if
    end function dk_ker
end module RKHS

!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains routines to read files and store them in arrays
!
!///////////////////////////////////////////////////////////////////////////////
module fileInOut
    implicit none
    private
    public :: readData, readData2, writeData, writeData2

    contains
    !subroutine to read in files
    subroutine readData(unitnr,ios,array, inputfile,dim1)
        implicit none
        !declaration of variables
        integer :: unitnr,dim1
        integer :: i,ios
        real*8 :: array(dim1)
        character(len=*) :: inputfile

        open(unitnr,file=inputfile, status='old', action='read', &
                 iostat=ios)
        !error check
        if (ios /= 0) then
            return
        end if 
        do i=1,dim1
            read(unitnr,*) array(i)
        end do
        close(unit=unitnr)    
    end subroutine readData

    subroutine readData2(unitnr,ios,array, inputfile,dim1,dim2)
        implicit none
        !declaration of variables
        integer :: unitnr,dim1,dim2
        integer :: i,ios
        real*8 :: array(dim1,dim2)
        character(len=*) :: inputfile

        open(unitnr,file=inputfile, status='old', action='read', &
                 iostat=ios)
        !error check
        if (ios /= 0) then
            return
        end if 
        do i=1,dim1
            read(unitnr,*) array(i,:)
        end do
        close(unit=unitnr)    
    end subroutine readData2

    subroutine writeData(unitnr,ios,array,outputfile,dim1)
        implicit none
        integer :: unitnr,dim1
        integer :: i,ios
        real*8 :: array(dim1)
        character(len=*) :: outputfile
        
        open(unitnr, file=outputfile, status='new', action='write', &
                iostat=ios)
        if (ios /= 0) then
            return
        end if
        do i=1,dim1
            write(unitnr,*) array(i)
        end do
        close(unit=unitnr)
    end subroutine writeData

    subroutine writeData2(unitnr,ios,array,outputfile,dim1,dim2)
        implicit none
        integer :: unitnr,dim1,dim2
        integer :: i,ios
        real*8 :: array(dim1,dim2)
        character(len=*) :: outputfile
        
        open(unitnr, file=outputfile, status='new', action='write', &
                iostat=ios)
        if (ios /= 0) then
            return
        end if
        do i=1,dim1
            write(unitnr,*) array(i,:)
        end do
        close(unit=unitnr)
    end subroutine writeData2

end module fileInOut

!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains all the routines necessary to return the custom energies
!    and forces from the RKHS interpolation
!
!///////////////////////////////////////////////////////////////////////////////
module CustomEnergy
    use RKHS
    implicit none
    private
    public :: initialize_PES, Epot, dEdR, dEdVdW, dEdalpha, adEdR, adEdVdW, adEdalpha
    
    
!-------------------------------------CHANGE HERE-------------------------------------      
! The dimensions of the PES grid is specified here (number of points for the three coordinates)
    integer, parameter ::   nAlpha = 11  ,& !number of gridpoints for alpha
                            nVdW   = 18  ,& !25,&! !number of gridpoints for bigR
                            nR     = 12     !number of gridpoints for smallR
!-------------------------------------------------------------------------------------  
                            
    real*8, parameter :: PI = 3.141592653589793238462643383279502884197d0
   
    !Array to store the grid data
    real*8 :: dataArray(nAlpha*nVdW*nR,4)  !array that the ab initio points are read into
    real*8 :: r_grid(nR)                   !store r_grid for asymptote correction
    real*8, save :: dataArray1(nAlpha*nVdW*nR,4), dataArray2(nAlpha*nVdW*nR,4), dataArray3(nAlpha*nVdW*nR,4)  
    real*8, save :: r_grid1(nR), r_grid2(nR), r_grid3(nR)                         
      
    !Arrays for the weights (coefficients)
    real*8 :: w_tot(nAlpha*nR*nVdW), w_R(nR)
    real*8, save :: w_tot1(nAlpha*nR*nVdW), w_R1(nR)
    real*8, save :: w_tot2(nAlpha*nR*nVdW), w_R2(nR)
    real*8, save :: w_tot3(nAlpha*nR*nVdW), w_R3(nR) 

contains
    !Reads in the weights for the RKHS interpolation. This is invoked at the startup of CHARMM
    subroutine initialize_PES()
        use fileInOut 
        implicit none
        integer i
        integer :: IOS, IOS1, IOS2 !used to check file status
        call readData2(20,IOS ,dataArray1,"data/data1.txt",nAlpha*nVdW*nR,4)
        call readData2(20,IOS1,dataArray2,"data/data2.txt",nAlpha*nVdW*nR,4)
        call readData2(20,IOS2,dataArray3,"data/data3.txt",nAlpha*nVdW*nR,4)
        if ((IOS /= 0).or.(IOS1 /= 0).or.(IOS2 /= 0)) then
            print*, "One of the data files that contains the grid in the data folder could not be read properly."
            print*, "program terminated."
            stop
        end if 
                
        !change alpha to new coordinate y=(1-cos(alpha))/2
        do i=1,nR*nVdW*nAlpha
            dataArray1(i,1)=y(dataArray1(i,1))
            dataArray2(i,1)=y(dataArray2(i,1))
            dataArray3(i,1)=y(dataArray3(i,1))
        end do
        !read in r_grid (needed for asymptote correction)
        do i=1,nR
            r_grid1(i)=dataArray1(1+(i-1)*nVdW,2)
            r_grid2(i)=dataArray2(1+(i-1)*nVdW,2)
            r_grid3(i)=dataArray3(1+(i-1)*nVdW,2)
        end do
        
        !read weight/coefficient data for surface 1 from files
        call readData(21,IOS1,w_R1,"data/w_r1.dat",nR)
        call readData(22,IOS2,w_tot1,"data/w_tot1.dat",nAlpha*nR*nVdW)    
        !Error check
        if ((IOS1 /= 0).or.(IOS2 /= 0).or.(IOS /= 0)) then
            print*, "One of the weight-files in the data folder for surface 1 could not be read properly."
            print*, "program terminated."
            stop
        end if    

        !read weight/coefficient data for surface 2 from files
        call readData(21,IOS1,w_R2,"data/w_r2.dat",nR)
        call readData(22,IOS2,w_tot2,"data/w_tot2.dat",nAlpha*nR*nVdW)
        !Error check
        if ((IOS1 /= 0).or.(IOS2 /= 0).or.(IOS /= 0)) then
            print*, "One of the weight-files in the data folder for surface 2 could not be read properly."
            print*, "program terminated."
            stop
        end if 
        
        !read weight/coefficient data for surface 3 from files
        call readData(21,IOS1,w_R3,"data/w_r3.dat",nR)
        call readData(22,IOS2,w_tot3,"data/w_tot3.dat",nAlpha*nR*nVdW)
        !Error check
        if ((IOS1 /= 0).or.(IOS2 /= 0).or.(IOS /= 0)) then
            print*, "One of the weight-files in the data folder for surface 3 could not be read properly."
            print*, "program terminated."
            stop
        end if   
        
        return
    end subroutine initialize_PES


    !transforms alpha to the new variable y=(1-cos(alpha))/2
    real*8 function y(alpha) 
        implicit none
        real*8,intent(in) :: alpha
        y = (1d0-DCOS(alpha*PI/180d0))/2d0
    end function y
    
    real*8 function dy(alpha) !used for calculating the analytical derivative with respect to alpha (chainrule)
        implicit none
        real*8,intent(in) :: alpha
        dy = (PI/360d0)*DSIN(alpha*PI/180d0)
    end function dy

    !calculates the new asymptote depending for r'
    real*8 function asymp(R,surf)
        implicit none
        real*8, intent(in) :: R
        !which surface shall be calculated
        integer, intent(in) :: surf
        integer :: counter
        asymp = 0.d0
        
        !read appropriate weights, depending on which surface is selected
        if      (surf == 1) then
            w_R = w_R1
            r_grid = r_grid1
        else if (surf == 2) then
            w_R = w_R2
            r_grid = r_grid2
        else if (surf == 3) then
            w_R = w_R3
            r_grid = r_grid3
        else
            print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
            stop
        end if 
        do counter=1,nR
            asymp = asymp + w_R(counter)*q_ker(R,r_grid(counter))
        end do        
    end function asymp
    
    real*8 function dasymp(R,surf)
        implicit none
        real*8, intent(in) :: R
        !which surface shall be calculated
        integer, intent(in) :: surf
        integer :: counter
        dasymp = 0.d0
        
                !read appropriate weights, depending on which surface is selected
        if      (surf == 1) then
            w_R = w_R1
            r_grid = r_grid1
        else if (surf == 2) then
            w_R = w_R2
            r_grid = r_grid2
        else if (surf == 3) then
            w_R = w_R3
            r_grid = r_grid3
        else
            print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
            stop
        end if 
        do counter=1,nR
            dasymp = dasymp + w_R(counter)*dq_ker(R,r_grid(counter))
        end do        
    end function dasymp

    !calculates the energy V(alpha',r',R') for an off-grid point
    real*8 function Epot(R,VdW,alpha,surf)
        implicit none
        !input coordinates
        real*8, intent(in) :: R,VdW,alpha
        real*8 :: y_val !coordinate transformation
        !which surface shall be calculated?
        integer, intent(in) :: surf
        integer :: counter
        y_val = y(alpha)    !transforms alpha to y
        Epot = 0d0
        
        !read appropriate weights, depending on which surface is selected
        if      (surf == 1) then
            w_tot = w_tot1
            dataArray = dataArray1
        else if (surf == 2) then
            w_tot = w_tot2
            dataArray = dataArray2
        else if (surf == 3) then
            w_tot = w_tot3
            dataArray = dataArray3
        else
            print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
            stop
        end if 
        
        
        do counter=1,nAlpha*nR*nVdW
            Epot = Epot + w_tot(counter)*k_ker(y_val,dataArray(counter,1)) &
                *q_ker(R,dataArray(counter,2))*q_ker(VdW,dataArray(counter,3))
        end do
        Epot = Epot + asymp(R,surf)
    end function Epot
                
    !calculates the partial derivative of V(alpha,r,R) with respect to VdW numerically
    real*8 function dEdVdW(R,VdW,alpha,surf)
        implicit none
        real*8, intent(in) :: alpha,R,VdW
        integer, intent(in) :: surf    
        real*8 :: eps
        real*8, dimension(4) :: stencil
        !construct epsilon
        !eps = dabs(dsqrt(epsilon(0.d0))*VdW)
        !if(eps < dsqrt(epsilon(0.d0))) eps = dsqrt(epsilon(0.d0))
        eps = 1d-5
        !build stencil
        stencil(1) = Epot(R, VdW +2*eps, alpha, surf)
        stencil(2) = Epot(R, VdW +1*eps, alpha, surf)
        stencil(3) = Epot(R, VdW -1*eps, alpha, surf)
        stencil(4) = Epot(R, VdW -2*eps, alpha, surf)
        dEdVdW = derivative_FivePointStencil(stencil,eps)   
    end function dEdVdW

    !calculates the partial derivative of V(alpha,r,R) with respect to r numerically
    real*8 function dEdR(R,VdW,alpha,surf)
        implicit none
        real*8, intent(in) :: alpha,R,VdW  
        integer, intent(in) :: surf 
        real*8 :: eps
        real*8, dimension(4) :: stencil
        !construct epsilon
        !eps = dabs(dsqrt(epsilon(0.d0))*R)
        !if(eps < dsqrt(epsilon(0.d0))) eps = dsqrt(epsilon(0.d0))
        eps = 1d-5
        !build stencil
        stencil(1) = Epot(R +2*eps, VdW, alpha, surf)
        stencil(2) = Epot(R +1*eps, VdW, alpha, surf)
        stencil(3) = Epot(R -1*eps, VdW, alpha, surf)
        stencil(4) = Epot(R -2*eps, VdW, alpha, surf)
        dEdR = derivative_FivePointStencil(stencil,eps)   
    end function dEdR

    !calculates the partial derivative of V(alpha,r,R) with respect to alpha numerically
    real*8 function dEdalpha(R,VdW,alpha,surf)
        implicit none
        real*8, intent(in) :: alpha,R,VdW
        integer, intent(in) :: surf 
        real*8 :: eps, tempy, yplusdy, tempalpha
        real*8, dimension(4) :: stencil
        !construct epsilon
        !eps = dabs(dsqrt(epsilon(0.d0))*alpha)
        !if(eps < dsqrt(epsilon(0.d0))) eps = dsqrt(epsilon(0.d0))
        eps = 1d-5        
        !This is a long block to construct an appropriate displacement for alpha
        tempy = y(alpha)
        yplusdy = tempy + eps
        if (yplusdy > 1.d0) yplusdy = tempy - eps 
        tempalpha = dacos(1-2*yplusdy)
        eps = dabs(alpha-tempalpha)
        !build stencil            
        stencil(1) = Epot(R, VdW, alpha +2*eps, surf)
        stencil(2) = Epot(R, VdW, alpha +1*eps, surf)
        stencil(3) = Epot(R, VdW, alpha -1*eps, surf)
        stencil(4) = Epot(R, VdW, alpha -2*eps, surf)
        dEdalpha = derivative_FivePointStencil(stencil,eps)
    end function dEdalpha
        
    real*8 function derivative_FivePointStencil(stencil,dx)
        implicit none
        real*8, dimension(4), intent(in) :: stencil
        real*8,               intent(in) :: dx
        derivative_FivePointStencil = (-stencil(1)+8*stencil(2)-8*stencil(3)+stencil(4))/(12*dx)
        return
    end function derivative_FivePointStencil
    
    !calculates the partial derivative of V(alpha,r,R) with respect to VdW analytically
    real*8 function adEdVdW(R,VdW,alpha,surf)
        implicit none
        real*8, intent(in) :: alpha,R,VdW
        real*8 :: y_val
        integer, intent(in) :: surf 
        integer :: counter
        y_val = y(alpha)    !transforms alpha to y
        adEdVdW = 0d0
        !read appropriate weights, depending on which surface is selected
        if      (surf == 1) then
            w_tot = w_tot1
            dataArray = dataArray1
        else if (surf == 2) then
            w_tot = w_tot2
            dataArray = dataArray2
        else if (surf == 3) then
            w_tot = w_tot3
            dataArray = dataArray3
        else
            print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
            stop
        end if 
        
        do counter=1,nAlpha*nR*nVdW
            adEdVdW = adEdVdW + w_tot(counter)*k_ker(y_val,dataArray(counter,1)) &
                    *q_ker(R,dataArray(counter,2))*dq_ker(VdW,dataArray(counter,3))
        end do
    end function adEdVdW

    !calculates the partial derivative of V(alpha,r,R) with respect to r analytically
    real*8 function adEdR(R,VdW,alpha,surf)
        implicit none
        real*8, intent(in) :: alpha,R,VdW
        real*8 :: y_val
        integer, intent(in) :: surf 
        integer :: counter
        y_val = y(alpha)    !transforms alpha to y
        adEdR = 0d0
        !read appropriate weights, depending on which surface is selected
        if      (surf == 1) then
            w_tot = w_tot1
            dataArray = dataArray1
        else if (surf == 2) then
            w_tot = w_tot2
            dataArray = dataArray2
        else if (surf == 3) then
            w_tot = w_tot3
            dataArray = dataArray3
        else
            print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
            stop
        end if 

        do counter=1,nAlpha*nR*nVdW
            adEdR = adEdR + w_tot(counter)*k_ker(y_val,dataArray(counter,1)) &
                        *dq_ker(R,dataArray(counter,2))*q_ker(VdW,dataArray(counter,3))
        end do
        adEdR = adEdR + dasymp(R,surf)
    end function adEdR

    !calculates the partial derivative of V(alpha,r,R) with respect to alpha analytically
    real*8 function adEdalpha(R,VdW,alpha,surf)
        implicit none
        real*8, intent(in) :: alpha,R,VdW
        real*8 :: y_val,dy_val
        integer, intent(in) :: surf 
        integer :: counter
        y_val = y(alpha)    !transforms alpha to y
        dy_val = dy(alpha)
        adEdalpha = 0d0
        
        !read appropriate weights, depending on which surface is selected
        if      (surf == 1) then
            w_tot = w_tot1
            dataArray = dataArray1
        else if (surf == 2) then
            w_tot = w_tot2
            dataArray = dataArray2
        else if (surf == 3) then
            w_tot = w_tot3
            dataArray = dataArray3
        else
            print*, "ERROR: Only surface numbers 1-3 are supported. Check call for Epot "
            stop
        end if 
        
        do counter=1,nAlpha*nR*nVdW
            adEdalpha = adEdalpha + w_tot(counter)*dk_ker(y_val,dataArray(counter,1))*dy_val &
                        *q_ker(R,dataArray(counter,2))*q_ker(VdW,dataArray(counter,3))
        end do
    end function adEdalpha

end module CustomEnergy

!///////////////////////////////////////////////////////////////////////////////
!
!    This module contains functions and subroutines needed for unit or coordinate
!    conversion or similar things. 
!
!///////////////////////////////////////////////////////////////////////////////
module utility_conversion
use chm_kinds
    implicit none
    private
    public :: calc_Jacobi_coordinates, calc_weights

contains
subroutine calc_Jacobi_coordinates(X,Y,Z,jacCoord,mass,dJacdCart)
!This subroutine calculates the Jacobi coordinates from the cartesian coordinates
!The indices of the Jacobi coordinate arrays coorespond to different surfaces:
!(the numbers here are the atom labels)
!index 1 (surface 1) : 1=2---3 
!index 2 (surface 2) : 1=3---2
!index 3 (surface 3) : 2=3---1
    implicit none
    
    !Input coordinates (cartesian)
    real(chm_real), dimension(3), intent(in) :: X, Y, Z
    
    !Input masses (needed for some centre of mass)
    real*8, dimension(3), intent(in) :: mass
    
    !Output coordinates (Jacobi)
    !1st index: surface
    !2nd index: Jacobi coordinate (r, bigR, alpha)
    real*8, dimension(3,3), intent(out) :: jacCoord
    
    !Output derivatives of Jacobi coordinates with respect to cartesian coordinates
    !1st index: surface
    !2nd index: Jacobi coordinate (r, bigR, alpha)
    !3rd index: cartesian coordinate (x1,y1,z1,x2,y2,z2,x3,y3,z3)
    real*8, intent(out) :: dJacdCart(3,3,9)
    
    !temporary array to store Jacobi coordinates. We need this for swapping back!
    real*8 :: dJdC_tmp(3,9)
    
    !Utility variable for convenience that stores all the coordinates in one array
    real*8, dimension(9) :: coord
    
    !temporary mass array (gets reordered depending on which surface is calculated)
    real*8, dimension(3) :: mass_tmp
    
    !loop variables
    integer i
    
    
    !calculate the Jacobi coordinates for surface 1 (1=2- -3)
        !store coordinates appropriately
    do i=1,3
        coord(1+3*(i-1)) = real(X(i),8)
        coord(2+3*(i-1)) = real(Y(i),8)
        coord(3+3*(i-1)) = real(Z(i),8)
    end do
    call cartesian2jacobi(coord,jacCoord(1,:),dJacdCart(1,:,:),mass)
    
    !calculate the Jacobi coordinates for surface 2 (1=3- -2)
        !store coordinates appropriately
    coord(1) = real(X(1),8)
    coord(2) = real(Y(1),8)
    coord(3) = real(Z(1),8)
    coord(4) = real(X(3),8)
    coord(5) = real(Y(3),8)
    coord(6) = real(Z(3),8)
    coord(7) = real(X(2),8)
    coord(8) = real(Y(2),8)
    coord(9) = real(Z(2),8)
    !build a temporary mass array with swapped around masses
    mass_tmp(1) = mass(1)
    mass_tmp(2) = mass(3)
    mass_tmp(3) = mass(2)
    call cartesian2jacobi(coord,jacCoord(2,:),dJacdCart(2,:,:),mass_tmp)
    
    !calculate the Jacobi coordinates for surface 3 (2=3- -1)
        !store coordinates appropriately
    coord(1) = real(X(2),8)
    coord(2) = real(Y(2),8)
    coord(3) = real(Z(2),8)
    coord(4) = real(X(3),8)
    coord(5) = real(Y(3),8)
    coord(6) = real(Z(3),8)
    coord(7) = real(X(1),8)
    coord(8) = real(Y(1),8)
    coord(9) = real(Z(1),8)
    !build a temporary mass array with swapped around masses
    mass_tmp(1) = mass(2)
    mass_tmp(2) = mass(3)
    mass_tmp(3) = mass(1)
    call cartesian2jacobi(coord,jacCoord(3,:),dJacdCart(3,:,:),mass_tmp) 
    
    !Swap back coordinates for surface 2
    dJdC_tmp = dJacdCart(2,:,:)
    dJacdCart(2,:,4:6) = dJdC_tmp(:,7:9)
    dJacdCart(2,:,7:9) = dJdC_tmp(:,4:6)
    
    !Swap back coordinates for surface 3
    dJdC_tmp = dJacdCart(3,:,:)
    dJacdCart(3,:,1:3) = dJdC_tmp(:,7:9)
    dJacdCart(3,:,4:6) = dJdC_tmp(:,1:3)
    dJacdCart(3,:,7:9) = dJdC_tmp(:,4:6)

    return
end subroutine calc_Jacobi_coordinates


subroutine cartesian2jacobi(x,xi,dxidx,m)
    !This subroutine expects a 9-vector as input that contains the cartesian coordinates,        
    !x -> cartesian coordinates (x1,y1,z1,x2,y2,z2,x3,y3,z3)

    !a 3-vector as output that contains the Jacobi coordinates,
    !xi -> 3 Jacobi coordinates (r,R,alpha)

    !an array to store the derivatives of the Jacobi coordinates with respect to the
    !cartesian coordinates. The 1st index corresponds to the Jacobi coordinate (r,R,alpha),
    !while the 2nd index corresponds to the cartesian coordinate (x1,y1,z1,x2,y2,z2,x3,y3,z3).
    !dxidx -> derivatives of the Jacobi coordinates with respect to cartesians

    !Also, a 3-vector containing the masses of the atoms is required. This is because the Jacobi
    !coordinates reference the centre of mass between atom1 and atom2.
    !m -> masses of the 3 atoms
              
    implicit none
    real*8, parameter   :: PI = 3.141592653589793238462643383279502884197d0
    real*8, intent(in)  :: x(9)
    real*8, intent(in)  :: m(3)
    real*8, intent(out) :: xi(3)
    real*8, intent(out) :: dxidx(3,9)
    real*8::  x12cm(3),dalphadx2(3),dalphadx12cm(3),dalphadx3(3)
                     
    !distance of atom1 and atom2 (Jacobi coordinate r)
    xi(1)  = distance(x(1:3),x(4:6))
            
    !centre of mass of atom1 and atom2 (needed to calculate bigR)
    x12cm  = (m(1)*x(1:3)+m(2)*x(4:6))/(m(1)+m(2))
                     
    !distance of atom3 and centre of mass of atom1 and atom2 (Jacobi coordinate bigR)
    xi(2)  = distance(x12cm,x(7:9))

    !Derivatives of the Jacobi coordinates with respect to the cartesian coordinates
    !dr/dx
    dxidx(1,1:3) =  (x(1:3)-x(4:6))/xi(1)
    dxidx(1,4:6) = -dxidx(1,1:3)
    dxidx(1,7:9) =  0.d0
    
    !dR/dx       
    dxidx(2,1:3) = (x12cm(1:3)-x(7:9))/xi(2)*m(1)/(m(1)+m(2)) 
    dxidx(2,4:6) = (x12cm(1:3)-x(7:9))/xi(2)*m(2)/(m(1)+m(2))   
    dxidx(2,7:9) = (x(7:9)-x12cm(1:3))/xi(2)     

    !calculate Jacobi coordinate alpha and derivatives needed to get the derivatives of
    !alpha with respect to cartesian coordinates (all in radians!!!)
    call angle(x(4:6),x12cm(1:3),x(7:9),xi(3),dalphadx2,dalphadx12cm,dalphadx3)
            
    !calculate the derivatives of alpha with respect to cartesian coordinates
    dxidx(3,1:3) = dalphadx12cm*m(1)/(m(1)+m(2))
    dxidx(3,4:6) = dalphadx12cm*m(2)/(m(1)+m(2)) + dalphadx2
    dxidx(3,7:9) = dalphadx3  
            
    !transform the derivatives of alpha and alpha itself from radians 
    !to degrees
    dxidx(3,1:9) = dxidx(3,1:9)*180.d0/PI
    xi(3)        = xi(3)*180.d0/PI
            
    return     
end subroutine cartesian2jacobi
        
subroutine angle(a,b,c,phi,dphida,dphidb,dphidc)
!This subroutine calculates the defined by the points a,b,c and the derivative of the angle
!with respect to these three points (all in radians)
    implicit none
    real*8,intent(in) ::a(3),b(3),c(3)
    real*8,intent(out)::phi
    real*8,intent(out),optional::dphida(3),dphidb(3),dphidc(3)
    real*8 :: y(3),z(3),y2,z2,ny,nz,yz,nynz,y0z0
    real*8 :: dcosphidy(3),dcosphidz(3),dcosdphi,dphidcos
    real*8 :: dphidy(3),dphidz(3)    
    ! y=a-b, z=c-b
    y   =  a-b
    z   =  c-b
    ! y2, z2
    y2  = sum(y**2) 
    z2  = sum(z**2) 
    ! abs(y), abs(z)
    ny  = dsqrt(y2)
    nz  = dsqrt(z2)
    ! yz
    yz  = sum(y*z)
    ! abs(y)abs(z)
    nynz =ny*nz   
    ! y0*z0

    y0z0=   yz/nynz
    
    y0z0=   min(1d0,max(-1d0,y0z0))         !=cosphi
    phi =   dacos(y0z0)
          
    if (.not.present(dphida).or..not.present(dphidb).or..not.present(dphidc)) return
          
    ! calculate derivatives
    ! y=a-b, z=c-b
    dcosphidy(1:3)=(z(1:3)-y(1:3)/y2*yz)/nynz 
    dcosphidz(1:3)=(y(1:3)-z(1:3)/z2*yz)/nynz 
    ! dcosdphi
    dcosdphi = -dsin(phi)                    
    ! dphidcos 

    !THIS WAS ADDED BY OLIVER TO CIRCUMVENT DIVISION BY ZERO!!! IMPORTANT!!!
    if(dabs(dcosdphi) < epsilon(0.d0)) dcosdphi = epsilon(0.d0)
    
    dphidcos =1d0/dcosdphi 
    !print*, "dphidcos", dphidcos
    
    ! dphidy dphidz
    dphidy(1:3) = dphidcos*dcosphidy
    dphidz(1:3) = dphidcos*dcosphidz   
    ! dphida= dphidy*dyda = dphidy
    dphida =  dphidy
    ! dphidc= dphidz*dzdc = dphidz
    dphidc =  dphidz
    ! dphida= dphidy*dydb + dphidz*dzdb = - dphidy - dphidz 
    dphidb = -dphida-dphidc
    
    return
end subroutine angle  
        
real*8 function distance(a,b)
    implicit none 
    !real*8 ::distance
    real*8,intent(in) ::a(3),b(3)
    distance=dsqrt(sum((a-b)*(a-b)))
end function distance
        
subroutine calc_weights(w,w0,dw0dr,dwdr,dwdw0,r)
    implicit none
    real*8, dimension(3),   intent(in)  :: r
    real*8, dimension(3),   intent(out) :: w,w0,dw0dr
    real*8, dimension(3,3), intent(out) :: dwdr
    real*8, dimension(3,3), intent(out) :: dwdw0
    
    real*8, dimension(3), save :: w_save = (/1.d0,0.d0,0.d0/) !this is a safety measure, so that we never get w=0 for all weights!
    
    real*8 :: sumw0
    
    !for loops
    integer i
    
    do i=1,3
        call switchfac(r(i),w0(i),dw0dr(i))
    end do
    sumw0 = sum(w0)
    !this is implemented to prevent division by zero in case we get full atomization
    if(sumw0 < epsilon(0.d0)) sumw0 = 10*epsilon(0.d0)
    
    !calculate the normalized weights
    w = w0/sumw0
    
    !check if we have the unfortunate case of all weights being 0. If yes, we load the weights from the previous step.
    if((w(1) == 0.d0).and.(w(2) == 0.d0).and.(w(3) == 0.d0)) then
        w = w_save
    end if
    
    !store old weights
    w_save = w
    
    !calculate derivatives of the normalized weights with respect to all different r
    dwdr(1,1)  = (dw0dr(1)*sumw0 - w0(1)*dw0dr(1))/(sumw0**2)
    dwdr(1,2)  = (-w0(1)*dw0dr(2))/(sumw0**2)
    dwdr(1,3)  = (-w0(1)*dw0dr(3))/(sumw0**2)
    
    dwdr(2,1)  = (-w0(2)*dw0dr(1))/(sumw0**2)
    dwdr(2,2)  = (dw0dr(2)*sumw0 - w0(2)*dw0dr(2))/(sumw0**2)
    dwdr(2,3)  = (-w0(2)*dw0dr(3))/(sumw0**2)
    
    dwdr(3,1)  = (-w0(3)*dw0dr(1))/(sumw0**2)
    dwdr(3,2)  = (-w0(3)*dw0dr(2))/(sumw0**2)
    dwdr(3,3)  = (dw0dr(3)*sumw0 - w0(3)*dw0dr(3))/(sumw0**2)
    
    dwdw0(1,1) = (sumw0-w0(1))/(sumw0**2)
    dwdw0(1,2) = (-w0(1))/(sumw0**2)
    dwdw0(1,3) = (-w0(1))/(sumw0**2)
    
    dwdw0(2,1) = (-w0(2))/(sumw0**2)
    dwdw0(2,2) = (sumw0-w0(2))/(sumw0**2)
    dwdw0(2,3) = (-w0(2))/(sumw0**2)
    
    dwdw0(3,1) = (-w0(3))/(sumw0**2)
    dwdw0(3,2) = (-w0(3))/(sumw0**2)
    dwdw0(3,3) = (sumw0-w0(3))/(sumw0**2)

    return
end subroutine calc_weights   

subroutine switchfac(x,f,dfdx)
!This is the subroutine to calculate the switching factors
  implicit none

  !input
  real*8, intent(in) :: x

  !output
  real*8, intent(out) :: f
  real*8, intent(out), optional :: dfdx

  !local
  !####CHANGE CUTOFF FOR SWITCHING FUNCTION HERE####
  real*8, parameter :: dx = 1.2d0 !determines the "cutoff" of the switching
  real*8, parameter :: a  = 4.d0 !determines plateau of switching function

!HERE THE SWITCHING FUNCTION IS SPECIFIED. BE SURE TO ALSO CHANGE THE DERIVATIVE
!SHOULD YOU CHANGE IT!!!

!new function
  f= dexp(-(x/dx)**a)

  if (present(dfdx))  &
        dfdx= -(a/x)*f*(x/dx)**a
         

end subroutine      
        
end module utility_conversion

!*****************************************************************************************
!    EVERYTHING BETWEEN THESE *** MARKERS WAS ADDED! EVERYTHING ELSE IN THIS SUBROUTINE IS
!    VANILLA CHARMM CODE! YOU SHOULD BE ABLE TO RETURN TO VANILLA CODE BY COMMENTING ALL
!    THE STUFF BETWEEN THESE MARKERS OUT. IMPORTANT: UNCOMMENT EU=0.0 LINE RIGHT AFTER
!    THE VARIABLE DECLARATIONS IF YOU WANT VANILLA CHARMM CODE (THIS WAS THERE ORIGINALLY)
!*****************************************************************************************
module usermod
  logical :: printe_flag = .false., printe_off = .false.
  integer :: iuene

contains

SUBROUTINE USERSB
  !
  !     THIS DOES NOTHING SAVE TELL USER THAT HE CALLED A NULL PROGRAM.
  !
  !     Author: Robert Bruccoleri
  !
 !!! use chm_kinds
 !!! use stream
 !!! use comand
 !!! use string
  implicit none
  !
!!!  CALL WRNDIE(1,'<USERSB>','USERSB CALLED. NO ACTION TAKEN.')
  RETURN
END SUBROUTINE USERSB

SUBROUTINE USERE(EU,X,Y,Z,DX,DY,DZ,QECONT,ECONT,NATOMX)
  !
  !     THE USER ENERGY ROUTINE. DOES NOTHING IN THE NORMAL SYSTEM
  !     EXCEPT SET EU TO ZERO. QECONT IS A FLAG WHICH WILL BE SET TO 1
  !     WHEN THE ROUTINE IS CALLED IN THE QECONTIS SECTION. INDIVIDUAL
  !     ATOM USER ENERGIES ARE THEN RETURNED IN ECONT, AND THE DERIVATIVE
  !     ARRAYS SHOULD NOT BE ACCESSED.
  !
  !     EU   - User energy to be returned
  !     X,Y,Z - Coordinates for energy evaluation
  !     DX,DY,DZ - Forces. Must be modified. Append (dE/dX,...)
  !     QECONT - Flag for analysis (0-no analysis,>0=analysis)
  !     ECONT(NATOMX) - Analysis array to be filled if QECONT>0.
  !     NATOMX - Number of atoms
  !
  !     Author: Robert Bruccoleri
  !
  use chm_kinds
!*****************************************************************************************
! BEGIN OF MODIFICATIONS
!*****************************************************************************************
  use CustomEnergy ,&
      only: Epot, dEdR, dEdVdW, dEdalpha, adEdR, adEdVdw, adEdalpha
  use utility_conversion  
  
!*****************************************************************************************
! END OF MODIFICATIONS
!*****************************************************************************************  
  implicit none
  real(chm_real) EU
  INTEGER NATOMX
  real(chm_real) X(NATOMX),Y(NATOMX),Z(NATOMX)
  real(chm_real) DX(NATOMX),DY(NATOMX),DZ(NATOMX)
  LOGICAL QECONT
  real(chm_real) ECONT(NATOMX)
  !
  INTEGER I
  !
  !EU=0.0  !***THIS WAS ORIGINALLY IN CHARMM, BUT I COMMENTED OUT TO CONTINUE DECLARATIONS
  
!*****************************************************************************************
! BEGIN OF MODIFICATIONS
!*****************************************************************************************

  !loop variable (integer i is already defined by standard charmm code)
  integer j
  
  !conversion from eV to kcal/mol
  real*8, parameter :: ev2kcalmol = 23.06035d0 
  
  !Array that stores masses. This is needed for some things like calculating the centre of
  !mass, needed to calculate Jacobi coordinates. YOU NEED TO CHANGE THIS WHERE IT'S DEFINED!
  real*8, dimension(NATOMX)  :: mass
  
  !Output coordinates (Jacobi)
  !1st index: surface
  !2nd index: Jacobi coordinate (r, bigR, alpha)
  real*8, dimension(3,3) :: jacCoord
  
  !Array for the derivatives of Jacobi coordinates with respect to the cartesian coordinates
  !1st index: surface
  !2nd index: Jacobi coordinate (r,bigR,alpha)
  !3rd index: cartesian coordinate (x1,y1,z1,x2,y2,z2,x3,y3,z3)
  real*8, dimension(3,3,9) :: dJacdCart
  
  !Array for the derivatives of the potential with respect to the Jacobi coordinates 
  !1st index: surface
  !2nd index: Jacobi coordinate (r,bigR,alpha)
  real*8, dimension(3,3) :: dEdJac
  
  !All this is just needed for the PES switching routine
  !weights for the 3 surfaces, w0: unnormalized, w: normalized
  real*8, dimension(3)   :: w0,w  
  !derivatives of raw (unnormalized) weights with respect to r
  real*8, dimension(3)   :: dw0dr 
  !derivatives of weight j with respect to r_i (1st index j, 2nd index i)
  real*8, dimension(3,3) :: dwdr
  !derivatives of weight j with respect to raw weight i (1st index j, 2nd index i)
  real*8, dimension(3,3) :: dwdw0
  !derivatives of weight j with respect to cart. coord. i (1st index j, 2nd index i)
  real*8, dimension(3,9) :: dwjdxi 
  !derivatives of potential j with respect to cart. coord. i (1st index j, 2nd index i)
  real*8, dimension(3,9) :: dEjdxi
  !temporary array that stores either dwjdxi or dEjdxi, because I need to swap back coordinates
  !for surface 2 and 3!
  real*8, dimension(9) :: ddxi
  !potential energies on the "pure" surfaces
  real*8, dimension(3)   :: Esurf
                                            
!-------------------------------------CHANGE HERE-------------------------------------  
! If this value is set to .true., reactive MD simulations are performed (the PESs are switched) 
! dynamically. If it is set to .false., the simulations will be performed on PES 1 only (bond 
! dissociation of the diatomic is then not possible).
  logical :: use_surface_switching = .true.
!-------------------------------------------------------------------------------------  
  
  !####CHOOSE HERE HOW THE DERIVATIVES OF THE PES SHOULD BE CALCULATED!####
  !####Analytical derivatives were found to be superior in this implementation####
  logical :: use_analytical_derivatives = .true.
  
  !THIS IS AN ERROR CHECK. IF WE DON'T HAVE 3 ATOMS, THIS CODE MAKES NO SENSE AT ALL!
  if(NATOMX /= 3) then
      write(*,*) "WARNING: Number of atoms is not 3."
      write(*,*) "This code is designed for 3-atomic systems ONLY."
      write(*,*) "program was terminated."
      stop
  end if
    
!-------------------------------------CHANGE HERE-------------------------------------        
! Masses of the atoms need to be specified here, because we need to calculate the centre of mass.
! The unit is not important, as long as all masses are given in the same unit.
! The reordering of the masses for the other surfaces is performed automatically.
  mass(1) = 14.007d0  !this corresponds to atom 1 
  mass(2) = 15.9994d0 !this corresponds to atom 2 
  mass(3) = 14.007d0 !this corresponds to atom 3 
!-------------------------------------------------------------------------------------     
  
  call calc_Jacobi_coordinates(X,Y,Z,jacCoord,mass,dJacdCart)
!DEBUG: Jacobi coordinates for all 3 surfaces are calculated fine (was tested)
!and deemed bug-free
!  print*, "Jacobi coordinates: ", "    r    ", "    R    ", "    alpha    "
!  print*, "(surface 1)       : ", jacCoord(1,1), jacCoord(1,2), jacCoord(1,3)
!  print*, "(surface 2)       : ", jacCoord(2,1), jacCoord(2,2), jacCoord(2,3)
!  print*, "(surface 3)       : ", jacCoord(3,1), jacCoord(3,2), jacCoord(3,3)
  
  !print*, "dJacdCart:", dJacdCart
 

if(use_surface_switching) then
!calculate weights and their respective derivatives etc.
    call calc_weights(w,w0,dw0dr,dwdr,dwdw0,jacCoord(:,1))
else
!give surface 1 the full weight
    w    = 0.d0
    w(1) = 1.d0
end if
 
!calculate the derivatives of the potential surfaces 1-3 with respect to the Jacobi
!coordinates (r,R,alpha).
 do i=1,3
     !if a weight is 0, we can just set all the derivatives to 0 and skip expensive
     !calculations
     if(w(i) == 0.d0) then
         dEdJac(i,:) = 0.d0
         cycle 
     end if
     
     if(.not.use_analytical_derivatives) then
     !do it numerically
         dEdJac(i,1) = dEdR(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
         dEdJac(i,2) = dEdVdW(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
         dEdJac(i,3) = dEdalpha(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)    
     else
     !do it analytically
         dEdJac(i,1) = adEdR(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
         dEdJac(i,2) = adEdVdW(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)
         dEdJac(i,3) = adEdalpha(jacCoord(i,1),jacCoord(i,2),jacCoord(i,3),i)   
     end if 
 end do

 
  !if we need to switch surfaces, we do it like this
  if(use_surface_switching) then
      !calculate energies on "pure" surfaces
      do j=1,3
          !if a weight is 0, we can just set the energy to 0 and skip expensive calculations
          if(w(j) == 0.d0) then
              Esurf(j) = 0.d0
              cycle 
          end if
          Esurf(j) = Epot(jacCoord(j,1),jacCoord(j,2),jacCoord(j,3),j)
      end do
      !Assign mixed energy value:
      EU = real((Esurf(1)*w(1) + Esurf(2)*w(2) + Esurf(3)*w(3))*ev2kcalmol,chm_real)
         
      !calculate all 18 dwjdxi and dEjdxi
      do j=1,3
          do i=1,9             
              dEjdxi(j,i) = sum(dEdJac(j,:)*dJacdCart(j,:,i))
              
              
              !calculate dwjdxi by using the raw weight derivatives
              dwjdxi(j,i) = dwdw0(j,1)*dw0dr(1)*dJacdCart(1,1,i) + dwdw0(j,2)*dw0dr(2)*dJacdCart(2,1,i) &
                           +dwdw0(j,3)*dw0dr(3)*dJacdCart(3,1,i)
              !print*, "dEjdxi: ", dEjdxi(j,i)
          end do
      end do
            
     !Assign derivatives:
     do i=1,3
         DX(i) = real(sum(w*dEjdxi(:,1+(i-1)*3) + dwjdxi(:,1+(i-1)*3)*Esurf)*ev2kcalmol,chm_real)
         DY(i) = real(sum(w*dEjdxi(:,2+(i-1)*3) + dwjdxi(:,2+(i-1)*3)*Esurf)*ev2kcalmol,chm_real)
         DZ(i) = real(sum(w*dEjdxi(:,3+(i-1)*3) + dwjdxi(:,3+(i-1)*3)*Esurf)*ev2kcalmol,chm_real)
     end do
                      
  !if not, we just consider the first surface and that's it  
  else
      !Assign energy:
      EU = real(Epot(jacCoord(1,1),jacCoord(1,2),jacCoord(1,3),1)*ev2kcalmol,chm_real)
      
      !Assign derivatives:
      do i=1,3
          DX(i) = real((dEdJac(1,1)*dJacdCart(1,1,1+(i-1)*3) + dEdJac(1,2)*dJacdCart(1,2,1+(i-1)*3) &
                   + dEdJac(1,3)*dJacdCart(1,3,1+(i-1)*3))*ev2kcalmol,chm_real)
          DY(i) = real((dEdJac(1,1)*dJacdCart(1,1,2+(i-1)*3) + dEdJac(1,2)*dJacdCart(1,2,2+(i-1)*3) &
                   + dEdJac(1,3)*dJacdCart(1,3,2+(i-1)*3))*ev2kcalmol,chm_real)
          DZ(i) = real((dEdJac(1,1)*dJacdCart(1,1,3+(i-1)*3) + dEdJac(1,2)*dJacdCart(1,2,3+(i-1)*3) &
                   + dEdJac(1,3)*dJacdCart(1,3,3+(i-1)*3))*ev2kcalmol,chm_real)
      end do 
  end if
     
!*****************************************************************************************
! END OF MODIFICATIONS
!*****************************************************************************************  

  IF(QECONT) THEN
     DO I=1,NATOMX
        ECONT(I)=0.0
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE USERE

SUBROUTINE USRACM(NATOMX,X,Y,Z,DX,DY,DZ,QECONT,ECONT,ICALL)
  !
  !     THE USER ACCUMULATION ROUTINE. DOES NOTHING IN THE NORMAL SYSTEM.
  !     This routine is called every time the energy is evaluated.
  !     The ICALL integer tells whether this is a "real" step for
  !     minimization or dynamics.
  !
  !     X,Y,Z    - Coordinates for energy evaluation
  !     DX,DY,DZ - Forces. Must be modified. Append (dE/dX,...)
  !     ICALL    - ECALLS increment 
  !
  !     Author: B. Brooks
  !
  use chm_kinds
  implicit none
  INTEGER NATOMX
  real(chm_real) X(NATOMX),Y(NATOMX),Z(NATOMX)
  real(chm_real) DX(NATOMX),DY(NATOMX),DZ(NATOMX)
  LOGICAL QECONT
  real(chm_real) ECONT(*)
  INTEGER ICALL
  !
  RETURN
END SUBROUTINE USRACM

SUBROUTINE USERNM(IUSER,X,Y,Z,XNORMA,YNORMA,ZNORMA,NATOMX)
  !
  !     THIS ROUTINE IS USED TO SETUP A GUESS NORMAL MODE, OR A MOTION
  !     OF INTEREST TO BE ANALYSED, OR APPENDED TO THE COORDINATE SET.
  !     THE CALLING SEQUENCE IS 3 ARRAYS WHICH ARE TO BE FILLED
  !     WITH COORDINATE DISPLACEMENTS. (I.E. DONT MASS WEIGHT THEM)
  !     SEE THE DOCUMENTATION OF VIBRAN FOR FURTHER DETAILS
  !
  !     IUSER  - Code for which function to use. This is only important
  !               if more than one user mode is specified.
  !     X,Y,Z   - Reference coordinates
  !     XNORMA,YNORMA,ZNORMA - Mode vector to be filled
  !     NATOMX  - Number of atoms
  !
  !     Author: Bernie Brooks
  !
  use chm_kinds
  implicit none
  INTEGER IUSER
  INTEGER NATOMX
  real(chm_real) X(NATOMX),Y(NATOMX),Z(NATOMX)
  real(chm_real) XNORMA(NATOMX),YNORMA(NATOMX),ZNORMA(NATOMX)
  !
  INTEGER I
  !
  !     The default user mode is the radius of gyration motion about
  !     the origin.
  !
  DO I=1,NATOMX
     XNORMA(I)=X(I)
     YNORMA(I)=Y(I)
     ZNORMA(I)=Z(I)
  ENDDO
  !
 !!! CALL WRNDIE(1,'<USERNM>', &
 !!!      'NO ROUTINE PROVIDED, RGYR MODE RETURNED')
  RETURN
END SUBROUTINE USERNM

SUBROUTINE USRSEL(NTAGS,FLAGS,X,Y,Z,QCOOR)
  !
  !     THIS ROUTINE ALLOWS A USER TO SELECT ATOMS.
  !
  !     NTAGS - number of entries (atoms)
  !     FLAGS(NTAGS) - array to hold selection values
  !     X,Y,Z  - coordinates
  !     QCOOR  - logical flag specifying that coordinates are present
  !              If .FALSE., then the coordinates should not be accessed.
  !
  !     If an atom is selected, set FLAGS value to be nonzero.
  !
  !
  use chm_kinds
  implicit none
  INTEGER NTAGS
  INTEGER FLAGS(NTAGS)
  real(chm_real) X(*),Y(*),Z(*)
  LOGICAL QCOOR
  !
  INTEGER I
  !
  !!! CALL WRNDIE(1,'<USRSEL>','NO USER SELECTION SPECIFIED')
  !
  !     THE DEFAULT USER SELECTION IS TO INCLUDE EVERYTHING
  DO I=1,NTAGS
     FLAGS(I)=1
  ENDDO
  RETURN
END SUBROUTINE USRSEL

SUBROUTINE USRTIM(SERVAL,QAT,NQ,ITIME, &
     NATOMX,X,Y,Z,XREF,YREF,ZREF,NSKIP,DELTA,TVAL,ISLCT)
  !
  !     This routine allows the user to specify a time series
  !     for use in CORREL.
  !
  !     SERVAL  - Value for time series specified in the ENTER command.
  !                May be used for distinguishing time series, or for any
  !                other purpose.
  !     QAT(NQ)- List of atoms that were specified in the ENTER command
  !     NQ     - Number of atoms specified in the ENTER command
  !     ITIME  - Which time step is currenlty being processed
  !     NATOMX - Number of atoms
  !     X,Y,Z  - Coordinates for current time step
  !     XREF,YREF,ZREF - Reference coordinates for time series
  !     DELTA  - Time interval between steps in picoseconds.
  !     TVAL   - Time series value returned
  !     ISLCT  - The first atom selection array from the TRAJ command
  !
  use chm_kinds
  implicit none
  !
  real(chm_real) SERVAL
  INTEGER ITIME,NATOMX,NQ
  INTEGER QAT(NQ)
  real(chm_real) X(NATOMX),Y(NATOMX),Z(NATOMX)
  real(chm_real) XREF(NATOMX),YREF(NATOMX),ZREF(NATOMX)
  INTEGER NSKIP
  real(chm_real) DELTA,TVAL
  INTEGER ISLCT(*)
  !
  !!! CALL WRNDIE(1,'<USECOR>','NO USER TIME SERIES SPECIFIED')
  !
  !     THE DEFAULT USER TIME SERIES IS SIMPLY THE TIME
  TVAL=ITIME/DELTA
  RETURN
END SUBROUTINE USRTIM

SUBROUTINE USRINI
  !
  !     This routine allows the user to specify a startup procedure.
  !     This routine is invoked at the beginning of each CHARMM run.
  !
!BEGIN************************************************************************************
  use CustomEnergy, &
      only : initialize_PES
!END**************************************************************************************
  use chm_kinds
  implicit none
  !
!BEGIN************************************************************************************
  call initialize_PES()
!END**************************************************************************************  
  RETURN
END SUBROUTINE USRINI

SUBROUTINE USPAFL(I,XX,YY,ZZ,XY,XZ,YZ, &
     XNORM,YNORM,ZNORM,NATOM,ISLCT)
  !
  !  This routine processes the user specific principal axis
  !  fluctuation code.  See documentation for details.
  !
  use chm_kinds
  implicit none
  INTEGER I,NATOM
  real(chm_real) XX(*),YY(*),ZZ(*),XY(*),XZ(*),YZ(*)
  real(chm_real) XNORM(*),YNORM(*),ZNORM(*)
  INTEGER ISLCT(*)
  !
  RETURN
END SUBROUTINE USPAFL

SUBROUTINE USER_RNGSPEC
  ! specification for USER random number generator

  IMPLICIT NONE

  RETURN
END SUBROUTINE USER_RNGSPEC


end module usermod


