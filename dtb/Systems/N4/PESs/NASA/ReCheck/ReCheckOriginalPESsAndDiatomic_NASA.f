program ComputeDiat

  implicit none 
      
  double precision :: X(4),   Y(4),  Z(4)
  double precision :: E
  double precision :: dEdX(4), dEdY(4), dEdZ(4)
  double precision :: Rp(6)
  double precision :: dVdR(6)
  double precision :: Qp(12)
  double precision :: dVdQ(12)
  
  integer                    :: Unit, Status
  character(200)             :: FileName
  integer                    :: i, j, iA
  double precision           :: RStart, REnd, hGrid
  double precision           :: RStart1, REnd1, hGrid1
  double precision           :: RStart2, REnd2, hGrid2
  integer                    :: NPoints, NPoints1, NPoints2
  double precision           :: alpha, Theta4
  character(3)               :: alphaChar
  double precision           :: alphaVec(5)

  double precision,parameter :: Pi  = acos( - 1.d0 )
  double precision,parameter :: Cconv = 0.52917721092d0
  double precision,parameter :: Econv = 0.159360144d-2
  double precision,parameter :: RConverter = 1.d0

 
!  
!  X    = [0.0d0,     0.0d0,   1.d10]
!  Y    = [0.0d0,     0.0d0,   1.d10]    
!  Z    = [0.0d0,     0.0d0,   1.d10]  
!  
!  NPoints = 10000
!  RStart  = 1.5d0
!  REnd    = 10.d0
!  hGrid = (REnd - RStart) / (NPoints-1)
!  X(2)  = RStart
!  
!  FileName = './DiatPot.dat'
!  open( File=trim(adjustl(FileName)), NewUnit=Unit, status='REPLACE', iostat=Status )
!  
!    do i=1,NPoints
!    
!      call pot(X,Y,Z,E,dEdX,dEdY,dEdZ)
!      
!      write(Unit,'(3d20.10)') X(2), E, dEdX(2)
!      
!      X(2) = X(2) + hGrid
!      
!    end do
!    
!  close(Unit)
!  
  
  NPoints1 = 200
  NPoints2 = 200
  RStart1  = 1.5d0
  REnd1    = 15.d0
  RStart2  = 1.5d0
  REnd2    = 15.d0
  hGrid1   = (REnd1 - RStart1) / (NPoints1-1)
  hGrid2   = (REnd2 - RStart2) / (NPoints2-1)

  
  alphaVec = [30.0, 60.0, 90.0, 120.0, 150.0]
  
  do iA=1,size(alphaVec,1)
    alpha = alphaVec(iA)
    write(alphaChar,'(I3)') int(alpha)
  
    FileName = './PESFromGrid.csv.' // trim(adjustl(alphaChar))
    open( File=trim(adjustl(FileName)), NewUnit=Unit, status='REPLACE', iostat=Status )
      write(Unit,'(A)') 'Variables = "r1", "r2", "r3", "r4", "r5", "r6", "E"'
      
      Rp    = 1.e3
      Rp(1) = RStart1
      do i=1,NPoints1
        Rp(2) = RStart2
        do j=1,NPoints2
          Theta4 = alpha / 180.d0 * Pi
          Rp(4)  = sqrt( Rp(1)**2 + Rp(2)**2 - 2.d0 * Rp(1) * Rp(2) * dcos(Theta4) ) * RConverter 

          Qp =[0.0, 0.0, 0.0, Rp(1), 0.0, 0.0, Rp(2)*cos(Theta4), Rp(2)*sin(Theta4), 0.0, 0.0, RpInf, 0.0] 
          
          R(1) = Sqrt( (  Qp(4) -  Qp(1) )**2 + (  Qp(5) -  Qp(2) )**2 + (  Qp(6) -  Qp(3) )**2)
          R(2) = Sqrt( (  Qp(7) -  Qp(1) )**2 + (  Qp(8) -  Qp(2) )**2 + (  Qp(9) -  Qp(3) )**2)
          R(3) = Sqrt( ( Qp(10) -  Qp(1) )**2 + ( Qp(11) -  Qp(2) )**2 + ( Qp(12) -  Qp(3) )**2)
          R(4) = Sqrt( (  Qp(4) -  Qp(7) )**2 + (  Qp(5) -  Qp(8) )**2 + (  Qp(6) -  Qp(9) )**2)
          R(5) = Sqrt( (  Qp(4) - Qp(10) )**2 + (  Qp(5) - Qp(11) )**2 + (  Qp(6) - Qp(12) )**2)
          R(6) = Sqrt( (  Qp(7) - Qp(10) )**2 + (  Qp(8) - Qp(11) )**2 + (  Qp(9) - Qp(12) )**2)

          E    = 0.0
          !call n4fitd(Qp(1:3), Qp(4:6), Qp(7:9), Qp(10:12), V, dVdQ)
          write(Unit,'(e20.10,6(A,e20.10))') Rp(1), ',', Rp(2), ',', Rp(3), ',', Rp(4), ',', Rp(5), ',', Rp(6), ',', E
          
          Rp(2) = Rp(2) + hGrid2
        end do
        Rp(1) = Rp(1) + hGrid1
      end do
      
    close(Unit)
    
  end do
  
  
end program ComputeDiat
