! -*-f90-*-
!!-------------------------------------------------------------------------bl-! 
!! 
!! High-fidElity tool for maGnEto-gasdynamic simuLations (HEGEL) 
!! 
!! Copyright (C) 2016 Alessandro Munafo' (University of Illinois at Urbana-Champaign) 
!! 
!! This program is free software; you can redistribute it and/or 
!! modify it under the terms of the Version 2.1 GNU Lesser General 
!! Public License as published by the Free Software Foundation. 
!! 
!! This program is distributed in the hope that it will be useful, 
!! but WITHOUT ANY WARRANTY; without even the implied warranty of 
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
!! Lesser General Public License for more details. 
!! 
!! You should have received a copy of the GNU Lesser General Public 
!! License along with this library; if not, write to the Free Software 
!! Foundation, Inc. 51 Franklin Street, Fifth Floor, 
!! 
!! Boston, MA  02110-1301  USA 
!! 
!!-------------------------------------------------------------------------el-!
!-----------------------------------------------------------------------------!
!> This module provides useful subroutines for writing Tecplot files in binary format
!! for structured grids storing variables at cell-centroids (e.g., Cell-Centered Finite-Volume method)
   module utils_tecplot_binary

     use, intrinsic :: iso_c_binding
#ifdef TECIOMPI           
     use mpi  
#endif

     implicit none

     include 'tecio.f90'

     !----------------------!
     !> Default: make everything private
     private
     public :: open_tecplot_binary, close_tecplot_binary, write_field_tecplot_binary, &
             & write_field_line_tecplot_binary, write_nodes_tecplot_binary,           & 
             & create_zone_tecplot_binary, write_face_conn_tecplot_binary,            & 
             & tecplot_file_binary
#ifdef TECIOMPI
     public :: open_tecplot_parallel_binary, create_par_zone_tecplot_binary
#endif

     !----------------------!
     !> Precision parameters
     integer, parameter             :: hpc_4 = 4
     integer, parameter             :: hpc_8 = 8
     integer(kind=hpc_4), parameter :: s_str = 40 
     integer(kind=hpc_4), parameter :: l_str = 1000

     !----------------------!
     !> Tecplot parameters
     !> Debug flags
     integer(kind=hpc_4), parameter :: NO_DEBUG = 0
     integer(kind=hpc_4), parameter :: DEBUG    = 1

     !> File formats
     integer(kind=hpc_4), parameter :: PLT_FORMAT   = 0   !< Tecplot binary .plt
     integer(kind=hpc_4), parameter :: SZPLT_FORMAT = 1   !< Tecplot subzone .szplt

     !> File types
     integer(kind=hpc_4), parameter :: FULL     = 0  !< Grid + solution 
     integer(kind=hpc_4), parameter :: GRID     = 1  !< Grid only
     integer(kind=hpc_4), parameter :: SOLUTION = 2  !< Solution only

     !> Variable types
     integer(kind=hpc_4), parameter :: SINGLE = 0
     integer(kind=hpc_4), parameter :: DOUBLE = 1

     !> Value location
     integer(kind=hpc_4), parameter :: CELL_CENTERED = 0 
     integer(kind=hpc_4), parameter :: NODE_CENTERED = 1

     !> Face neighbor mode
     integer(kind=hpc_4), parameter :: LOCAL_ONE_TO_ONE   = 0
     integer(kind=hpc_4), parameter :: LOCAL_ONE_TO_MANY  = 1 
     integer(kind=hpc_4), parameter :: GLOBAL_ONE_TO_ONE  = 2 
     integer(kind=hpc_4), parameter :: GLOBAL_ONE_TO_MANY = 3

     !> Block or non-Block data
     integer(kind=hpc_4), parameter :: NO_BLOCK = 0
     integer(kind=hpc_4), parameter :: IS_BLOCK = 1

     !> Zone tyes
     integer(kind=hpc_4), parameter :: ORDERED = 0

     !> Variable types
     integer(kind=hpc_4), parameter :: ACTIVE  = 0
     integer(kind=hpc_4), parameter :: PASSIVE = 1

     !> Variable sharing
     integer(kind=hpc_4), parameter :: NO_SHARED = 0 
     integer(kind=hpc_4), parameter :: IS_SHARED = 1

     contains 

       !--------------------------------------------------!
       !> This subroutine enables to switch between different Tecplot files
       subroutine tecplot_file_binary(f)

         integer(kind=hpc_4), intent(in) :: f   !< Tecplot file ID 

         integer(kind=hpc_4)             :: ierr

         ierr = TECFIL142(f)

       end subroutine tecplot_file_binary

       !--------------------------------------------------!
       !> This subroutine opens a serial Tecplot binary file 
       subroutine open_tecplot_binary(tecplot_file, tecplot_fields)

         character(len=*), intent(in)               :: tecplot_file     !< Tecplot file name
         character(len=*), dimension(:), intent(in) :: tecplot_fields   !< Tecplot field variable names

         integer(kind=hpc_4)                        :: i, l1, l2, n
         integer(kind=hpc_4)                        :: ierr
         character(len=l_str)                       :: tecplot_var_list

         !------------------!
         ! Create var list 
         n = size(tecplot_fields)

         tecplot_var_list = trim(tecplot_fields(1))
         do i = 2,n
            l1 = len_trim(tecplot_var_list)
            l2 = len_trim(tecplot_fields(i))
            tecplot_var_list(l1 + 1:l1 + l2 + 1) = ' '//trim(tecplot_fields(i))
         enddo

         !------------------!
         ! Open Tecplot file - grid + solution are written 
         ierr = TECINI142('T'//C_NULL_CHAR,                         &
                          trim(tecplot_var_list)//C_NULL_CHAR,      &
                          trim(tecplot_file)//'.plt'//C_NULL_CHAR,  &
                          '.'//C_NULL_CHAR,                         &
                          PLT_FORMAT,                               &
                          FULL,                                     &
                          NO_DEBUG,                                 &
                          DOUBLE)

       end subroutine open_tecplot_binary

       !--------------------------------------------------!
       !> This subroutine opens a parallel Tecplot binary file 
#ifdef TECIOMPI       
       subroutine open_tecplot_parallel_binary(tecplot_file, tecplot_fields, master_rank, comm)

         character(len=*), intent(in)               :: tecplot_file     !< Tecplot file name
         character(len=*), dimension(:), intent(in) :: tecplot_fields   !< Tecplot field variable names
         integer(kind=hpc_4), intent(in)            :: master_rank      !< ID of MPI master rank/process
         integer, intent(in)                        :: comm             !< MPI communicator

         integer(kind=hpc_4)                        :: i, l1, l2, n
         integer(kind=hpc_4)                        :: ierr
         character(len=l_str)                       :: tecplot_var_list

         !------------------!
         ! Create var list 
         n = size(tecplot_fields)

         tecplot_var_list = ''
         tecplot_var_list = trim(tecplot_fields(1))
         do i = 2,n
            l1 = len_trim(tecplot_var_list)
            l2 = len_trim(tecplot_fields(i))
            tecplot_var_list(l1 + 1:l1 + l2 + 1) = ' '//trim(tecplot_fields(i))
         enddo

         !------------------!
         ! Open Tecplot file - grid + solution are written
         ierr = TECINI142('T'//C_NULL_CHAR,                          &
                          trim(tecplot_var_list)//C_NULL_CHAR,       &
                          trim(tecplot_file)//'.szplt'//C_NULL_CHAR, &
                          '.'//C_NULL_CHAR,                          &
                          SZPLT_FORMAT,                              &
                          FULL,                                      &
                          NO_DEBUG,                                  &
                          DOUBLE)      

         !------------------!
         ! Initialize the MPI writer
         ierr = TECMPIINIT142(comm, master_rank)    
      
       end subroutine open_tecplot_parallel_binary
#endif
       !--------------------------------------------------!
       !> This subroutine creates a Tecplot zone for serial output 
       subroutine create_zone_tecplot_binary(strand_id, var_loc, n_fields, n_face_conn, time, in_min, in_max, jn_min, jn_max, & 
                                           & kn_min, kn_max)

          integer(kind=hpc_4), intent(in)          :: strand_id     !< strand ID
          integer(kind=hpc_4), intent(in)          :: var_loc       !< variable location (node cenetered or cell-centered)
          integer(kind=hpc_4), intent(in)          :: n_fields      !< number of field variables
          integer(kind=hpc_4), intent(in)          :: n_face_conn   !< number of face connectivities
          real(kind=hpc_8), intent(in)             :: time          !< time [s]
          integer(kind=hpc_4), intent(in)          :: in_min        !< miminum global node index along x direction 
          integer(kind=hpc_4), intent(in)          :: in_max        !< maximum global node index along x direction
          integer(kind=hpc_4), intent(in)          :: jn_min        !< miminum global node index along y direction
          integer(kind=hpc_4), intent(in)          :: jn_max        !< maximum global node index along y direction
          integer(kind=hpc_4), intent(in)          :: kn_min        !< miminum global node index along z direction
          integer(kind=hpc_4), intent(in)          :: kn_max        !< maximum global node index along z direction

          integer(kind=hpc_4), parameter           :: zero     = 0
          integer(kind=hpc_4), parameter           :: IcellMax = zero
          integer(kind=hpc_4), parameter           :: JcellMax = zero
          integer(kind=hpc_4), parameter           :: KCellMax = zero
          integer(kind=hpc_4)                      :: ierr
          integer(kind=hpc_4)                      :: Nx, Ny, Nz
          integer(kind=hpc_4)                      :: ParentZn
          integer(kind=hpc_4)                      :: ShrConn
          integer(kind=hpc_4), dimension(n_fields) :: var_location, passive_var, shared_var

          !-----------------!
          ! Variable location
          ! - node centered for nodes
          ! - cell/node centered for other variables
          var_location(1:3) = NODE_CENTERED

          ! Safety check 
          if ((var_loc.ne.NODE_CENTERED).and.(var_loc.ne.CELL_CENTERED)) then
              write(*,5)'Tec360_binary_utils::ERROR!!'
              write(*,5)'Tec360_binary_utils:: variable location can be only CELL(=0) or NODE(=1) centered...'
              write(*,5)'Tec360_binary_utils:: please correct the variable location ID!'
              print*
              stop
          endif
          var_location(4:n_fields) = var_loc

          ! Passive variables 
          passive_var = ACTIVE

          ! Shared variables 
          shared_var = NO_SHARED

          ! Number of nodes along x, y and z directions
          Nx = in_max - in_min + 1
          Ny = jn_max - jn_min + 1
          Nz = kn_max - kn_min + 1

          ParentZn = 0; ShrConn  = 0

          ierr = TECZNE142('ZONE'//C_NULL_CHAR, &
                          FULL,                 &
                          Nx,                   &
                          Ny,                   &
                          Nz,                   &
                          ICellMax,             & 
                          JCellMax,             &
                          KCellMax,             &
                          time,                 &
                          strand_id,            &
                          ParentZn,             &
                          IS_BLOCK,             &
                          n_face_conn,          &
                          GLOBAL_ONE_TO_ONE,    &
                          zero,                 &
                          zero,                 &
                          zero,                 &
                          passive_var,          &
                          var_location,         &
                          shared_var,           &
                          ShrConn)
        
5 format(a)
 
       end subroutine create_zone_tecplot_binary
#ifdef TECIOMPI
       !--------------------------------------------------!
       !> This subroutine creates a Tecplot zone for parallel output 
       subroutine create_par_zone_tecplot_binary(proc_id, nb_proc, var_loc, n_fields, strand_id, time, Nx, Ny, Nz, & 
                                               & in_min, in_max, jn_min, jn_max, kn_min, kn_max)

          integer(kind=hpc_4), intent(in)          :: proc_id       !< processor ID
          integer(kind=hpc_4), intent(in)          :: nb_proc       !< number of processors within the communicator
          integer(kind=hpc_4), intent(in)          :: var_loc       !< variable location (node cenetered or cell-centered)
          integer(kind=hpc_4), intent(in)          :: n_fields      !< number of field variables
          integer(kind=hpc_4), intent(in)          :: strand_id     !< strand ID (0 for static zones)
          real(kind=hpc_8), intent(in)             :: time          !< time [s]
          integer(kind=hpc_4), intent(in)          :: Nx            !< number of global nodes along x direction
          integer(kind=hpc_4), intent(in)          :: Ny            !< number of global nodes along y direction 
          integer(kind=hpc_4), intent(in)          :: Nz            !< number of global nodes along z direction
          integer(kind=hpc_4), intent(in)          :: in_min        !< miminum local node index along x direction 
          integer(kind=hpc_4), intent(in)          :: in_max        !< maximum local node index along x direction
          integer(kind=hpc_4), intent(in)          :: jn_min        !< miminum local node index along y direction
          integer(kind=hpc_4), intent(in)          :: jn_max        !< maximum local node index along y direction
          integer(kind=hpc_4), intent(in)          :: kn_min        !< miminum local node index along z direction
          integer(kind=hpc_4), intent(in)          :: kn_max        !< maximum local node index along z direction

          integer(kind=hpc_4), parameter           :: zero     = 0
          integer(kind=hpc_4), parameter           :: IcellMax = zero
          integer(kind=hpc_4), parameter           :: JcellMax = zero
          integer(kind=hpc_4), parameter           :: KCellMax = zero
          integer(kind=hpc_4)                      :: p, partition
          integer(kind=hpc_4)                      :: ierr
          integer(kind=hpc_4)                      :: ParentZn
          integer(kind=hpc_4)                      :: ShrConn
          integer(kind=hpc_4)                      :: n_face_conn
          integer(kind=hpc_4), dimension(nb_proc)  :: Partition_Ranks
          integer(kind=hpc_4), dimension(n_fields) :: var_location, passive_var, shared_var

          !-----------------!
          ! Variable location
          ! - node centered for nodes
          ! - cell/node centered for other variables
          var_location(1:3) = NODE_CENTERED

          ! Safety check 
          if ((var_loc.ne.NODE_CENTERED).and.(var_loc.ne.CELL_CENTERED)) then
              write(*,5)'Tec360_binary_utils::ERROR!!'
              write(*,5)'Tec360_binary_utils:: variable location can be only CELL(=0) or NODE(=1) centered...'
              write(*,5)'Tec360_binary_utils:: please correct the variable location ID!'
              print*
              stop
          endif
          var_location(4:n_fields) = var_loc
           
          ! Passive variables 
          passive_var = ACTIVE

          ! Shared variables 
          shared_var = NO_SHARED

          ParentZn = zero;  ShrConn = zero;  n_face_conn = zero

          ierr = TECZNE142('ZONE'//C_NULL_CHAR, &
                          FULL,                 &
                          Nx,                   &
                          Ny,                   &
                          Nz,                   &
                          ICellMax,             & 
                          JCellMax,             &
                          KCellMax,             &
                          time,                 &
                          strand_id,            &
                          ParentZn,             &
                          IS_BLOCK,             &
                          n_face_conn,          &
                          LOCAL_ONE_TO_ONE,     &
                          zero,                 &
                          zero,                 &
                          zero,                 &
                          passive_var,          &
                          var_location,         &
                          shared_var,           &
                          ShrConn)
      
          ! Fill partition ranks with processors IDs
          partition = proc_id + 1
          do p = 0,nb_proc - 1
             Partition_Ranks(p + 1) = p 
          enddo
           
          ! Create partition maps - each rank writes its portion of the domain (grid + solution in general)
          ierr = TECZNEMAP142(nb_proc, Partition_Ranks)
          ierr = TECIJKPTN142(partition, in_min, jn_min, kn_min, in_max, jn_max, kn_max)
         
5 format(a)
 
       end subroutine create_par_zone_tecplot_binary
#endif
       !--------------------------------------------------!
       !> This subroutines writes a field to a Tecplot binary file 
       subroutine write_field_tecplot_binary(offset, field_data)

         integer(kind=hpc_4), dimension(:), intent(in)  :: offset       !< offset array      
         real(kind=hpc_8), dimension(:,:,:), intent(in) :: field_data   !< field data

         integer(kind=hpc_4)                            :: astat
         integer(kind=hpc_4)                            :: nx, ny, nz, ox, oy, oz
         integer(kind=hpc_4)                            :: i, j, k
         integer(kind=hpc_4)                            :: ierr
         real(kind=hpc_8), dimension(:), allocatable    :: buffer

         ! Limits and offsets
         ox = offset(1)
         oy = offset(2)
         oz = offset(3)

         nx = size(field_data,1) - 2*ox
         ny = size(field_data,2) - 2*oy
         nz = size(field_data,3) - 2*oz

         allocate(buffer(nx),stat=astat)
         if (astat.ne.0) stop "Error while allocating: 'buffer'"
        
         do k = 1,nz
            do j = 1,ny
               do i = 1,nx
                  buffer(i) = field_data(i + ox,j + oy,k + oz)
               enddo
               ierr = TECDATD142(nx, buffer(1:nx))
            enddo
         enddo

         deallocate(buffer,stat=astat)
         if (astat.ne.0) stop "Error while deallocaing: 'buffer'"

       end subroutine write_field_tecplot_binary

       !--------------------------------------------------!
       !> This subroutines writes a field (store in a one-dimensional array) to a Tecplot binary file 
       subroutine write_field_line_tecplot_binary(field_data)

         real(kind=hpc_8), dimension(:), intent(in) :: field_data   !< scalar field

         integer(kind=hpc_4)                        :: ierr

         ierr = TECDATD142(size(field_data), field_data)

       end subroutine write_field_line_tecplot_binary

       !--------------------------------------------------!
       !> This subroutine writes the local nodes on a Tecplot ascii file 
       subroutine write_nodes_tecplot_binary(rn)

         real(kind=hpc_8), dimension(:,:,:,:), intent(in) :: rn   !< node coordinates

         integer(kind=hpc_4)                              :: astat
         integer(kind=hpc_4)                              :: i_max, j_max, k_max
         integer(kind=hpc_4)                              :: d, i, j, k
         integer(kind=hpc_4)                              :: ierr
         real(kind=hpc_8), dimension(:), allocatable      :: buffer

         ! Limits
         i_max = size(rn,1);  j_max = size(rn,2);  k_max = size(rn,3)
          
         allocate(buffer(i_max),stat=astat)
         if (astat.ne.0) stop "Error while allocating: 'buffer'"

         ! Write x, y and z node coordinates
         do d = 1,3

            do k = 1,k_max
               do j = 1,j_max
                  do i = 1,i_max
                     buffer(i) = rn(i,j,k,d)
                  enddo
                  ierr = TECDATD142(i_max, buffer)
               enddo
            enddo

         enddo

         deallocate(buffer,stat=astat)
         if (astat.ne.0) stop "Error while deallocating: 'buffer'"

       end subroutine write_nodes_tecplot_binary

       !--------------------------------------------------!
       !> This subroutines write the face connectivity information to a Tecplot binary file 
       subroutine write_face_conn_tecplot_binary(face_conn)

         integer(kind=hpc_4), dimension(:), intent(in) :: face_conn    !< face_connectivity 

         integer(kind=hpc_4)                           :: ierr

         ierr = TECFACE142(face_conn)

       end subroutine write_face_conn_tecplot_binary

       !--------------------------------------------------!
       !> This subroutine closes the Tecplot binary file writer
       subroutine close_tecplot_binary()

         integer(kind=hpc_4) :: ierr

         ierr = TECEND142()
         
       end subroutine close_tecplot_binary

  end module utils_tecplot_binary
!-----------------------------------------------------------------------------!
