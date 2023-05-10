module gol_common
!---------------------------------------------------------------------
!
!  Common routines and functions for Conway's Game of Life
!
!---------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
#ifdef _OPENMP
    use omp_lib
#endif
#ifdef _MPI 
    include "mpif.h"
#endif
    integer, parameter :: NUMSTATES = 4
    integer, parameter :: CellState_ALIVE = 0
    integer, parameter :: CellState_DEAD = 1
    integer, parameter :: CellState_DYING = 2
    integer, parameter :: CellState_BORN = 3
    
    integer, parameter :: NUMVISUAL = 4
    integer, parameter :: VisualiseType_VISUAL_ASCII = 0
    integer, parameter :: VisualiseType_VISUAL_PNG = 1
    integer, parameter :: VisualiseType_VISUAL_OPENGL = 2 
    integer, parameter :: VisualiseType_VISUAL_NONE = 3

    integer, parameter :: NUMICS = 2
    integer, parameter :: ICType_IC_RAND = 0
    integer, parameter :: ICType_IC_FILE = 1

    integer, parameter :: NUMRULES = 3
    integer, parameter :: RuleType_RULE_STANDARD = 0
    integer, parameter :: RuleType_RULES_EXTENDED = 1
    integer, parameter :: RuleType_RULES_PROB = 2

    integer, parameter :: NUMNEIGHBOURCHOICES = 2
    integer, parameter :: NeighbourType_NEIGHBOUR_STANDARD = 0
    integer, parameter :: NeighbourType_NEIGHBOUR_EXTENDED = 1
    
    integer, parameter :: NUMBOUNDARYCHOICES = 4
    integer, parameter :: BoundaryType_BOUNDARY_HARD = 0
    integer, parameter :: BoundaryType_BOUNDARY_TORAL = 1
    integer, parameter :: BoundaryType_BOUNDARY_TORAL_X_HARD_Y = 2
    integer, parameter :: BoundaryType_BOUNDARY_TORAL_Y_HARD_X = 3
    
    type Options 
        integer :: n, m, nsteps 
        integer :: iictype
        integer :: ivisualisetype
        integer :: iruletype
        integer :: ineighbourtype
        integer :: iboundarytype
        character(len=2000) :: statsfile
    end type Options 

contains 

    ! get the memory used as reported by system for process
    subroutine report_memory_usage()
        implicit none
        character(len=2000) :: stat_file
        integer :: ios 
        integer, parameter :: read_unit = 99
        character(len=200) :: line
        character(len=200), dimension(4) :: searchstr
        integer :: i, ii
        integer, dimension(4) :: mem_usage
    ! apple does not provide /proc/self/status file so 
    ! just ignore mem usage 
#ifndef __APPLE__
        write(*,*) "Mem reporting still in dev"
        stat_file = "/proc/self/status";
        searchstr(1) = "VmSize:"
        searchstr(2) = "VmPeak:"
        searchstr(3) = "VmRSS:"
        searchstr(4) = "VmHWM:"
        open (unit = read_unit, file = stat_file, iostat = ios, status = "old")
        if ( ios /= 0 ) then 
            write(*,*) "========================================================= "
            write(*,*) "Couldn't open", stat_file, " for memory usage reading"
            write(*,*) "========================================================= "
            return 
        end if
        do
            read(read_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            do i=1,4
                ii = index(line, searchstr(i))
                ! would need to add parsing of line to extract value of memory 
                if (ii > 0) then
                    mem_usage(ii) = 0
                end if 
            end do 
        end do
        close(read_unit)
        write(*,*) "========================================================= "
        write(*,*) "Memory Usage (VM, Peak VM, RSS, Peak RSS) = \n"
        write(*,*) mem_usage(1), mem_usage(2), mem_usage(3), mem_usage(4), "kilobytes"
        write(*,*) "========================================================= "
#endif
    end subroutine


    !   ascii visualisation
    subroutine visualise_ascii(step, grid, n, m)
        implicit none 
        integer, intent(in) :: step, n, m
        integer, dimension(:,:), intent(in) :: grid
        character :: cell 
        integer :: i, j

        write(*,*) "Game of Life"
        write(*,*) "Step ", step
        do i = 1, n
            do j = 1, m
                cell = ' '
                ! could use where 
                if (grid(i,j) .eq. CellState_ALIVE) cell = '*'
                write(*,"(A)", advance="no") cell
            end do 
            write(*,*) ""
        end do 
    end subroutine 

    ! png visualisation
    subroutine visualise_png(step, grid, n, m)
        implicit none 
        integer, intent(in) :: step, n, m
        integer, dimension(:,:), intent(in) :: grid

    end subroutine 

    ! no visualisation
    subroutine visualise_none(step)
        implicit none 
        integer, intent(in) :: step
        write(*,*) "Game of Life, Step ", step
    end subroutine 

    ! visualisation routine
    subroutine visualise(ivisualisechoice, step, grid, n, m)
        implicit none 
        integer, intent(in) :: ivisualisechoice 
        integer, intent(in) :: step, n, m
        integer, dimension(:,:), intent(inout) :: grid
        if (ivisualisechoice .eq. VisualiseType_VISUAL_ASCII) then 
            call visualise_ascii(step, grid, n, m)
        else if (ivisualisechoice .eq. VisualiseType_VISUAL_PNG) then 
            call visualise_png(step, grid, n, m)
        else  
            call visualise_none(step)
        end if 

    end subroutine

    ! generate random IC
    subroutine generate_rand_IC(grid, n, m)
        implicit none 
        integer, intent(in) :: n, m 
        integer, dimension(:,:), intent(inout) :: grid
        real :: xrand, rand
        integer :: i, j
        do i = 1, n
            do j = 1, m
#if defined(_CRAYFTN) || defined(_INTELFTN)
                call RANDOM_NUMBER(xrand)
#else
                xrand=rand()
#endif
                if (xrand .lt. 0.4) then 
                    grid(i,j) = CellState_DEAD
                else 
                    grid(i,j) = CellState_ALIVE
                end if 
            end do 
        end do 

    end subroutine 

    ! generate IC
    subroutine generate_IC(ic_choice, grid, n, m)
        implicit none 
        integer, intent(in) :: ic_choice
        integer, intent(in) :: n, m 
        integer, dimension(:,:), intent(inout) :: grid
        if (ic_choice .eq. ICType_IC_RAND) then 
            call generate_rand_IC(grid, n, m)
        end if 
    end subroutine 


    ! Report run state 
    subroutine run_state(opt)
        implicit none 
        type(Options), intent(inout) :: opt
        integer :: count
        integer*8 :: nbytes
        real*4 :: memfootprint
        integer :: mpi_size, omp_size, ierror

        nbytes = sizeof(opt%n) * opt%n * opt%m
        memfootprint = real(nbytes)/1024.0/1024.0/1024.0
        write(*,*) "========================================================= "
        write(*,*) " GOL Running with following "
        write(*,*) "========================================================= "
        write(*,*) "Requesting grid size of ", opt%n, opt%m
        write(*,*) " which requires", memfootprint, " GB "
        write(*,*) "Number of steps ", opt%nsteps
        write(*,*) "IC type ",opt%iictype
        write(*,*) "Visualization type ",opt%ivisualisetype
        write(*,*) "Neighbour rule type ",opt%ineighbourtype
        write(*,*) "Boundary rule type ",opt%iboundarytype
        write(*,*) "========================================================= "
#ifndef USEPNG
        if (opt%ivisualisetype .eq. VisualiseType_VISUAL_PNG) then 
            write(*, *) "PNG visualisation not enabled at compile time,"
            write(*, *) "turning off visualisation from now on."
        end if
#endif
    
#ifdef _MPI
        call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierror)
        write(*,*)a "========================================================= "
        write(*,*) "Running with MPI"
        write(*,*) "Comm World contains ", mpi_size, " processes "
        write(*,*) "========================================================= \n"
#endif
#ifdef _OPENMP
        omp_size = omp_get_max_threads()
        write(*, *) "========================================================= "
        write(*, *) "Running with OpenMP, version ", _OPENMP
        write(*, *) "Maximum number of threads %d \n", omp_size
        write(*, *) "========================================================= \n"
#endif
    end subroutine

    ! UI
    subroutine getinput(opt)
        implicit none 
        type(Options), intent(inout) :: opt
        character(len=2000) :: cmd
        character(len=32) :: arg 
        character(len=2000) :: statsfilename
        integer :: count
        integer*8 :: nbytes
        real*4 :: memfootprint
        ! get the commands passed and the number of args passed 
        call get_command(cmd)
        count = command_argument_count()
        if (count .lt. 2) then 
            write(*,*) "Usage: <grid height> <grid width> "
            write(*,*) "[<nsteps> <IC type> <Visualisation type> <Rule type> <Neighbour type> "
            write(*,*) "<Boundary type> <stats filename> ]"
            call exit();
        end if 
        
        statsfilename = "GOL-stats.txt"
        call get_command_argument(1,arg)
        read(arg,*) opt%n
        call get_command_argument(2,arg)
        read(arg,*) opt%m
        opt%nsteps = -1
        opt%ivisualisetype = VisualiseType_VISUAL_ASCII
        opt%iruletype = RuleType_RULE_STANDARD
        opt%iictype = ICType_IC_RAND
        opt%ineighbourtype = NeighbourType_NEIGHBOUR_STANDARD
        opt%iboundarytype = BoundaryType_BOUNDARY_HARD
        if (count .ge. 3) then 
            call get_command_argument(3,arg)
            read(arg,*) opt%nsteps
        end if 
        if (count .ge. 4) then 
            call get_command_argument(4,arg)
            read(arg,*) opt%iictype
        end if 
        if (count .ge. 5) then 
            call get_command_argument(5,arg)
            read(arg,*) opt%ivisualisetype
        end if 
        if (count .ge. 6) then
            call get_command_argument(6,arg)
            read(arg,*) opt%iruletype
        end if 
        if (count .ge. 7) then
            call get_command_argument(7,arg)
            read(arg,*) opt%ineighbourtype
        end if 
        if (count .ge. 8) then
            call get_command_argument(8,arg)
            read(arg,*) opt%iboundarytype
        end if 
        if (count .ge. 9) then
            call get_command_argument(9,arg)
            read(arg,*) statsfilename
        end if 
        if (opt%n .le. 0 .or. opt%m .le. 0) then
            write(*,*) "Invalid grid size."
            call exit(1)
        end if 
        opt%statsfile = statsfilename
        call run_state(opt)
    end subroutine 

    ! get some basic timing info
    real*8 function init_time()
        integer, dimension(8) :: value
        call date_and_time(VALUES=value)
        init_time = value(5)*3600.0+value(6)*60.0+value(7)+value(8)/1000.0
        return 
    end function
    ! get the elapsed time relative to start
    subroutine get_elapsed_time(start)
        real*8, intent(in) :: start
        real*8 :: finish, delta
        integer, dimension(8) :: value
        call date_and_time(VALUES=value)
        finish = value(5)*3600.0+value(6)*60.0+value(7)+value(8)/1000.0
        delta = finish - start
        write(*,*) "Elapsed time is ", delta, "s"
    end subroutine
    
end module
