PROGRAM smesh_run

    USE smesh_io
    USE smesh

    IMPLICIT NONE

    INTEGER                              :: problem_type
    CHARACTER(len=200)                   :: filename_points
    REAL(8)                              :: xmin, xmax, ymin, ymax, dphi, drho_k, drho_m
    INTEGER                              :: nx, ny, nnode, nelem, n, voronoi_mesh_type
    LOGICAL                              :: shuffle
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: pt
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ve
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ne
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: dualb
    INTEGER, DIMENSION(:),   ALLOCATABLE :: dual
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: edge
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: vor_pt
    INTEGER, DIMENSION(:),   ALLOCATABLE :: vor_ve
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: vor_veb
    CHARACTER(len=200) :: filename_config, output_path
    LOGICAL :: using_default_config_file

    CALL parse_command_line(filename_config, using_default_config_file)

    ! --- Read settings from a text file ---
    OPEN(unit=21, file=trim(filename_config), status="old", access="sequential", action="read")
    READ(21,*) ! output path
    READ(21,*) output_path
    READ(21,*) ! input data
    READ(21,*) problem_type
    READ(21,*) filename_points
    READ(21,*) ! mesh
    READ(21,*) shuffle
    READ(21,*) ! voronoi
    READ(21,*) voronoi_mesh_type
    READ(21,*) ! simple grid settings
    READ(21,*) xmin
    READ(21,*) xmax
    READ(21,*) ymin
    READ(21,*) ymax
    READ(21,*) nx
    READ(21,*) ny
    READ(21,*) ! flower settings
    READ(21,*) dphi
    READ(21,*) drho_k
    READ(21,*) drho_m
    READ(21,*) n
    CLOSE(21)
     
    IF (problem_type .eq. -1) THEN
        CALL mesh_bisected_rectangle(pt, ve, xmin, xmax, ymin, ymax, nx, ny)
    ELSE IF (problem_type .eq. 0) THEN
        CALL load(filename_points, pt, .true.)
        ! CALL load_bin(filename_points, pt)
    ELSE IF (problem_type .eq. 1) THEN
        CALL mesh_bisected_rectangle(pt, ve, xmin, xmax, ymin, ymax, nx, ny)
        DEALLOCATE(ve)
    ELSE IF (problem_type .eq. 2) THEN
        CALL mesh_basic(pt, xmin, xmax, ymin, ymax, nx, ny)
    ELSE
        CALL mesh_flower(pt, dphi, drho_k, drho_m, n)
    END IF

    IF (problem_type .ne. -1) THEN
        CALL build_delaunay_triangulation(ve, pt, shuffle, verbose=.true.)
    END IF

! --- Output ---

    ! get the number of elements
    nelem = size(ve,2)
    ! get the number of nodes
    nnode = size(pt,2)
    ! Compute memory distance between elements
    CALL delaunay_compute_dual_grid(dual, dualb, ve, nnode)
    CALL delaunay_compute_neighbors(ne, ve, dual, dualb)
    CALL compute_unique_edges(edge, ve, ne)
    ! Compute voronoi control volumes
    CALL build_polygon_mesh(vor_pt, vor_ve, vor_veb, pt, ve, mesh_type=voronoi_mesh_type)
    CALL prt_bin(pt,      trim(output_path)//"/dt_pt.dat")
    CALL prt_bin(ve,      trim(output_path)//"/dt_ve.dat")
    CALL prt_bin(edge,    trim(output_path)//"/dt_edge.dat")
    CALL prt_bin(vor_pt,  trim(output_path)//"/dt_vor_pt.dat")
    CALL prt_bin(vor_ve,  trim(output_path)//"/dt_vor_ve.dat")
    CALL prt_bin(vor_veb, trim(output_path)//"/dt_vor_veb.dat")
    
CONTAINS

    SUBROUTINE parse_command_line(filename_config, using_default_config_file)
        CHARACTER(LEN=*), INTENT(OUT) :: filename_config
        LOGICAL,          INTENT(OUT) :: using_default_config_file
        INTEGER                       :: n_arg     
        n_arg = COMMAND_ARGUMENT_COUNT()
        using_default_config_file = .false.
        IF (n_arg .eq. 0) THEN
            WRITE(*,"(A)") "Using default config file smesh_test.cfg" 
            filename_config = "smesh_test.cfg"
            using_default_config_file = .true.
        ELSE IF (n_arg .eq. 1) THEN
            CALL get_command_argument(1, filename_config) 
        ELSE
            write(*,"(A)") "ERROR, syntax error in command line"
            STOP
        END IF
    END SUBROUTINE parse_command_line

    SUBROUTINE mesh_basic(points, xmin, xmax, ymin, ymax, nx, ny)
        REAL(8),                              INTENT(IN)  :: xmin, xmax, ymin, ymax
        INTEGER,                              INTENT(IN)  :: nx, ny
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: points
        REAL(8)                                           :: dx, dy
        INTEGER                                           :: i, j, k
        ALLOCATE(points(2, nx*ny + (ny - mod(ny, 2))/2))
        dx = real(xmax - xmin, 8)/real(nx - 1, 8)
        dy = real(ymax - ymin, 8)/real(ny - 1, 8)
        k = 1
        DO j = 0,ny-1
            DO i = 0,nx-1
                points(:,k) = [xmin + real(i, 8)*dx, ymin + real(j, 8)*dy]
                IF (mod(j, 2) == 1 .and. i /= 0) THEN
                    points(1,k) = points(1,k) - 0.5d0*dx
                END IF
                k = k + 1
                IF (mod(j, 2) == 1 .and. i == nx - 1) THEN
                    points(1,k) = points(1,k-1) + 0.5d0*dx
                    points(2,k) = points(2,k-1) 
                    k = k + 1
                END IF
            END DO
        END DO
    END SUBROUTINE mesh_basic

    SUBROUTINE mesh_flower(points, dphi, drho_k, drho_m, n)
        ! 0.6180339887498948482 ! dphi
        ! 1.0                   ! drho_k
        ! 0.5                   ! drho_m
        ! 5000                  ! n
        REAL(8),                              INTENT(IN)  :: dphi, drho_k, drho_m
        INTEGER,                              INTENT(IN)  :: n
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: points
        INTEGER                                           :: i
        REAL(8)                                           :: rho, phi
        REAL(8), PARAMETER                                :: twopi = 8.0d0*atan(1.0d0)
        ALLOCATE(points(2, n))
        DO i = 1, n
            phi = real(i-1, 8)*dphi*twopi
            rho = drho_k*(real(i, 8)*dphi)**drho_m
            points(:,i) = [rho*cos(phi), rho*sin(phi)]
        END DO
    END SUBROUTINE mesh_flower

    SUBROUTINE mesh_bisected_rectangle(pt, ve, xmin, xmax, ymin, ymax, nx, ny)
        REAL(8),                              INTENT(IN)  :: xmin, xmax, ymin, ymax
        INTEGER,                              INTENT(IN)  :: nx, ny
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: pt
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ve
        REAL(8)                                           :: dx, dy
        INTEGER                                           :: i, j, k
        ALLOCATE(pt(2, (nx+1)*(ny+1)))
        ALLOCATE(ve(3, 2*nx*ny))
        dx = (xmax - xmin)/real(nx, 8)
        dy = (ymax - ymin)/real(ny, 8)
        DO j = 0, ny
            DO i = 0, nx
                k = j*(nx+1) + i + 1
                pt(:,k) = [xmin + real(i, 8)*dx, ymin + real(j, 8)*dy]
            END DO
        END DO
        k = 0
        DO j = 0, ny - 1
            DO i = 0, nx - 1
                k = k + 1
                ve(:,k) = [j*(nx+1) + i + 1, j*(nx+1) + i + 2, (j+1)*(nx+1) + i + 1]
                k = k + 1
                ve(:,k) = [j*(nx+1) + i + 2, (j+1)*(nx+1) + i + 2, (j+1)*(nx+1) + i + 1]
            END DO
        END DO
    END SUBROUTINE mesh_bisected_rectangle

    ! PURE SUBROUTINE mesh_circle(pt, xmin, xmax, ymin, ymax, IMAX, JMAX)
    !     REAL(8),                              INTENT(IN)  :: xmin, xmax, ymin, ymax
    !     INTEGER,                              INTENT(IN)  :: IMAX, JMAX
    !     REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: pt
    !     REAL(8)                                           :: drho, dtheta
    !     INTEGER                                           :: i, j, k
    !     drho = (xmax - xmin)/real(IMAX, 8)
    !     dtheta = (ymax - ymin)/real(JMAX, 8)
    !     ALLOCATE(pt(2,(IMAX+1)*(JMAX+1)))
    !     k = 0
    !     DO i = 1, IMAX
    !         DO j = 1, JMAX
    !             k = k + 1
    !         END DO
    !     END DO

    ! END SUBROUTINE mesh_circle

! ---- settings for problem type 10: flower ---- !


! def square_point_list():
!     x = np.linspace(-1, 1, num=11, endpoint=True)
!     y = np.linspace(-1, 1, num=11, endpoint=True)
!     input_points = np.zeros((len(x)*len(y), 2))
!     k = 0
!     dx = x[1]-x[0]
!     dy = y[1]-y[0]
!     for j in range(len(y)):
!         for i in range(len(x)):
!             input_points[k,0] = x[0] + (i-1)*dx
!             input_points[k,1] = y[0] + (j-1)*dy
!             k += 1
!     return input_points


! def triangle_point_list():
!     x = np.linspace(-1, 1, num=11, endpoint=True)
!     y = np.linspace(-1, 1, num=11, endpoint=True)
!     input_points = np.zeros((len(x)*len(y) + (len(y)-len(y)%2)//2, 2))
!     k = 0
!     dx = x[1]-x[0]
!     dy = y[1]-y[0]
!     for j in range(len(y)):
!         for i in range(len(x)):
!             input_points[k,0] = x[0] + (i-1)*dx
!             input_points[k,1] = y[0] + (j-1)*dy
!             if j%2 == 1 and i != 0:
!                 input_points[k,0] = input_points[k,0] - 0.5*dx
!             k += 1
!             if j%2 == 1 and i == len(x)-1: # warning python
!                 input_points[k,0] = input_points[k-1,0] + 0.5*dx
!                 input_points[k,1] = input_points[k-1,1]
!                 k += 1
!     return input_points


! def circle_point_list():
!     points = np.zeros((31*31,2))
!     k = 0
!     for i in range(31):
!         for j in range(31):
!             r = np.random.random()
!             theta = 2*np.pi*np.random.random()
!             points[k,:] = (r*np.cos(theta), r*np.sin(theta))
!             k += 1
!     return points

END PROGRAM smesh_run



