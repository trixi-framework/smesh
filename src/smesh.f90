! Compute the Delaunay triangulation of a given set of points in 2d.
!
!     The subroutine is build_delaunay_triangulation(ve, input_points, shuffle, verbose)
!
!     Input:
!         input_points: a 2 X npoints matrix
!         shuffle: .true. or .false. (optional, defaults to .false.)
!         verbose: .true. or .false. (optional, defaults to .false.)
!
!     Output:
!         ve: an allocatable matrix (the triangulation matrix)
!

MODULE smesh

    ! USE smesh_io, only: prt_bin, prt

    IMPLICIT NONE

    PRIVATE

    ! correct calling convention CALL subroutine_name(outputs, inputs)
    PUBLIC :: build_delaunay_triangulation
    PUBLIC :: duplicate_cleanup
    PUBLIC :: build_voronoi_control_volumes
    PUBLIC :: physical_space_to_reference_space_build_map
    PUBLIC :: physical_space_to_reference_space
    PUBLIC :: reference_space_to_physical_space
    PUBLIC :: compute_neighbors_delaunay
    PUBLIC :: sort_dual_grid
    PUBLIC :: compute_dual_grid
    PUBLIC :: compute_unique_edges
    PUBLIC :: compute_unique_edges_voronoi

    ! PUBLIC :: segments_intersect
    ! PUBLIC :: locate_vertex
    ! PUBLIC :: circumcenter

    ! old calling convention
    ! PUBLIC :: onedge
    ! PUBLIC :: intriangle

CONTAINS

    ! ---- The delaunay triangulation routine ---------------------------------------------------- !

    SUBROUTINE build_delaunay_triangulation(ve_out, data_points, shuffle, verbose)
        ! Please do not input overlapping points. 
        !     You can use duplicate_cleanup to remove the duplicates beforehand
        !     
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ve_out
        REAL(8), DIMENSION(:,:),              INTENT(IN)  :: data_points
        LOGICAL, OPTIONAL,                    INTENT(IN)  :: shuffle
        LOGICAL, OPTIONAL,                    INTENT(IN)  :: verbose
        INTEGER, DIMENSION(:,:), ALLOCATABLE              :: ve ! Vertices
        INTEGER, DIMENSION(:,:), ALLOCATABLE              :: ne ! Neighbors
        INTEGER, DIMENSION(:,:), ALLOCATABLE              :: ov ! Opposing vertices in neighbors
        REAL(8), DIMENSION(:,:), ALLOCATABLE              :: pt ! Points
        REAL(8), DIMENSION(:,:), ALLOCATABLE              :: bc ! Centroids
        INTEGER, DIMENSION(:,:), ALLOCATABLE              :: edge_stack ! stack for edges
        INTEGER, DIMENSION(:),   ALLOCATABLE              :: insertion_order, search_steps
        LOGICAL, DIMENSION(:),   ALLOCATABLE              :: copy_triangle
        INTEGER                                           :: npmax, ntmax, nemax
        INTEGER                                           :: count_flips, edge_stack_i, k_out
        INTEGER                                           :: i, j, k, k1, k2, ii, jj, i0, i1, i2, i3
        INTEGER                                           :: v1, v2, v3, v4
        INTEGER                                           :: j1, j2, j3, m2, m3
        INTEGER                                           :: n1, n2, n3, ae
        LOGICAL                                           :: triangle_inside, triangle_edge, be_verbose
        INTEGER                                           :: search_index
        REAL(8), DIMENSION(2)                             :: bbox_center 
        REAL(8)                                           :: bbox_scale
        REAL(8), PARAMETER                                :: oot = 1.0d0/3.0d0 ! one over three
        INTEGER, PARAMETER, DIMENSION(2,3)                :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)
        REAL(8), DIMENSION(:),   ALLOCATABLE              :: search_sample_point_distance
        INTEGER, DIMENSION(:),   ALLOCATABLE              :: search_sample_point_triangle
        REAL(8), PARAMETER                                :: search_sample_size_exponent = 1.0d0/3.0d0
        REAL(8), PARAMETER                                :: search_sample_size_scale    = 2.0d0
        INTEGER                                           :: search_sample_size
        REAL(8)                                           :: search_sample_step

        IF (present(verbose)) THEN
            be_verbose = verbose
        ELSE
            be_verbose = .true.
        END IF

        j2 = 0

        ! Estimate the dimensions of the dataset
        npmax = 4 + size(data_points,2) ! The number of points provided plus a bounding triangle
        ntmax = 2*npmax - 2 - 4 ! We expect ntmax = 2*npmax - 2 - npmax_convex_hull
        nemax = 3*npmax - 3 - 4 ! We expect nemax = 3*npmax - 3 - npmax_convex_hull
        
        ! allocate sampling arrays for point location
        search_sample_size = min(ceiling(search_sample_size_scale*real(ntmax, 8)**search_sample_size_exponent), ntmax)
        ALLOCATE(search_sample_point_distance(search_sample_size + 1))
        ALLOCATE(search_sample_point_triangle(search_sample_size + 1))
        
        IF (be_verbose) THEN
            WRITE(*,"(A)") "Computing Delaunay triangulation."
            ! WRITE(*,"(A, I0, A, I0, A, I0)") "Computing Delaunay triangulation. npmax ", &
            !     npmax, ", ntmax ", ntmax, ", nemax ", nemax
        END IF

        ! Allocate the triangulation tables
        ALLOCATE(ve(3, ntmax), ne(3, ntmax), ov(3, ntmax))
        ALLOCATE(edge_stack(2, nemax*2))
        ALLOCATE(pt(2, npmax), bc(2, ntmax))
        ve = 0 ! Vertices
        ne = 0 ! Neighbors
        ov = 0 ! Opposing vertices
        bc = 0.0d0
        ! Initialize an edge stack, which stores 
        !     the triangle name and the pivot location local index
        edge_stack = 0

        ! Build the bounding box with two triangles
        CALL setup_computational_domain(data_points, pt, ve, ne, ov, bbox_center, bbox_scale)
        bc(:,1) = sum(pt(:,ve(:,1)),2)*oot
        bc(:,2) = sum(pt(:,ve(:,2)),2)*oot

        ! Choose the insertion order of the points
        ALLOCATE(insertion_order(npmax), search_steps(npmax))
        search_steps = 0
        DO i = 5,npmax
            insertion_order(i) = i
        END DO

        ! One can decide to pre-shuffle the points or not; if the input points have 
        !     some degree of spatial organization, apparently not shuffling them leads to
        !     very short search times when starting the search at the last inserted triangle, 
        !     while shuffling renders very few edge flips necessary in orger to 
        !     maintain the delaunay property. 

        IF (present(shuffle)) THEN
            IF (shuffle) THEN
                ! Fisher–Yates-Durstenfeld shuffle
                DO i = npmax, 2+4, -1
                    j = not_quite_random_integer(5, i, i)
                    v1 = insertion_order(i)
                    insertion_order(i) = insertion_order(j)
                    insertion_order(j) = v1
                END DO
            END IF 
        END IF

        ! Set the index of the last triangle in the list, 
        !     in k .eq. 1 and k .eq. 2 at first are stored the bounding box triangles
        k = 2 
        count_flips = 0 ! A counter for the edges
        edge_stack_i = 0 ! The index of the edge on top of the stack
        DO i0 = 5, npmax
            i = insertion_order(i0)
            ! Loop over all triangles built up until now 
            !     to find the one containing the current point i
            ! j = k ! begin with the last inserted triangle (not anymore: use random sampling instead)
            search_sample_size = min(ceiling(search_sample_size_scale*real(k, 8)**search_sample_size_exponent), k)
            search_sample_step = real(max(1, k-1), 8)/real(max(1, search_sample_size-1), 8)
            DO ii = 1, search_sample_size
                jj = nint(real(ii-1, 8)*search_sample_step + 1.0d0)
                search_sample_point_triangle(ii) = jj
                search_sample_point_distance(ii) = sum((bc(:,jj) - pt(:,i))**2)
            END DO
            j = search_sample_point_triangle(&
                minloc(search_sample_point_distance(:search_sample_size), 1))
            DO search_index = 1, k
                ! first we perform the edge test since we have seen that
                !     the triangle inside test validates some edges as well
                CALL onedge(pt(:,i), pt(:,ve(1,j)), pt(:,ve(2,j)), &
                    pt(:,ve(3,j)), triangle_edge, i1)
                IF (triangle_edge) THEN
                    ! The point has been located on an edge of triangle j
                    k1 = k + 1
                    k2 = k + 2
                    ! Set a frame of reference
                    i2 = ed(1,i1) ! mod(i1,3)+1
                    i3 = ed(2,i1) ! mod(i1+1,3)+1
                    ! Set names
                    v1 = ve(i1,j)
                    v2 = ve(i2,j)
                    v3 = ve(i3,j)
                    n1 = ne(i1,j)
                    n2 = ne(i2,j)
                    n3 = ne(i3,j)
                    ! If neighbors exist, update them and their neighbors
                    IF (n2 .ne. 0) THEN
                        ov(locate_vertex(ve(:,n2), ov(i2,j)),n2) = i
                    END IF
                    IF (n3 .ne. 0) THEN
                        ii = locate_vertex(ve(:,n3), ov(i3,j))
                        ov(ii,n3) = i
                        ne(ii,n3) = k1
                        ov( 3,k1) = ve(ii,n3)
                    END IF
                    ! Build k1
                    ve(:,k1) = [v1, v2, i]
                    ne(:,k1) = [n1, j, n3]
                    ov(2,k1) = v3
                    ! After, resize j 
                    ve(i2,j) = i
                    ne(i3,j) = k1
                    ov(i3,j) = v2
                    ! If n1 exists:
                    IF (n1 .ne. 0) THEN
                        ! Set up a frame of reference on the neighbor
                        j1 = locate_vertex(ve(:,n1), ov(i1,j))
                        j2 = ed(1,j1) ! mod(j1,3)+1
                        j3 = ed(2,j1) ! mod(j1+1,3)+1
                        v4 = ve(j1,n1)
                        m2 = ne(j2,n1)
                        m3 = ne(j3,n1)
                        ! Update k1
                        ov(1,k1) = v4
                        ! Update j
                        ne(i1,j) = k2
                        ov(i1,j) = v4
                        ! Build k2:
                        ve(:,k2) = [v4, v3,  i]
                        ne(:,k2) = [ j, n1, m3]
                        ov(1,k2) = v1
                        ov(2,k2) = v2
                        ! If neighbors exist, update them and their neighbors
                        IF (m2 .ne. 0) THEN
                            ov(locate_vertex(ve(:,m2), ov(j2,n1)),m2) = i
                        END IF
                        IF (m3 .ne. 0) THEN
                            ii = locate_vertex(ve(:,m3), ov(j3,n1))
                            ov(ii,m3) = i
                            ne(ii,m3) = k2
                            ov( 3,k2) = ve(ii,m3)
                        END IF
                        ! Resize the neighbor n1
                        ve(j2,n1) = i
                        ne(j1,n1) = k1
                        ne(j3,n1) = k2
                        ov(j3,n1) = v3
                    END IF
                    ! Update Centroids
                    bc(:, j) = sum(pt(:,ve(:, j)),2)*oot
                    bc(:,k1) = sum(pt(:,ve(:,k1)),2)*oot
                    ! Push the new edges on the stack and set the top edge index      
                    IF (n1 .ne. 0) THEN
                        !  WRITE(*,"(A, I6, A, I0, A, I6, A, I6, A, I6, A, I6, A, I6)") "Point ", &
                        !      i, " located on edge ", i1, " of triangle ", j, &
                        !      " | Building triangles ", j, ", ", n1, ", ", k1, ", ", k2
                        edge_stack(:,1) = [ j, i2]
                        edge_stack(:,2) = [k1,  3]
                        edge_stack(:,3) = [n1, j2]
                        edge_stack(:,4) = [k2,  3]
                        edge_stack_i = 4
                        ! Update Centroids
                        bc(:,k2) = sum(pt(:,ve(:,k2)), 2)*oot
                        bc(:,n1) = sum(pt(:,ve(:,n1)), 2)*oot
                        k = k2 ! 4 triangles added - 2 deleted
                    ELSE
                        !  WRITE(*,"(A, I6, A, I0, A, I6, A, I6, A, I6)") "Point ", i, &
                        !      " located on edge ", i1, " of triangle ", j, &
                        !      " | Building triangles ", j, ", ", k1
                        edge_stack(:,1) = [ j, i2]
                        edge_stack(:,2) = [k1,  3]
                        edge_stack_i = 2
                        k = k1 ! 2 triangles added - 1 deleted
                    END IF
                ELSE
                    ! Check if the point lies inside or on the edge of the current triangle
                    CALL intriangle(pt(:,i), pt(:,ve(1,j)), &
                        pt(:,ve(2,j)), pt(:,ve(3,j)), triangle_inside)
                    IF (triangle_inside) THEN
                        ! The point has been located inside the triangle j
                        ! Save the information about the old triangle j and build k1 and k2
                        k1 = k + 1
                        k2 = k + 2
                        ! WRITE(*,"(A, I6, A, I6, A, I6, A, I6, A, I6)") "Point ", i, &
                        !     " located in triangle ", j, " | Building triangles ", j, &
                        !     ", ", k1, ", ", k2
                        v1 = ve(1,j)
                        v2 = ve(2,j)
                        v3 = ve(3,j)
                        n1 = ne(1,j)
                        n2 = ne(2,j)
                        n3 = ne(3,j)
                        ! Store [T1 in j], [T2 in k1], [T3 in k2]
                        ve(3, j) = i ! Only the last of ve(j,:) needs to be updated
                        !          [v1, v2, i]
                        ve(:,k1) = [v2, v3, i]
                        ve(:,k2) = [v3, v1, i]
                        ! Neighbors of the new triangles
                        ne(1:2,j) = [k1, k2] ! T2, T3, N3 ! Only the first two need to be updated
                        !           [k1, k2, n3]
                        ne(:, k1) = [k2,  j, n1] ! T3, T1, N1
                        ne(:, k2) = [ j, k1, n2] ! T1, T2, N2
                        ! If neighbors exists, update them and their neighbors
                        IF (n3 .ne. 0) THEN
                            ii = locate_vertex(ve(:,n3), ov(3,j))
                            ov( 3, j) = ve(ii,n3)
                            ne(ii,n3) = j
                            ov(ii,n3) = i
                        END IF
                        IF (n1 .ne. 0) THEN
                            ii = locate_vertex(ve(:,n1), ov(1,j))
                            ov( 3,k1) = ve(ii,n1)
                            ne(ii,n1) = k1
                            ov(ii,n1) = i
                        END IF
                        IF (n2 .ne. 0) THEN
                            ii = locate_vertex(ve(:,n2), ov(2,j))
                            ov( 3,k2) = ve(ii,n2)
                            ne(ii,n2) = k2
                            ov(ii,n2) = i
                        END IF
                        ! Opposing vertices to the new triangles
                        ov(1:2, j) = [v3, v3]
                        ov(1:2,k1) = [v1, v1]
                        ov(1:2,k2) = [v2, v2]
                        ! Push the new edges on the stack
                        edge_stack(:,1) = [ j,3]
                        edge_stack(:,2) = [k1,3]
                        edge_stack(:,3) = [k2,3]
                        edge_stack_i = 3
                        ! Update Centroids
                        bc(:, j) = sum(pt(:,ve(:, j)),2)*oot
                        bc(:,k1) = sum(pt(:,ve(:,k1)),2)*oot
                        bc(:,k2) = sum(pt(:,ve(:,k2)),2)*oot
                        k = k2 ! 3 triangles added - 1 deleted
                    END IF
                END IF
                IF (triangle_inside .or. triangle_edge) THEN
                    ! If any location for the point has been found,
                    ! enforce the delaunay property of the triangulation
                    DO WHILE (edge_stack_i .gt. 0)
                        ! If point i has been found inside a triangle or on an edge,
                        !    until the stack is empty, check the active triangle/edge
                        !    and flip if non locally Delaunay
                        ! Copy triangle and pivot index as ae, i1
                        ae = edge_stack(1, edge_stack_i) ! The active element (triangle)
                        i1 = edge_stack(2, edge_stack_i) ! The pivot index
                        i2 = ed(1,i1) ! mod(i1, 3) + 1 ! The second counterclockwise index
                        i3 = ed(2,i1) ! mod(i1+1, 3) + 1 ! The third counterclockwise index
                        n1 = ne(i1,ae) ! The neighbor name associated to the pivot index
                        edge_stack_i = edge_stack_i - 1 ! Pop an edge from the stack
                        ! The first column of the edge stack contains the main triangle name,
                        !     we access the vertex ov[ae,i1] from the index in the second column
                        IF (n1 .ne. 0) THEN ! If the neighbor exists we are not opposing a border
                            ! WRITE(*,"(A, I5)") "    checking triangle ", ae
                            IF (not_delaunay(pt(:,ve(1,ae)), pt(:,ve(2,ae)), pt(:,ve(3,ae)), &
                                pt(:,ov(i1,ae)))) THEN
                                ! WRITE(*,"(A, I5, A, I0, A, I0, A, I0)") " --> flip triangle ", &
                                !     ae, ", ov ", i1, ", vertex ", ov(ae,i1), " on triangle ", n1
                                count_flips = count_flips + 1
                                ! Set auxiliary names
                                n2 = ne(i2,ae)
                                n3 = ne(i3,ae)
                                v1 = ve(i1,ae)
                                v2 = ve(i2,ae)
                                v3 = ve(i3,ae)
                                v4 = ov(i1,ae)
                                ! Set up a frame of reference on the neighbor
                                j1 = locate_vertex(ve(:,n1), ov(i1,ae))
                                j2 = ed(1,j1) ! mod(j1,3)+1
                                j3 = ed(2,j1) ! mod(j1+1,3)+1
                                m2 = ne(j2,n1)
                                m3 = ne(j3,n1)
                                ! If neighbors exists, update them and their neighbors,
                                !     and push new edges on the stack
                                IF (n3 .ne. 0) THEN
                                    ii = locate_vertex(ve(:,n3), ov(i3,ae))
                                    ov(ii,n3) = v4
                                    edge_stack_i = edge_stack_i + 1
                                    edge_stack(:, edge_stack_i) = [ae, i3]
                                END IF
                                IF (m3 .ne. 0) THEN
                                    ii = locate_vertex(ve(:,m3), ov(j3,n1))
                                    ov(ii,m3) = v1
                                    edge_stack_i = edge_stack_i + 1
                                    edge_stack(:, edge_stack_i) = [n1, j3]
                                END IF
                                IF (n2 .ne. 0) THEN
                                    ii = locate_vertex(ve(:,n2), ov(i2,ae))
                                    ne(ii,n2) = n1
                                    ov(ii,n2) = v4
                                    ov(j1,n1) = ve(ii,n2)
                                    edge_stack_i = edge_stack_i + 1
                                    edge_stack(:, edge_stack_i) = [n1, j1]
                                END IF
                                IF (m2 .ne. 0) THEN
                                    ii = locate_vertex(ve(:,m2), ov(j2,n1))
                                    ne(ii,m2) = ae
                                    ov(ii,m2) = v1
                                    ov(i1,ae) = ve(ii,m2)
                                    edge_stack_i = edge_stack_i + 1
                                    edge_stack(:, edge_stack_i) = [ae, i1]
                                END IF
                                ! Rearrange triangles
                                ve(i3,ae) = v4
                                ve(j3,n1) = v1
                                ne(i1,ae) = m2
                                ne(j1,n1) = n2
                                ne(i2,ae) = n1
                                ne(j2,n1) = ae
                                ov(i2,ae) = v3
                                ov(j2,n1) = v2
                                ! Update Centroids
                                bc(:,ae) = sum(pt(:,ve(:,ae)),2)*oot
                                bc(:,n1) = sum(pt(:,ve(:,n1)),2)*oot
                            END IF
                        END IF
                    END DO
                    search_steps(i) = search_index
                    EXIT ! We have found a triangle or edge, do not scan further
                ELSE
                    ! get a new triangle j to be tested
                    ! loop over each edge
                    DO ii = 1, 3
                        ! check if edge intersects the line connecting the 
                        !     centroid of the triangle and the destination point
                        IF (segments_intersect(pt(:,ve(ed(1,ii),j)), pt(:,ve(ed(2,ii),j)), &
                            bc(:,j), pt(:,i))) THEN
                            j = ne(ii, j)
                            EXIT
                        END IF
                    END DO
                END IF
            END DO
        END DO

        ! Write the triangulation on the output matrix
        k_out = 0
        ALLOCATE(copy_triangle(k))
        DO i = 1, k ! k being the index of the last inserted triangle
            copy_triangle(i) = .true.
            DO ii = 1, 3
                DO i0 = 1, 4
                    copy_triangle(i) = copy_triangle(i) .and. ve(ii,i) .ne. i0
                END DO
            END DO
            IF (copy_triangle(i)) THEN
                k_out = k_out + 1
            END IF
        END DO
        ALLOCATE(ve_out(3,k_out))
        i = 0
        DO j = 1, k ! k being the index of the last inserted triangle
            IF (copy_triangle(j)) THEN
                i = i + 1
                ve_out(:,i) = ve(:,j) - 4
            END IF
        END DO
        IF (be_verbose) THEN
            ! WRITE(*,"(A, I0, A, I0, A, I0)") "Delaunay triangulation: npmax ", npmax, &
            !     ", ntmax ", ntmax, ", nemax ", nemax
            WRITE(*,"(A, I10)") "Triangulation elements: ", k_out
            WRITE(*,"(A, I10)") "Total flipped edges:    ", count_flips
            WRITE(*,"(A, F12.2)") "Average search time:  ", REAL(sum(search_steps),8)/REAL(npmax-4,8)
            WRITE(*,"(A, F12.2)") "Flips/triangle:       ", REAL(count_flips,8)/REAL(k_out,8)
            WRITE(*,"(A, F12.2)") "Flips/node:           ", REAL(count_flips,8)/REAL(npmax,8)
        END IF
    END SUBROUTINE build_delaunay_triangulation

    SUBROUTINE duplicate_cleanup(points, segments)
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: points_clean
        INTEGER, DIMENSION(:), ALLOCATABLE :: replace_map
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: points
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: segments
        REAL(8), PARAMETER :: tol = 1.0d-14**2
        INTEGER :: i, j, k, np, ns, n_replaced, c, v1
        REAL(8) :: sqd
        LOGICAL :: remove_point

        np = size(points, 2)
        ns = size(segments, 2)
        n_replaced = 0

        ALLOCATE(points_clean(2,np))
        ALLOCATE(replace_map(np))
        replace_map = 0
        DO i = 1, np
            remove_point = .false.
            DO j = 1, np
                IF (i .ne. j .and. replace_map(j) .eq. 0 .and. replace_map(i) .eq. 0) THEN
                    sqd = sum((points(:,i) - points(:,j))**2)
                    IF (sqd .lt. tol) THEN
                        replace_map(j) = i
                        n_replaced = n_replaced + 1
                        write(*,"(A, I10, I10, 1P 2 E27.16)") "Removing duplicate nodes! ", i, j, points(:,i)
                        remove_point = .true.
                    END IF
                END IF
            END DO
        END DO
        points_clean = points
        DEALLOCATE(points)
        ALLOCATE(points(2,np-n_replaced))
        c = 0
        DO i = 1, np
            IF (replace_map(i) .eq. 0) THEN
                c = c + 1
                points(:,c) = points_clean(:,i)
                DO j = 1, ns
                    DO k = 1, 2
                        v1 = segments(k,j)
                        IF (v1 .eq. i .or. replace_map(v1) .eq. i) THEN
                            segments(k,j) = c
                        END IF
                    END DO
                END DO
            END IF
        END DO
    END SUBROUTINE duplicate_cleanup

    ! ---- Walking procedures for location of elements ------------------------------------------- !

    PURE FUNCTION locate_segment(pt, ve, ne, bc, v1, v2, i0, i00)
        ! returns the index of the triangle containing v1 and v2 as an edge or 
        !     zero if no such triangle is found
        REAL(8), DIMENSION(:,:), INTENT(IN) :: pt, bc
        INTEGER, DIMENSION(:,:), INTENT(IN) :: ve, ne
        INTEGER,                 INTENT(IN) :: v1, v2 ! input verticese
        INTEGER,                 INTENT(IN) :: i0, i00 ! first try indices
        INTEGER                             :: locate_segment ! output index
        INTEGER                             :: i, k, n, v3, v4
        LOGICAL                             :: found
        REAL(8), DIMENSION(2)               :: x
        REAL(8), DIMENSION(2)               :: p1, p2, pm
        INTEGER, PARAMETER, DIMENSION(2,3)  :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)
        x = 0.5d0*(pt(:,v1) + pt(:,v2))
        locate_segment = 0
        i = i0
        IF (i .eq. 0) THEN
            i = i00
        END IF
        found = .false.
        DO n = 1, size(bc,2)
            DO k = 1, 3
                v3 = ve(ed(1,k),i)
                v4 = ve(ed(2,k),i)
                IF ((v1 .eq. v3 .and. v2 .eq. v4) .or. (v1 .eq. v4 .and. v2 .eq. v3)) THEN
                    found = .true.
                    EXIT
                END IF
            END DO
            IF (found) THEN
                locate_segment = i
                EXIT
            ELSE
                DO k = 1, 3
                    ! check if edge intersects the line connecting the barycenter 
                    !     of the triangle and the destination point
                    p1 = pt(:,ve(ed(1,k),i))
                    p2 = pt(:,ve(ed(2,k),i))
                    pm = 0.49d0*(p1 + p2) + 0.02d0*bc(:,i)
                    IF (segments_intersect(p1, p2, pm, x)) THEN
                        i = ne(k,i)
                        EXIT
                    END IF
                END DO
                IF (i .eq. 0) THEN
                    EXIT
                END IF
            END IF
        END DO
        IF (locate_segment .eq. 0) THEN
            i = i00
            DO n = 1, size(bc,2)
                DO k = 1, 3
                    v3 = ve(ed(1,k),i)
                    v4 = ve(ed(2,k),i)
                    IF ((v1 .eq. v3 .and. v2 .eq. v4) .or. (v1 .eq. v4 .and. v2 .eq. v3)) THEN
                        found = .true.
                        EXIT
                    END IF
                END DO
                IF (found) THEN
                    locate_segment = i
                    EXIT
                ELSE
                    DO k = 1, 3
                        ! check if edge intersects the line connecting the barycenter of the triangle and the destination point
                        p1 = pt(:,ve(ed(1,k),i))
                        p2 = pt(:,ve(ed(2,k),i))
                        pm = 0.49d0*(p1 + p2) + 0.02d0*bc(:,i)
                        IF (segments_intersect(p1, p2, pm, x)) THEN
                            i = ne(k, i)
                            EXIT
                        END IF
                    END DO
                    IF (i .eq. 0) THEN
                        EXIT
                    END IF
                END IF
            END DO
        END IF
    END FUNCTION locate_segment

    PURE FUNCTION locate_point(pt, ve, ne, bc, x, i0, i00)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: pt, bc
        INTEGER, DIMENSION(:,:), INTENT(IN) :: ve, ne
        INTEGER,                 INTENT(IN) :: i0, i00 ! first try indices
        INTEGER                             :: i, k, n
        INTEGER                             :: locate_point
        LOGICAL                             :: triangle_edge, triangle_inside
        REAL(8), DIMENSION(2), INTENT(IN)   :: x
        REAL(8), DIMENSION(2)               :: p1, p2
        INTEGER, PARAMETER, DIMENSION(2,3)  :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)
        locate_point = 0
        i = i0
        IF (i .eq. 0) THEN
            i = i00
        END IF
        DO n = 1, size(bc,2)
            CALL onedge_light(x, pt(:,ve(1,i)), pt(:,ve(2,i)), pt(:,ve(3,i)), triangle_edge)
            CALL intriangle(x, pt(:,ve(1,i)), pt(:,ve(2,i)), pt(:,ve(3,i)), triangle_inside)
            IF (triangle_edge .or. triangle_inside) THEN
                locate_point = i
                EXIT
            ELSE
                DO k = 1, 3
                    ! check if edge intersects the line connecting the barycenter 
                    !     of the triangle and the destination point
                    p1 = pt(:,ve(ed(1,k),i))
                    p2 = pt(:,ve(ed(2,k),i))
                    IF (segments_intersect(p1, p2, 0.49d0*(p1 + p2) + 0.02d0*bc(:,i), x)) THEN
                        i = ne(k,i)
                        EXIT
                    END IF
                END DO
                IF (i .eq. 0) THEN
                    EXIT
                END IF
            END IF
        END DO
        IF (locate_point .eq. 0) THEN
            i = i00
            DO n = 1, size(bc,2)
                CALL onedge_light(x, pt(:,ve(1,i)), pt(:,ve(2,i)), pt(:,ve(3,i)), triangle_edge)
                CALL intriangle(x, pt(:,ve(1,i)), pt(:,ve(2,i)), pt(:,ve(3,i)), triangle_inside)
                IF (triangle_edge .or. triangle_inside) THEN
                    locate_point = i
                    EXIT
                ELSE
                    DO k = 1, 3
                        p1 = pt(:,ve(ed(1,k),i))
                        p2 = pt(:,ve(ed(2,k),i))
                        ! check if edge intersects the line connecting the barycenter of the triangle and the destination point
                        IF (segments_intersect(p1, p2, 0.49d0*(p1 + p2) + 0.02d0*bc(:,i), x)) THEN
                            i = ne(k, i)
                            EXIT
                        END IF
                    END DO
                    IF (i .eq. 0) THEN
                        EXIT
                    END IF
                END IF
            END DO
        END IF
    END FUNCTION locate_point

    ! ---- Initialization of reference space for the delaunay triangulation procedure ------------ !

    PURE SUBROUTINE setup_computational_domain(points, pt, ve, ne, ov, &
        bbox_center, bbox_scale, bx, by, reference_diameter, ntmax_estimate)
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: points 
        REAL(8), OPTIONAL,       INTENT(IN)  :: reference_diameter 
        REAL(8), DIMENSION(:,:), INTENT(OUT) :: pt
        INTEGER, DIMENSION(:,:), INTENT(OUT) :: ve, ne, ov
        REAL(8), DIMENSION(2),   INTENT(OUT) :: bbox_center
        REAL(8),                 INTENT(OUT) :: bbox_scale
        REAL(8), OPTIONAL,       INTENT(OUT) :: bx, by
        INTEGER, OPTIONAL,       INTENT(OUT) :: ntmax_estimate
        REAL(8)                              :: xmin, xmax, ymin, ymax, width, height, l
        REAL(8)                              :: element_l, element_a
        INTEGER                              :: i
        ! Build the bounding box with two triangles
        xmin = minval(points(1,:))
        xmax = maxval(points(1,:))
        ymin = minval(points(2,:))
        ymax = maxval(points(2,:))
        width = xmax - xmin
        height = ymax - ymin
        ! give a rough estimate of the number of triangles to be generated
        IF (present(reference_diameter)) THEN
            element_a = 3.0d0*sqrt(3.0d0)*reference_diameter**2/16.0d0
            element_l = sqrt(3.0d0)*reference_diameter/2.0d0
            ntmax_estimate = nint(width*height/element_a + 2.0d0*(width + height)/element_l)  
        END IF
        ! set the coordinates of the center of the unscaled bounding box
        bbox_center = [(xmin + xmax), (ymin + ymax)]/2.0d0
        ! Make a square bounding box
        l = max(width, height)
        ! scale to unity coordinates
        bbox_scale = 3.0d0*l/2.0d0
        ! Initialize a new point list
        DO i = 1, size(points,2)
            pt(:,i+4) = (points(:,i) - bbox_center)/bbox_scale
        END DO
        ! bx and by are the bounds beyond which there is no need to enforce quality
        ! maybe valid for convex boxes only?
        IF (present(bx)) THEN
            bx = 0.5d0*width/bbox_scale
        END IF
        IF (present(by)) THEN
            by = 0.5d0*height/bbox_scale
        END IF
        ! Set the vertices of the bounding triangles
        pt(:,1) = [-1.0d0, -1.0d0]
        pt(:,2) = [ 1.0d0, -1.0d0]
        pt(:,3) = [ 1.0d0,  1.0d0]
        pt(:,4) = [-1.0d0,  1.0d0]
        ve(:,1) = [1, 2, 3]
        ne(2,1) = 2
        ov(2,1) = 4
        ve(:,2) = [3, 4, 1]
        ne(2,2) = 1
        ov(2,2) = 2
    END SUBROUTINE setup_computational_domain

    ! ---- Basic geometrical predicates and constructors ----------------------------------------- !

    PURE FUNCTION not_quite_random_integer(minimum, maximum, seed)
        ! We generate a random integer for shuffle operations
        !     A seed is required so that the result exactly depends on the input
        ! Unfortunately we must use 64bit integers in order to compute
        !     a 32 bit integer (without exploiting modular representation)
        INTEGER, INTENT(IN)   :: minimum, maximum, seed
        INTEGER               :: not_quite_random_integer
        INTEGER(8), PARAMETER :: a = 1103515245_8, b = 12345_8
        INTEGER(8)            :: x
        x = mod(a*int(seed, 8) + b, 2147483647_8)
        not_quite_random_integer = minimum + mod(int(x, 4), maximum-minimum+1)
    END FUNCTION not_quite_random_integer

    PURE FUNCTION cross2d(u, v)
        ! Returns the third component of the cross product of two vertices lying on a plane
        REAL(8), DIMENSION(2), INTENT(IN) :: u, v
        REAL(8)                           :: cross2d
        cross2d = u(1)*v(2) - u(2)*v(1)
    END FUNCTION cross2d

    PURE FUNCTION segments_intersect(p1, p2, p3, p4)
        ! Used for point localization
        !     Returns .true. if the two segments p1--p2 and p3--p4 intersect
        !     (or even just come close to intersecting?)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3, p4
        REAL(8), DIMENSION(2)             :: v21, v43, v32, v42, v14, v24
        LOGICAL                           :: segments_intersect
        REAL(8), PARAMETER                :: DPmacheps = 2.22044604926d-16
        REAL(8), PARAMETER                :: epsilon = DPmacheps**2
        v21 = p2 - p1
        v43 = p4 - p3
        v32 = p3 - p2
        v42 = p4 - p2
        v14 = p1 - p4
        v24 = p2 - p4
        segments_intersect = cross2d(v21, v32)*cross2d(v21, v42) .lt. epsilon .and. &
            cross2d(v43, v14)*cross2d(v43, v24) .lt. epsilon
    END FUNCTION segments_intersect

    PURE FUNCTION onsegment(p, a, b)
        ! Tests if point p lies between point a and point b.
        !     The function is used by the routine for point localization at edges 
        !     and to detect edges to be constrained in ruppert's algorithm
        REAL(8), DIMENSION(2), INTENT(IN) :: p, a, b
        REAL(8), DIMENSION(2)             :: v1, v2
        LOGICAL                           :: onsegment
        REAL(8)                           :: v1v2
        REAL(8), PARAMETER                :: DPmacheps = 2.22044604926d-16
        REAL(8), PARAMETER                :: epsilon = 2.0d0*DPmacheps
        ! REAL(8), PARAMETER                :: epsilon = 1.0d-14
        v1 = p - a
        v2 = b - a
        v1v2 = sum(v1*v2)
        onsegment = v1v2 .ge. sum(v1**2) .and. &
            abs(v1(1)*v2(2) - v1(2)*v2(1)) .lt. max(1.0d-10*sqrt(sum(v2**2)), epsilon)
    END FUNCTION onsegment

    PURE SUBROUTINE fuzzy_test_segment(p, a, b, test, test_k)
        ! Tests if point p lies between point a and point b.
        !     The function is used by the routine for point localization at edges 
        !     and to detect edges to be constrained in ruppert's algorithm
        REAL(8), DIMENSION(2), INTENT(IN)  :: p, a, b
        REAL(8), DIMENSION(2)              :: v1, v2
        LOGICAL,               INTENT(OUT) :: test
        REAL(8),               INTENT(OUT) :: test_k
        REAL(8)                            :: v1v2
        REAL(8), PARAMETER                 :: DPmacheps = 2.22044604926d-16
        REAL(8), PARAMETER                 :: epsilon = 2.0d0*DPmacheps
        ! REAL(8), PARAMETER                 :: epsilon = 1.0d-14
        v1 = p - a
        v2 = b - a
        v1v2 = sum(v1*v2)
        test_k = abs(v1(1)*v2(2) - v1(2)*v2(1))
        test = v1v2 .ge. sum(v1**2) .and. test_k .lt. max(1.0d-9*sqrt(sum(v2**2)), epsilon)
    END SUBROUTINE fuzzy_test_segment

    PURE SUBROUTINE onedge(p, v1, v2, v3, edge, ii)
        ! Tests if p lies on an edge of triangle (v1, v2, v3).
        !     Returns the index of the edge yielding maximum collinearity
        REAL(8), DIMENSION(2), INTENT(IN)  :: p, v1, v2, v3
        LOGICAL,               INTENT(OUT) :: edge
        INTEGER,               INTENT(OUT) :: ii
        LOGICAL, DIMENSION(3)              :: test
        REAL(8), DIMENSION(3)              :: test_k
        CALL fuzzy_test_segment(p, v2, v3, test(1), test_k(1))
        CALL fuzzy_test_segment(p, v3, v1, test(2), test_k(2)) 
        CALL fuzzy_test_segment(p, v1, v2, test(3), test_k(3)) 
        edge = any(test)
        IF (edge) THEN
            ii = minloc(test_k, 1)
        ELSE
            ii = 1
        END IF
    END SUBROUTINE onedge

    PURE SUBROUTINE onedge_light(p, v1, v2, v3, edge)
        ! Tests if p lies on an edge of triangle (v1, v2, v3).
        REAL(8), DIMENSION(2), INTENT(IN)  :: p, v1, v2, v3
        LOGICAL,               INTENT(OUT) :: edge
        edge = onsegment(p, v2, v3) .or. onsegment(p, v3, v1) .or. onsegment(p, v1, v2)
    END SUBROUTINE onedge_light

    PURE SUBROUTINE intriangle(p, v1, v2, v3, inside)
        ! Tests if p is lies inside the triangle (v1, v2, v3).
        !     The vertices should be ordered counter-clockwise. 
        REAL(8), DIMENSION(2), INTENT(IN)  :: p, v1, v2, v3
        LOGICAL,               INTENT(OUT) :: inside
        REAL(8)                            :: dot33, dot32, dot30, dot22, dot20
        REAL(8)                            :: den, c1, c2
        REAL(8), DIMENSION(2)              :: u0, u2, u3
        REAL(8), PARAMETER                 :: tol = 1.0d0 + 1.0d-8
        ! Compute radius vectors with respect to v1
        u3 = v3 - v1
        u2 = v2 - v1
        u0 = p - v1
        ! Precompute the dot products of the radius vectors
        dot33 = sum(u3**2)
        dot22 = sum(u2**2)
        dot32 = sum(u3*u2)
        dot30 = sum(u3*u0)
        dot20 = sum(u2*u0)
        ! Precompute the denominator for the test
        !     The denominator should be always positive so absolute value should not not needed
        !     as division by den does not change the sign of c1 or c2
        den = dot33*dot22 - dot32**2
        ! Perform the intriangle test
        c1 = dot22*dot30 - dot32*dot20
        c2 = dot33*dot20 - dot32*dot30
        inside = c1 .gt. 0.0d0 .and. c2 .gt. 0.0d0 .and. c1 + c2 .lt. den*tol
    END SUBROUTINE intriangle
 
    PURE FUNCTION not_delaunay(v1, v2, v3, p)
        ! Returns .true. if triangle is NOT locally Delaunay.
        !     Assumes orientation of triangles to be counter-clockwise
        !     Strictly setting the threshold at zero instead of machine epsilon may
        !     unnecessarily flip some edges, but machine epsilon has shown to be too much
        !     in the strong gradation test case.
        REAL(8), DIMENSION(2), INTENT(IN) :: v1, v2, v3, p
        REAL(8), DIMENSION(2)             :: u, v, w
        LOGICAL                           :: not_delaunay
        u = v1 - p
        v = v2 - p
        w = v3 - p
        not_delaunay = sum(u**2)*(v(1)*w(2) - v(2)*w(1)) + &
                       sum(v**2)*(u(2)*w(1) - u(1)*w(2)) + &
                       sum(w**2)*(u(1)*v(2) - u(2)*v(1)) .gt. 0.0d0
    END FUNCTION not_delaunay

    PURE FUNCTION is_delaunay(v1, v2, v3, p)
        ! Returns .true. if triangle is locally Delaunay.
        REAL(8), DIMENSION(2), INTENT(IN) :: v1, v2, v3, p
        REAL(8), DIMENSION(2)             :: u, v, w
        LOGICAL                           :: is_delaunay
        u = v1 - p
        v = v2 - p
        w = v3 - p
        is_delaunay = sum(u**2)*(v(1)*w(2) - v(2)*w(1)) + &
                       sum(v**2)*(u(2)*w(1) - u(1)*w(2)) + &
                       sum(w**2)*(u(1)*v(2) - u(2)*v(1)) .le. 0.0d0
    END FUNCTION is_delaunay

    PURE FUNCTION locate_vertex(triangle, v)
        ! associates a local index in tringle such that it points to the global node v
        ! Returns the local vertex index such that the corresponding vertex in tv is ov
        ! tv = a vector containing the three global indices of the vertices of a triangle
        ! ov1 = the global index of the vertex
        INTEGER, DIMENSION(3), INTENT(IN) :: triangle 
        INTEGER,               INTENT(IN) :: v
        INTEGER                           :: locate_vertex
        DO locate_vertex = 1, 3
            IF (triangle(locate_vertex) .eq. v) THEN
                EXIT
            END IF
        END DO 
    END FUNCTION locate_vertex

    PURE FUNCTION circumcenter(p1, p2, p3)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8), DIMENSION(2)             :: circumcenter
        REAL(8), DIMENSION(2)             :: u, v
        REAL(8)                           :: uq, vq
        u = p1 - p3
        v = p2 - p3
        uq = sum(u**2)
        vq = sum(v**2)
        circumcenter = p3 + [uq*v(2) - vq*u(2), vq*u(1) - uq*v(1)]/(2.0d0*(u(1)*v(2) - u(2)*v(1)))
    END FUNCTION circumcenter

    PURE FUNCTION incenter(p1, p2, p3)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8), DIMENSION(2)             :: incenter
        REAL(8), DIMENSION(2)             :: u1, u2
        REAL(8), DIMENSION(3)             :: l
        u1 = p2 - p3
        u2 = p1 - p3
        l = [sqrt(sum(u1**2)), sqrt(sum(u2**2)), sqrt(sum((p1 - p2)**2))]
        incenter = p3 + (l(1)*u2 + l(2)*u1)/sum(l)
    END FUNCTION incenter

    PURE FUNCTION convcenter(p1, p2, p3)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8), DIMENSION(2)             :: convcenter
        REAL(8), DIMENSION(2)             :: u1, u2
        REAL(8), DIMENSION(3)             :: l
        u1 = p2 - p3
        u2 = p1 - p3
        l = [sqrt(sum(u1**2)), sqrt(sum(u2**2)), sqrt(sum((p1 - p2)**2))]
        convcenter = p3 + (l(1)*u2 + l(2)*u1)/sum(l)
    END FUNCTION convcenter

    PURE FUNCTION is_good_triangle(p1, p2, p3, q_min)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8),               INTENT(IN) :: q_min
        LOGICAL                           :: is_good_triangle
        REAL(8)                           :: u1, u2, v1, v2, uq, vq, wq
        u1 = p2(1) - p1(1)
        u2 = p2(2) - p1(2)
        v1 = p3(1) - p1(1)
        v2 = p3(2) - p1(2)
        uq = u1**2 + u2**2
        vq = v1**2 + v2**2
        wq = sum((p2 - p3)**2)
        is_good_triangle = 4.0d0*min(uq, vq, wq)*(u1*v2 - v1*u2)**2 .ge. q_min*((v2*uq - u2*vq)**2 + (u1*vq - v1*uq)**2)
    END FUNCTION is_good_triangle

    PURE FUNCTION is_bad_triangle(p1, p2, p3, q_min)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8),               INTENT(IN) :: q_min
        LOGICAL                           :: is_bad_triangle
        REAL(8)                           :: u1, u2, v1, v2, uq, vq, wq
        u1 = p2(1) - p1(1)
        u2 = p2(2) - p1(2)
        v1 = p3(1) - p1(1)
        v2 = p3(2) - p1(2)
        uq = u1**2 + u2**2
        vq = v1**2 + v2**2
        wq = sum((p2 - p3)**2)
        is_bad_triangle = 4.0d0*min(uq, vq, wq)*(u1*v2 - v1*u2)**2 .lt. q_min*((v2*uq - u2*vq)**2 + (u1*vq - v1*uq)**2)
    END FUNCTION is_bad_triangle

    PURE FUNCTION steiner_obtuse(p1, p2, p3)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8), DIMENSION(2)             :: u, v, w
        REAL(8)                           :: uq, vq, wq
        REAL(8), DIMENSION(2)             :: steiner_obtuse
        u = p3 - p2
        v = p1 - p3
        w = p2 - p1
        uq = sum(u**2)
        vq = sum(v**2)
        wq = sum(w**2)
        IF (uq .gt. vq .and. uq .gt. wq) THEN
            ! u
            steiner_obtuse = 0.5d0*(p2 + p3)
        ELSE IF (vq .gt. wq) THEN
            ! v
            steiner_obtuse = 0.5d0*(p3 + p1)
        ELSE
            ! w
            steiner_obtuse = 0.5d0*(p1 + p2)
        END IF
    END FUNCTION steiner_obtuse

    PURE FUNCTION obtuse_triangle(p1, p2, p3)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8), DIMENSION(2)             :: u, v, w
        LOGICAL                           :: obtuse_triangle
        u = p3 - p2
        v = p1 - p3
        w = p2 - p1
        obtuse_triangle = sum(u*v) .gt. 0.0d0 .or. sum(v*w) .gt. 0.0d0 .or. sum(w*u) .gt. 0.0d0
    END FUNCTION obtuse_triangle

    PURE FUNCTION obtuse(po, pe1, pe2)
        REAL(8), DIMENSION(2), INTENT(IN) :: po, pe1, pe2
        LOGICAL                           :: obtuse
        obtuse = sum((pe1 - po)*(pe2 - po)) .lt. 0.0d0
        ! in doubt we might say that pi/2 and slightly less is still
        ! obtuse and generate one more node, for robustness
        ! but maybe we introduce instability on bounding box borders?
    END FUNCTION obtuse

    PURE FUNCTION acute(po, pe1, pe2)
        REAL(8), DIMENSION(2), INTENT(IN) :: po, pe1, pe2
        LOGICAL                           :: acute
        acute = sum((pe1 - po)*(pe2 - po)) .gt. 0.0d0 
        ! in doubt we might say that pi/2 and slightly less is still
        ! obtuse and generate one more node, for robustness
        ! but maybe we introduce instability on bounding box borders?
    END FUNCTION acute

    PURE FUNCTION near_obtuse(po, pe1, pe2, cosmin)
        REAL(8), DIMENSION(2), INTENT(IN) :: po, pe1, pe2
        REAL(8),               INTENT(IN) :: cosmin
        LOGICAL                           :: near_obtuse
        near_obtuse = sum((pe1 - po)*(pe2 - po)) .lt. cosmin
    END FUNCTION near_obtuse

    PURE FUNCTION near_obtuse_triangle(p1, p2, p3, cosmin)
        REAL(8), DIMENSION(2), INTENT(IN) :: p1, p2, p3
        REAL(8),               INTENT(IN) :: cosmin
        REAL(8), DIMENSION(2)             :: u, v, w
        REAL(8)                           :: un, vn, wn
        LOGICAL                           :: near_obtuse_triangle
        u = p3 - p2
        v = p1 - p3
        w = p2 - p1
        un = sqrt(sum(u**2))
        vn = sqrt(sum(v**2))
        wn = sqrt(sum(w**2))
        near_obtuse_triangle = -sum(u*v) .lt. cosmin*un*vn .or. -sum(v*w) .lt. cosmin*vn*wn .or. &
            -sum(w*u) .lt. cosmin*wn*un
    END FUNCTION near_obtuse_triangle

    ! ---- Linear time operations on the entire mesh --------------------------------------------- !

    PURE SUBROUTINE compute_dual_grid(dual, dualb, ve, nnode)
        INTEGER, DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: dual
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: dualb 
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: ve
        INTEGER,                              INTENT(IN)  :: nnode
        INTEGER                                           :: i, j, nelem, ii, n
        nelem = size(ve, 2)
        IF (ALLOCATED(dual)) DEALLOCATE(dual) 
        IF (ALLOCATED(dualb)) DEALLOCATE(dualb) 
        ALLOCATE(dualb(2,nnode))
        dualb = 0
        ! Run over the triangles once in order to determine the size of the dual and dualn mtrices
        DO j = 1, nelem 
            DO ii = 1, 3
                i = ve(ii,j)
                dualb(2,i) = dualb(2,i) + 1
            END DO
        END DO
        dualb(1,1) = 1
        DO i = 2, nnode
            dualb(1,i) = dualb(1,i-1) + dualb(2,i-1) ! cumulative index
        END DO
        n = dualb(1,nnode) + dualb(2,nnode) - 1
        dualb(2,:) = dualb(1,:) - 1 ! set initially to the first position minus one for each node
        ALLOCATE(dual(n))
        ! worst case scenario we have that all the ndual_max triangles
        !     on a circle are disjointed which would mean we have to
        !     allocate for a number of vertices twice as big
        DO j = 1, nelem ! loop over all elements
            DO ii = 1, 3 ! loop over all vertices of an element
                ! get the global node number
                ! node_pointer(v1) gives the node number of v1 (old mesh) in the output mesh
                i = ve(ii,j)
                dualb(2,i) = dualb(2,i) + 1 ! count up one element
                dual(dualb(2,i)) = j
            END DO
        END DO
    END SUBROUTINE compute_dual_grid

    PURE SUBROUTINE sort_dual_grid(dual, dualb, ve, ne)
        ! Sorts the vertices of the dual grid dual(dualb(1,i),dualb(2,i)) in counter-clockwise order
        !     without crossing boundaries if any present
        INTEGER, DIMENSION(:,:), INTENT(IN)    :: ve
        INTEGER, DIMENSION(:,:), INTENT(IN)    :: ne
        INTEGER, DIMENSION(:),   INTENT(INOUT) :: dual
        INTEGER, DIMENSION(:,:), INTENT(IN)    :: dualb
        INTEGER, DIMENSION(:),   ALLOCATABLE   :: dual_column
        INTEGER                                :: i, j, k, nnode, ae, ae0, aep, u1, i1, i2, i3
        INTEGER                                :: dual_column_size
        LOGICAL                                :: counterclockwise, keep_spinning
        INTEGER, PARAMETER, DIMENSION(2,3)     :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)
        nnode = size(dualb, 2)
        dual_column_size = 3
        DO i = 1, nnode
            dual_column_size = max(dual_column_size, dualb(2,i) - dualb(1,i) + 1)
        END DO
        ALLOCATE(dual_column(dual_column_size))
        DO i = 1, nnode
            u1 = i
            ae = dual(dualb(2,i))
            ae0 = ae
            aep = ae
            k = 1
            dual_column = 0
            dual_column(k) = ae
            counterclockwise = .true.
            keep_spinning = .true.
            i1 = locate_vertex(ve(:,ae), u1) 
            DO WHILE (keep_spinning) ! spinning loop, as detailed above in generate mesh
                i2 = ed(1,i1)
                i3 = ed(2,i1)
                IF (counterclockwise) THEN
                    aep = ae
                    ae = ne(i2,ae)
                    IF (ae .eq. ae0) THEN
                        EXIT
                    ELSE
                        k = k + 1
                        dual_column(k) = ae
                        IF (ae .eq. 0) THEN
                            ! we crossed a boundary, start spinning clockwise from HERE
                            ae = ne(ed(2,locate_vertex(ve(:,aep), u1)),aep)
                            k = 1
                            dual_column(k) = aep
                            k = 2
                            dual_column(k) = ae
                            counterclockwise = .false.
                            IF (ae .eq. 0) THEN
                                EXIT
                            END IF
                        END IF
                        i1 = locate_vertex(ve(:,ae), u1) 
                        i2 = ed(1,i1)
                        i3 = ed(2,i1)
                    END IF
                ELSE
                    IF (ae .ne. 0) THEN
                        ae = ne(i3,ae)
                    ELSE
                        EXIT
                    END IF
                    IF (ae .eq. 0) THEN
                        EXIT
                    ELSE
                        k = k + 1
                        dual_column(k) = ae
                        i1 = locate_vertex(ve(:,ae), u1)
                        i2 = ed(1,i1)
                        i3 = ed(2,i1)
                    END IF
                END IF
            END DO
            IF (counterclockwise) THEN
                dual(dualb(1,i):dualb(2,i)) = dual_column(:dualb(2,i)-dualb(1,i)+1)
            ELSE
                k = dualb(2,i)-dualb(1,i)+1
                DO j = dualb(1,i), dualb(2,i)
                    dual(j) = dual_column(k)
                    k = k - 1
                END DO
            END IF
        END DO
    END SUBROUTINE sort_dual_grid

    PURE SUBROUTINE compute_neighbors_delaunay(ne, ve, dual, dualb)
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ne
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: ve
        INTEGER, DIMENSION(:),                INTENT(IN)  :: dual
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: dualb
        INTEGER                                           :: j, i1, j1, ii, ae, v2, v3, nelem
        INTEGER, PARAMETER, DIMENSION(2,3) :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)
        nelem = size(ve,2)
        IF(ALLOCATED(ne)) DEALLOCATE(ne)
        ALLOCATE(ne(3,nelem))
        ne = 0
        ! search for neighbors: loop over all triangles
        DO j = 1, nelem
            ! pivot on the i1_th vertex
            DO i1 = 1, 3
                ! get the edge opposing the pivot
                v2 = ve(ed(1,i1), j)
                v3 = ve(ed(2,i1), j)
                ! loop around one of the endpoints of the edge
                DO ii = dualb(1,v2), dualb(2,v2)
                    ae = dual(ii)
                    IF (ae .ne. j) THEN
                        DO j1 = 1, 3
                            IF (v2 .eq. ve(ed(2,j1),ae) .and. v3 .eq. ve(ed(1,j1),ae)) THEN
                                ! found a neighbor
                                ne(i1, j) = ae
                                EXIT
                            END IF
                        END DO
                        IF (ne(i1,j) .gt. 0) THEN
                            ! we found the neighbor opposing v1 = ve(i1,j)
                            EXIT
                        END IF
                    END IF
                END DO
            END DO
        END DO
    END SUBROUTINE compute_neighbors_delaunay

    PURE SUBROUTINE compute_unique_edges(edge, ve, ne)
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: ve
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: ne
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: edge
        LOGICAL, DIMENSION(:,:), ALLOCATABLE              :: edge_unique_marker
        INTEGER                                           :: j, nelem, i, n1, j1, i1, nedge
        INTEGER, PARAMETER, DIMENSION(2,3) :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)
        INTEGER, DIMENSION(2) :: edge1, edge2
        nelem = size(ve, 2)
        ALLOCATE(edge_unique_marker(3,nelem))
        edge_unique_marker = .false.
        DO j = 1, nelem
            DO i1 = 1, 3
                n1 = ne(i1,j)
                IF (n1 .ne. 0) THEN
                    edge1 = [ve(ed(1,i1),j), ve(ed(2,i1),j)]
                    edge1 = [minval(edge1), maxval(edge1)] 
                    DO j1 = 1, 3
                        edge2 = [ve(ed(1,j1),n1), ve(ed(2,j1),n1)]
                        edge2 = [minval(edge2), maxval(edge2)]
                        IF (all(edge1 .eq. edge2)) THEN
                            IF (.not. edge_unique_marker(j1,n1)) THEN
                                edge_unique_marker(i1,j) = .true.
                            END IF
                            EXIT
                        END IF
                    END DO
                ELSE
                    edge_unique_marker(i1,j) = .true.
                END IF
            END DO
        END DO
        nedge = 0
        DO j = 1, nelem
            DO i1 = 1, 3
                IF (edge_unique_marker(i1,j)) THEN
                    nedge = nedge + 1
                END IF
            END DO
        END DO
        IF (ALLOCATED(edge)) DEALLOCATE(edge)
        ALLOCATE(edge(2,nedge))
        i = 0
        DO j = 1, nelem
            DO i1 = 1, 3
                IF (edge_unique_marker(i1,j)) THEN
                    i = i + 1
                    edge(:,i) = [ve(ed(1,i1),j), ve(ed(2,i1),j)]
                END IF
            END DO
        END DO
        DEALLOCATE(edge_unique_marker)
    END SUBROUTINE compute_unique_edges

    PURE SUBROUTINE compute_unique_edges_voronoi(edge, edge_elements, ve, veb, ne)
        ! could be wrong, not tested
        INTEGER, DIMENSION(:),              INTENT(IN)  :: ve
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: veb
        INTEGER, DIMENSION(:),              INTENT(IN)  :: ne
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: edge
        INTEGER, DIMENSION(:), ALLOCATABLE              :: edge_found_marker
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)  :: edge_elements
        INTEGER, DIMENSION(:,:), ALLOCATABLE  :: edge_elements_temp
        INTEGER                                           :: nelem, i, n1, j1, i1, nedge, k
        INTEGER, PARAMETER, DIMENSION(2,3) :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)
        INTEGER, DIMENSION(2) :: edge1, edge2
        nelem = size(veb, 2)
        ALLOCATE(edge_found_marker(size(ve)))
        ALLOCATE(edge_elements_temp(2,size(ve)))
        edge_elements_temp = 0
        edge_found_marker = 0
        DO k = 1, nelem
            DO i1 = veb(1,k), veb(2,k)
                n1 = ne(i1)
                IF (n1 .ne. 0) THEN
                    IF (i1 .ne. veb(2,k)) THEN
                        edge1 = [ve(i1), ve(i1+1)]
                    ELSE
                        edge1 = [ve(i1), ve(veb(1,k))]
                    END IF
                    edge1 = [minval(edge1), maxval(edge1)] 
                    DO j1 = veb(1,n1), veb(2,n1)
                        IF (j1 .ne. veb(2,n1)) THEN
                            edge2 = [ve(j1), ve(j1+1)]
                        ELSE
                            edge2 = [ve(j1), ve(veb(1,n1))]
                        END IF
                        edge2 = [minval(edge2), maxval(edge2)] 
                        IF (all(edge1 .eq. edge2)) THEN
                            IF (edge_found_marker(i1) .eq. 0 .and. edge_found_marker(j1) .eq. 0) THEN
                                edge_found_marker(i1) = j1
                                edge_found_marker(j1) = -1
                                edge_elements_temp(:,i1) = [k, n1]
                            END IF
                        END IF
                    END DO
                ELSE
                    edge_found_marker(i1) = 1
                    edge_elements_temp(:,i1) = [k, -1]
                END IF
            END DO
        END DO
        nedge = 0
        DO i = 1, size(ve)
            IF (edge_found_marker(i) .gt. 0) THEN
                nedge = nedge + 1
            END IF
        END DO
        IF (ALLOCATED(edge)) DEALLOCATE(edge)
        ALLOCATE(edge(2,nedge))
        ALLOCATE(edge_elements(2,nedge))
        
        nedge = 0
        DO i = 1, size(ve)
            IF (edge_found_marker(i) .gt. 0) THEN
                nedge = nedge + 1
                edge(:,nedge) = [ve(i), ve(edge_found_marker(i))]
                edge_elements(:,nedge) = edge_elements_temp(:,i)
            END IF
        END DO

        DEALLOCATE(edge_found_marker)
        DEALLOCATE(edge_elements_temp)
    END SUBROUTINE compute_unique_edges_voronoi

    ! ---- Computation of Voronoi grid ----------------------------------------------------------- !

    PURE FUNCTION projection_on_segment(p0, p1, p2)
        ! projects p0 on the segment having p1 and p2 as endpoints
        REAL(8), DIMENSION(2), INTENT(IN) :: p0, p1, p2
        REAL(8), DIMENSION(2)             :: u0, u1
        REAL(8), DIMENSION(2)             :: projection_on_segment
        u0 = p0 - p2
        u1 = p1 - p2
        projection_on_segment = p2 + (sum(u0*u1)/sum(u1**2))*u1
    END FUNCTION projection_on_segment

    PURE SUBROUTINE build_voronoi_control_volumes(vor_pt, vor_ve, vor_veb, pt, ve, ne, mesh_type)
        ! mesh_type: 0-> voronoi if fully acute triangulation, 1-> centroid based polygons, 2-> funky incenter based polygons
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: vor_pt
        INTEGER, DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: vor_ve
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: vor_veb
        REAL(8), DIMENSION(:,:),              INTENT(IN)  :: pt
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: ve
        INTEGER, DIMENSION(:,:),              INTENT(IN)  :: ne
        INTEGER, OPTIONAL,                    INTENT(IN)  :: mesh_type
        INTEGER, DIMENSION(:,:), ALLOCATABLE              :: add_edge_midpoint
        INTEGER, DIMENSION(:),   ALLOCATABLE              :: add_vertex
        INTEGER, DIMENSION(:),   ALLOCATABLE              :: dual
        INTEGER, DIMENSION(:,:), ALLOCATABLE              :: dualb
        INTEGER                                           :: nelem, nnode, nnode_voronoi, nelem_voronoi
        INTEGER                                           :: i, j, ii, jj, kk, n, jl, jr, jn
        REAL(8), DIMENSION(2)                             :: p1, p2, p3
        REAL(8), PARAMETER                                :: oot = 1.0d0/3.0d0
        LOGICAL                                           :: inside, add_this_vertex
        INTEGER                                           :: mesh_type_selector
        INTEGER, PARAMETER, DIMENSION(2,3)                :: ed = reshape([2,3,3,1,1,2], [2,3]) ! mod(i+j, 3)

        nelem = size(ve, 2)
        nnode = size(pt, 2)
        CALL compute_dual_grid(dual, dualb, ve, nnode)
        CALL sort_dual_grid(dual, dualb, ve, ne) 
        nelem_voronoi = nnode
        IF (ALLOCATED(vor_veb)) DEALLOCATE(vor_veb)
        ALLOCATE(add_edge_midpoint(3,nelem), add_vertex(nnode), vor_veb(2,nelem_voronoi))

        ! first pass: count nnode_voronoi
        add_edge_midpoint = 0
        add_vertex = 0
        nnode_voronoi = nelem
        vor_veb(2,:) = dualb(2,:) - dualb(1,:) + 1 ! first we use the second row of vor_veb to count sizes
        ! loop over all elements
        DO j = 1, nelem
            ! for each edge
            DO jj = 1, 3
                ! if no neighbour
                IF (ne(jj,j) .eq. 0) THEN
                    ! add edge midpoint to the voronoi list
                    IF (add_edge_midpoint(jj,j) .eq. 0) THEN
                        nnode_voronoi = nnode_voronoi + 1
                        add_edge_midpoint(jj,j) = nnode_voronoi
                        ! no need to communicate to the neighbour that the edge midpoint has been 
                        ! added since this happened exactly because there is no neighbour
                    END IF
                    ! add edge nodes to the voronoi list
                    DO kk = 1, 2
                        i = ve(ed(kk,jj),j) ! global node number
                        vor_veb(2,i) = vor_veb(2,i) + 2
                        IF (add_vertex(i) .eq. 0) THEN
                            nnode_voronoi = nnode_voronoi + 1
                            add_vertex(i) = nnode_voronoi
                        END IF
                    END DO
                END IF
            END DO
        END DO

        ! now we compute cumulative indexing and size of vor_ve
        vor_veb(1,1) = 1
        DO i = 2, nelem_voronoi
            vor_veb(1,i) = vor_veb(1,i-1) + vor_veb(2,i-1) ! cumulative index
        END DO
        n = vor_veb(1,nnode) + vor_veb(2,nnode) - 1
        vor_veb(2,:) = vor_veb(1,:) - 1 ! set initially to the first position minus one for each node
        IF (ALLOCATED(vor_ve)) DEALLOCATE(vor_ve)
        IF (ALLOCATED(vor_pt)) DEALLOCATE(vor_pt)
        ALLOCATE(vor_ve(n))
        ALLOCATE(vor_pt(2,nnode_voronoi))
        ! ALLOCATE(movable_node(nnode_voronoi))
        ! movable_node = .true.

        IF (present(mesh_type)) THEN
            mesh_type_selector = mesh_type
        ELSE
            mesh_type_selector = 0
        END IF
        IF (mesh_type_selector .eq. 0) THEN
            ! standard voronoi, but use centroid if the circumcenter lies outside the triangle
            DO j = 1, nelem
                p1 = pt(:,ve(1,j))
                p2 = pt(:,ve(2,j))
                p3 = pt(:,ve(3,j))
                vor_pt(:,j) = circumcenter(p1, p2, p3)
                CALL intriangle(vor_pt(:,j), p1, p2, p3, inside)
                IF (.not. inside) THEN
                    vor_pt(:,j) = sum(pt(:,ve(:,j)), dim=2)*oot
                END IF
            END DO
        ELSE IF (mesh_type_selector .eq. 1) THEN
            ! not an actual voronoi, always use centroids and not circumcenters as vertices for the mesh
            DO j = 1, nelem
                vor_pt(:,j) = sum(pt(:,ve(:,j)), dim=2)*oot
            END DO
        ELSE IF (mesh_type_selector .eq. 2) THEN
            ! not an actual voronoi, always use incenters and not circumcenters as vertices for the mesh
            DO j = 1, nelem
                vor_pt(:,j) = incenter(pt(:,ve(1,j)), pt(:,ve(2,j)), pt(:,ve(3,j)))
            END DO
        ELSE IF (mesh_type_selector .eq. 3) THEN
            ! not an actual voronoi mode for experiments
            DO j = 1, nelem
                p1 = pt(:,ve(1,j))
                p2 = pt(:,ve(2,j))
                p3 = pt(:,ve(3,j))
                IF (obtuse_triangle(p1, p2, p3)) THEN
                    vor_pt(:,j) = incenter(p1, p2, p3)
                ELSE
                    vor_pt(:,j) = sum(pt(:,ve(:,j)), dim=2)*oot
                END IF
            END DO
        END IF
        DO i = 1, nnode
            DO jj = dualb(1,i), dualb(2,i)
                vor_veb(2,i) = vor_veb(2,i) + 1 ! increment the second bound
                vor_ve(vor_veb(2,i)) = dual(jj)         
            END DO
        END DO

        ! second pass: compute the voronoi node coordinates for the boundaries
        DO j = 1, nelem
            DO jj = 1, 3
                IF (ne(jj,j) .eq. 0) THEN
                    ! add edge midpoint to the voronoi list of 
                    ! instead of the mid point of the boundary edge take 
                    ! the barycenter of the triangle projected along the boundary edge or
                    ! the incenter of the triangle projected along the boundary edge
                    vor_pt(:,add_edge_midpoint(jj,j)) = projection_on_segment(&
                        vor_pt(:,j), pt(:,ve(ed(1,jj),j)), pt(:,ve(ed(2,jj),j)))
                    ! movable_node(add_edge_midpoint(jj,j)) = .false.
                    ! add edge nodes to the voronoi list
                    DO kk = 1, 2
                        i = ve(ed(kk,jj),j) ! global node number
                        ! this node will be on the boundary of its voronoi control volume, as will
                        ! the edge midpoint added above
                        add_this_vertex = .true.
                        n = dualb(2,i)-dualb(1,i)+1
                        DO ii = vor_veb(1,i)+n-1+1, vor_veb(2,i)
                            IF (vor_ve(ii) .eq. add_vertex(i)) THEN
                                add_this_vertex = .false.
                                EXIT
                            END IF
                        END DO
                        IF (add_this_vertex) THEN
                            vor_pt(:,add_vertex(i)) = pt(:,i)
                            ! movable_node(add_vertex(i)) = .false.
                            vor_veb(2,i) = vor_veb(2,i) + 1
                            vor_ve(vor_veb(2,i)) = add_vertex(i)
                        END IF
                        vor_veb(2,i) = vor_veb(2,i) + 1
                        vor_ve(vor_veb(2,i)) = add_edge_midpoint(jj,j)
                    END DO
                END IF
            END DO
        END DO

        ! insertion of boundary nodes in the correct order
        DO i = 1, nelem_voronoi
            ! if there are more voronoi vertices than dual elements for node i (or voronoi element i)
            n = dualb(2,i)-dualb(1,i)
            jl = vor_veb(1,i)
            jr = vor_veb(2,i)
            jn = jl + n
            IF (jr - jl .gt. n) THEN
                ! copy from sorted dual
                vor_ve(jl:jn) = dual(dualb(1,i):dualb(2,i))
                ! now zero out the last three spots (in order to detect corners made from a single element)
                vor_ve(jn+1:jn+3) = 0
                ! add the first edge midpoint to the boundary of the control volume, 
                ! immediately following the last triangle
                j = vor_ve(jn)
                jj = ed(1,locate_vertex(ve(:,j), i))
                vor_ve(jn+1) = add_edge_midpoint(jj,j)
                ! add the voronoi control volume vertex to its boundary, in the next spot
                vor_ve(jn+2) = add_vertex(i)
                ! add the second edge midpoint to the boundary of the control volume in the last spot
                j = vor_ve(jl)
                jj = ed(2,locate_vertex(ve(:,j), i))
                vor_ve(jn+3) = add_edge_midpoint(jj,j)
            END IF
        END DO
        DEALLOCATE(add_edge_midpoint, add_vertex)
    END SUBROUTINE build_voronoi_control_volumes

    ! ---- forward and back maps: reference space and physical space

    PURE SUBROUTINE physical_space_to_reference_space_build_map(bbox_center, bbox_scalelen, pt)
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: pt 
        REAL(8), DIMENSION(2),   INTENT(OUT) :: bbox_center
        REAL(8),                 INTENT(OUT) :: bbox_scalelen
        REAL(8)                              :: xmin, xmax, ymin, ymax
        xmin = minval(pt(1,:))
        xmax = maxval(pt(1,:))
        ymin = minval(pt(2,:))
        ymax = maxval(pt(2,:))
        bbox_scalelen = max(xmax - xmin, ymax - ymin)
        bbox_center = [(xmin + xmax), (ymin + ymax)]/2.0d0
    END SUBROUTINE physical_space_to_reference_space_build_map

    PURE SUBROUTINE physical_space_to_reference_space(pt, bbox_center, bbox_scalelen)
        REAL(8), DIMENSION(:,:), INTENT(INOUT)  :: pt 
        REAL(8), DIMENSION(2),   INTENT(IN)     :: bbox_center
        REAL(8),                 INTENT(IN)     :: bbox_scalelen
        INTEGER                                 :: i
        DO i = 1, size(pt,2)
            pt(:,i) = (pt(:,i) - bbox_center)/bbox_scalelen
        END DO
    END SUBROUTINE physical_space_to_reference_space

    PURE SUBROUTINE reference_space_to_physical_space(pt, bbox_center, bbox_scalelen)
        REAL(8), DIMENSION(:,:), INTENT(INOUT)  :: pt 
        REAL(8), DIMENSION(2),   INTENT(IN)     :: bbox_center
        REAL(8),                 INTENT(IN)     :: bbox_scalelen
        INTEGER                                 :: i
        DO i = 1, size(pt,2)
            pt(:,i) = pt(:,i)*bbox_scalelen + bbox_center
        END DO
    END SUBROUTINE reference_space_to_physical_space

END MODULE smesh
