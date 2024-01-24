MODULE input_output

    IMPLICIT NONE

    INTERFACE load
        MODULE PROCEDURE readmatrix
        MODULE PROCEDURE readmatrix_integer
        MODULE PROCEDURE readvector
        MODULE PROCEDURE readvector_integer
    END INTERFACE load

    ! BINARY IO
    INTERFACE prt_bin
        MODULE PROCEDURE store_binary_matrix_real
        MODULE PROCEDURE store_binary_matrix_integer
        MODULE PROCEDURE store_binary_vector_real
        MODULE PROCEDURE store_binary_vector_integer
        MODULE PROCEDURE store_binary_matrix3_real
        MODULE PROCEDURE store_binary_matrix3_integer
    END INTERFACE prt_bin

    INTERFACE load_bin
        MODULE PROCEDURE load_binary_matrix_real
        MODULE PROCEDURE load_binary_matrix3_real
        MODULE PROCEDURE load_binary_matrix3_integer
        MODULE PROCEDURE load_binary_matrix_integer
        MODULE PROCEDURE load_binary_vector_real
        MODULE PROCEDURE load_binary_vector_integer
    END INTERFACE load_bin

   
    PRIVATE
    ! READER
    PUBLIC :: load
    ! BINARY IO
    PUBLIC :: prt_bin
    PUBLIC :: load_bin

CONTAINS

! READER ------------------------------------------------------------------------------------------
! Legge da file di testo:
!
!     Una matrice o un vettore:  CALL load("matrice_o_vettore.txt", matrice_o_vettore_allocatable)
!                                     - La procedura \`e interfaccia alle due procedure di livello
!                                       pi\`u basso per matrici o vettori di seguito descritte
!
!     Un vettore: CALL readvector("vettore.txt", vettore_allocatable)
!                      - Gli elementi sono separati per riga
!                      - Tutte le celle devono essere riempite, NaN (non case-sensitive) per i valori non presenti
!
!     Una matrice: CALL readmatrix("matrice.txt", matrice_allocatable)
!                       - Le celle si separano con un qualsiasi numero di spazi o tabulazioni
!                       - Tutte le celle devono essere riempite, NaN (non case-sensitive) per i valori non presenti
!     

    SUBROUTINE readmatrix(filename, matrix, transpose, verbose)
        REAL(8), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: matrix
        CHARACTER(len=*), INTENT(IN)                         :: filename
        LOGICAL, OPTIONAL, INTENT(IN)                     :: transpose
        LOGICAL, OPTIONAL, INTENT(IN)                     :: verbose
        LOGICAL                                           :: do_transpose
        INTEGER                                           :: buffer, fu
        CHARACTER, DIMENSION(:), ALLOCATABLE              :: buffer_test_line
        CHARACTER(len=:), ALLOCATABLE                     :: line
        CHARACTER(len=:), ALLOCATABLE                     :: microline
        CHARACTER(len=100)                                   :: fmt
        INTEGER                                           :: rows, cols, cursor, status_code, read_chars
        INTEGER                                           :: i, j
        LOGICAL                                           :: there, be_verbose
        IF (present(verbose)) THEN
            be_verbose = verbose
        ELSE
            be_verbose = .true.
        END IF
        INQUIRE(file=filename, exist=there) 
        IF (there) THEN
            IF (PRESENT(transpose)) THEN
                do_transpose = transpose
            ELSE
                do_transpose = .false.
            END IF
            ALLOCATE(CHARACTER(1) :: microline)
            fu     = 6550
            rows   = 0
            cols   = 0
            cursor = 1
            ! open files
            OPEN(unit=fu, file=filename, status="old", access='sequential', action='read')
            ! count rows
            DO 
                READ(fu, '(A)', iostat=status_code) microline
                IF (status_code < 0) exit
                rows = rows + 1
            END DO
            IF (rows .ne. 0) THEN
                ! rewind file
                REWIND(unit=fu)
                ! count columns of first row and get buffer size
                buffer = 256
                ALLOCATE(buffer_test_line(buffer))
                read_chars = buffer
                DO ! equals should be enough
                    WRITE(fmt,"(A,I0,A)") "(", buffer, "a)"
                    READ(fu, fmt, advance="no", size=read_chars, iostat=status_code) buffer_test_line
                    ! WRITE(*,"(A,I0,A,I0)") "READ MATRIX: read_chars: ", read_chars, ", status_code: ", status_code
                    IF (buffer <= read_chars) THEN ! equals should be enough
                        buffer = buffer*2
                        ! WRITE(*,"(A,I0)") "READ MATRIX: extending buffer size to ", buffer
                        DEALLOCATE(buffer_test_line)
                        ALLOCATE(buffer_test_line(buffer))
                        REWIND(unit=fu)
                    ELSE
                        IF (be_verbose) THEN
                            WRITE(*,"(A, A, A, I0)") "Reading matrix from file ", trim(filename), ". Buffer size: ", buffer
                        END IF
                        EXIT
                    END IF
                END DO
                ! no not shrink the buffer size to fit the first line, the next lines could still be longer
                ! remember to throw an error if a line is too long, if it does not prove to be too much of an effort
                REWIND(unit=fu) ! might be not necessary due to advance no
                ALLOCATE(CHARACTER(buffer) :: line)
                READ(fu, '(A)') line
                DO WHILE(cursor /= 0)
                    cursor = 1
                    cols   = cols+1
                    DO WHILE(cursor == 1)
                        cursor = SCAN(trim(line), " "//CHAR(9), .FALSE.)
                        line = line(cursor+1:)
                    END DO
                END DO
                IF (rows == 0) THEN
                    IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
                    ALLOCATE(matrix(1,1))
                    matrix = 0.0d0
                    WRITE(*,*) "WARNING, EMPTY FILE: ", filename
                ELSE
                    ! rewind file
                    REWIND(unit=fu) ! maybe not needed
                    ! allocate memory
                    IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
                    IF (do_transpose) THEN
                        ALLOCATE(matrix(cols, rows))
                        ! read file
                        DO i=1,rows
                            READ(fu,*) (matrix(j,i), j=1,cols)
                        END DO
                    ELSE
                        ALLOCATE(matrix(rows, cols))
                        ! read file
                        DO i=1,rows
                            READ(fu,*) (matrix(i,j), j=1,cols)
                        END DO
                    END IF
                END IF
                CLOSE(unit=fu)
            ELSE
                IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
                ALLOCATE(matrix(0,0))
            END IF
        ELSE
            IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
            ALLOCATE(matrix(1,1))
            matrix = 0.0d0
            WRITE(*,*) "WARNING, FILE NOT FOUND: ", filename
        END IF
    END SUBROUTINE readmatrix

    SUBROUTINE readmatrix_integer(filename, matrix, transpose, verbose)
        ! should be just the same as above except for it reads integers, and the zero
        ! MAYBE more care should be given in reading a couple more lines in order to determine buffer size
        ! or read again longer lines?
        INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: matrix
        CHARACTER(len=*), INTENT(IN)                         :: filename
        LOGICAL, OPTIONAL, INTENT(IN)                     :: transpose
        LOGICAL, OPTIONAL, INTENT(IN)                     :: verbose
        LOGICAL                                           :: do_transpose
        INTEGER                                           :: buffer, fu
        CHARACTER, DIMENSION(:), ALLOCATABLE              :: buffer_test_line
        CHARACTER(len=:), ALLOCATABLE                     :: line
        CHARACTER(len=:), ALLOCATABLE                     :: microline
        CHARACTER(len=100)                                   :: fmt
        INTEGER                                           :: rows, cols, cursor, status_code, read_chars
        INTEGER                                           :: i, j
        LOGICAL                                           :: there, be_verbose
        IF (present(verbose)) THEN
            be_verbose = verbose
        ELSE
            be_verbose = .true.
        END IF
        INQUIRE(file=filename, exist=there) 
        IF (there) THEN
            IF (PRESENT(transpose)) THEN
                do_transpose = transpose
            ELSE
                do_transpose = .false.
            END IF
            ALLOCATE(CHARACTER(1) :: microline)
            fu     = 6550
            rows   = 0
            cols   = 0
            cursor = 1
            ! open files
            OPEN(unit=fu, file=filename, status="old", access='sequential', action='read')
            ! count rows
            DO 
                READ(fu, '(A)', iostat=status_code) microline
                IF (status_code < 0) exit
                rows = rows + 1
            END DO
            ! rewind file
            REWIND(unit=fu)
            ! count columns of first row and get buffer size
            buffer = 256
            ALLOCATE(buffer_test_line(buffer))
            read_chars = buffer
            DO ! equals should be enough
                WRITE(fmt,"(A,I0,A)") "(", buffer, "a)"
                READ(fu, fmt, advance="no", size=read_chars, iostat=status_code) buffer_test_line
                ! WRITE(*,"(A,I0,A,I0)") "READ MATRIX: read_chars: ", read_chars, ", status_code: ", status_code
                IF (buffer <= read_chars) THEN ! equals should be enough
                    buffer = buffer*2
                    ! WRITE(*,"(A,I0)") "READ MATRIX: extending buffer size to ", buffer
                    DEALLOCATE(buffer_test_line)
                    ALLOCATE(buffer_test_line(buffer))
                    REWIND(unit=fu)
                ELSE
                    IF (be_verbose) THEN
                        WRITE(*,"(A,I0)") "READ MATRIX: adequate buffer size ", buffer
                    END IF
                    EXIT
                END IF
            END DO
            ! no not shrink the buffer size to fit the first line, the next lines could still be longer
            ! remember to throw an error if a line is too long, if it does not prove to be too much of an effort
            REWIND(unit=fu) ! might be not necessary due to advance no
            ALLOCATE(CHARACTER(buffer) :: line)
            READ(fu, '(A)') line
            DO WHILE(cursor /= 0)
                cursor = 1
                cols   = cols+1
                DO WHILE(cursor == 1)
                    cursor = SCAN(trim(line), " "//CHAR(9), .FALSE.)
                    line = line(cursor+1:)
                END DO
            END DO
            IF (rows == 0) THEN
                IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
                ALLOCATE(matrix(1,1))
                matrix = 0
                WRITE(*,*) "WARNING, EMPTY FILE: ", filename
            ELSE
                ! rewind file
                REWIND(unit=fu) ! maybe not needed
                ! allocate memory
                IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
                IF (do_transpose) THEN
                    ALLOCATE(matrix(cols, rows))
                    ! read file
                    DO i=1,rows
                        READ(fu,*) (matrix(j,i), j=1,cols)
                    END DO
                ELSE
                    ALLOCATE(matrix(rows, cols))
                    ! read file
                    DO i=1,rows
                        READ(fu,*) (matrix(i,j), j=1,cols)
                    END DO
                END IF
            END IF
            CLOSE(unit=fu)
        ELSE
            IF (ALLOCATED(matrix)) DEALLOCATE(matrix)
            ALLOCATE(matrix(1,1))
            matrix = 0
            WRITE(*,*) "WARNING, FILE NOT FOUND: ", filename
        END IF
    END SUBROUTINE readmatrix_integer

    SUBROUTINE readvector(filename, vector)
        REAL(8), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: vector
        CHARACTER(len=*), INTENT(IN)                       :: filename
        INTEGER                                         :: buffer, fu
        CHARACTER(len=:), ALLOCATABLE                   :: line
        INTEGER                                         :: rows, status_code
        INTEGER                                         :: i
        LOGICAL                                         :: there 
        INQUIRE(file=filename, exist=there) 
        IF (there) THEN
            buffer = 100 ! how long can a vector be??
            ALLOCATE(CHARACTER(buffer) :: line)
            fu   = 6550
            rows = 0
            ! open files
            OPEN(unit=fu, file=filename, status="old", access='sequential', action='read')
            ! count rows
            DO 
                READ(fu, '(A)', iostat=status_code) line
                IF (status_code < 0) exit
                rows = rows + 1
            END DO
            IF (rows == 0) THEN
                IF (ALLOCATED(vector)) DEALLOCATE(vector)
                ALLOCATE(vector(1))
                vector = 0.0d0
                WRITE(*,*) "WARNING, FILE IS EMPTY: ", filename
            ELSE
                ! rewind file
                REWIND(unit=fu)
                ! allocate memory
                IF (ALLOCATED(vector)) DEALLOCATE(vector)
                ALLOCATE(vector(rows))
                ! read file
                DO i=1,rows
                    READ(fu,*) vector(i)
                END DO
            END IF
            CLOSE(unit=fu)
        ELSE
            IF (ALLOCATED(vector)) DEALLOCATE(vector)
            ALLOCATE(vector(1))
            vector = 0.0d0
            WRITE(*,*) "WARNING, FILE NOT FOUND: ", filename
        END IF
    END SUBROUTINE readvector

    SUBROUTINE readvector_integer(filename, vector)
        INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: vector
        CHARACTER(len=*), INTENT(IN)                       :: filename
        INTEGER                                         :: buffer, fu
        CHARACTER(len=:), ALLOCATABLE                   :: line
        INTEGER                                         :: rows, status_code
        INTEGER                                         :: i
        LOGICAL                                         :: there 
        INQUIRE(file=filename, exist=there) 
        IF (there) THEN
            buffer = 100 ! how long can a vector be??
            ALLOCATE(CHARACTER(buffer) :: line)
            fu   = 6550
            rows = 0
            ! open files
            OPEN(unit=fu, file=filename, status="old", access='sequential', action='read')
            ! count rows
            DO 
                READ(fu, '(A)', iostat=status_code) line
                IF (status_code < 0) exit
                rows = rows + 1
            END DO
            IF (rows == 0) THEN
                IF (ALLOCATED(vector)) DEALLOCATE(vector)
                ALLOCATE(vector(1))
                vector = 0
                WRITE(*,*) "WARNING, FILE IS EMPTY: ", filename
            ELSE
                ! rewind file
                REWIND(unit=fu)
                ! allocate memory
                IF (ALLOCATED(vector)) DEALLOCATE(vector)
                ALLOCATE(vector(rows))
                ! read file
                DO i=1,rows
                    READ(fu,*) vector(i)
                END DO
            END IF
            CLOSE(unit=fu)
        ELSE
            IF (ALLOCATED(vector)) DEALLOCATE(vector)
            ALLOCATE(vector(1))
            vector = 0
            WRITE(*,*) "WARNING, FILE NOT FOUND: ", filename
        END IF
    END SUBROUTINE readvector_integer

! BINARY INPUT/OUTPUT

   SUBROUTINE store_binary_vector_real(V, filename)
        REAL(8), DIMENSION(:), INTENT(IN) :: V
        CHARACTER(len=*), INTENT(IN)         :: filename
        INTEGER, DIMENSION(1)             :: V_shape
        V_shape = shape(V)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) V_shape
        WRITE(6550) V
        CLOSE(6550)
    END SUBROUTINE store_binary_vector_real

    SUBROUTINE store_binary_matrix_real(M, filename)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: M
        CHARACTER(len=*), INTENT(IN)           :: filename
        INTEGER, DIMENSION(2)               :: M_shape
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_matrix_real

    SUBROUTINE load_binary_matrix_real(filename, M)
        CHARACTER(len=*), INTENT(IN)                         :: filename
        REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(2)                             :: M_shape
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_matrix_real

    SUBROUTINE load_binary_vector_real(filename, V)
        CHARACTER(len=*), INTENT(IN)                       :: filename
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: V
        INTEGER, DIMENSION(1)                           :: V_shape
        IF (ALLOCATED(V)) DEALLOCATE(V)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) V_shape
        ALLOCATE(V(V_shape(1)))
        READ (6551) V
        CLOSE(6551)
    END SUBROUTINE load_binary_vector_real

! uguale ma per interi

    SUBROUTINE store_binary_vector_integer(V, filename)
        INTEGER, DIMENSION(:), INTENT(IN) :: V
        CHARACTER(len=*), INTENT(IN)         :: filename
        INTEGER, DIMENSION(1)             :: V_shape
        V_shape = shape(V)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) V_shape
        WRITE(6550) V
        CLOSE(6550)
    END SUBROUTINE store_binary_vector_integer

    SUBROUTINE store_binary_matrix_integer(M, filename)
        INTEGER, DIMENSION(:,:), INTENT(IN) :: M
        CHARACTER(len=*), INTENT(IN)           :: filename
        INTEGER, DIMENSION(2)               :: M_shape
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_matrix_integer

    SUBROUTINE load_binary_matrix_integer(filename, M)
        CHARACTER(len=*), INTENT(IN)                         :: filename
        INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(2)                             :: M_shape
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_matrix_integer

    SUBROUTINE load_binary_vector_integer(filename, V)
        CHARACTER(len=*), INTENT(IN)                       :: filename
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: V
        INTEGER, DIMENSION(1)                           :: V_shape
        IF (ALLOCATED(V)) DEALLOCATE(V)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) V_shape
        ALLOCATE(V(V_shape(1)))
        READ (6551) V
        CLOSE(6551)
    END SUBROUTINE load_binary_vector_integer

    SUBROUTINE store_binary_matrix3_real(M, filename)
        REAL(8), DIMENSION(:,:,:), INTENT(IN) :: M
        CHARACTER(len=*), INTENT(IN)           :: filename
        INTEGER, DIMENSION(3)               :: M_shape
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_matrix3_real

    SUBROUTINE load_binary_matrix3_real(filename, M)
        CHARACTER(len=*), INTENT(IN)                         :: filename
        REAL(8), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(3)                             :: M_shape
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2), M_shape(3)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_matrix3_real

    SUBROUTINE store_binary_matrix3_integer(M, filename)
        INTEGER, DIMENSION(:,:,:), INTENT(IN) :: M
        CHARACTER(len=*), INTENT(IN)           :: filename
        INTEGER, DIMENSION(3)               :: M_shape
        M_shape = shape(M)
        OPEN(unit=6550, file=filename, form="unformatted", access="sequential", action="write", status="replace")
        WRITE(6550) M_shape
        WRITE(6550) M
        CLOSE(6550)
    END SUBROUTINE store_binary_matrix3_integer

    SUBROUTINE load_binary_matrix3_integer(filename, M)
        CHARACTER(len=*), INTENT(IN)                         :: filename
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: M
        INTEGER, DIMENSION(3)                             :: M_shape
        IF (ALLOCATED(M)) DEALLOCATE(M)
        OPEN(unit=6551, file=filename, form="unformatted", access="sequential")
        READ (6551) M_shape
        ALLOCATE(M(M_shape(1), M_shape(2), M_shape(3)))
        READ (6551) M
        CLOSE(6551)
    END SUBROUTINE load_binary_matrix3_integer


END MODULE input_output

! for reading binary files in python
    ! def load_fortran_binary(filename):
    !     # Fortran data structure:
    !     #     4 byte head --> N
    !     #     N byte first record (int32, the shape of the object contained in the second record)
    !     #     4 byte foot == head
    !     #     4 byte head --> M
    !     #     M byte second record (int32 or float64, the matrix or vector)
    !     #     4 byte foot == head
    !     with open(filename) as f:
    !         # read shape
    !         record_length = np.fromfile(f, dtype=np.int32, count=1, sep="")[0]
    !         shape = np.fromfile(f, dtype=np.int32, count=record_length//4, sep="")
    !         dummy = np.fromfile(f, dtype=np.int32, count=1, sep="")
    !         # compute matrix size
    !         nelem = np.prod(shape)
    !         # here we use the record length to determine the type of data (8bit/element -> real, 4bit/element -> integer)
    !         record_length = np.fromfile(f, dtype=np.int32, count=1, sep="")[0]
    !         # select data type
    !         if record_length//nelem == 4:
    !             dtype = np.int32
    !         else:
    !             dtype = np.float64
    !         # read data
    !         if shape.size > 1:
    !             return np.fromfile(f, dtype=dtype, count=nelem, sep="").reshape(shape[::-1]).transpose()
    !         else:
    !             return np.fromfile(f, dtype=dtype, count=nelem, sep="").flatten()
    