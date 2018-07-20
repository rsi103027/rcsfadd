PROGRAM RCSFADD

!*********************************************************************** 
! This program adds the CSFs in the second list to the end of the first list.
! Three options are provided:
! 1) Add the CSFs in the second list to the end of the first list without checking duplicates
! 2) Only the CSFs in the second list have n bigger than nmax of the first list are added
! 3) Only the CSFs in the second list have different core electrons with the first list are added
! Output file: rcsf.out
!*********************************************************************** 

implicit none
integer :: I, IOS
character :: list1*100, list2*100
integer :: ADD, NUM, NBLOCK, BLOCK, NCORE
CHARACTER :: RECORD_1*15, RECORD_2*15
CHARACTER :: S_closed_1*181, S_closed_2*181
CHARACTER :: S_orbitals_1*1070, S_orbitals_2*1070
INTEGER ::  NCFBLK_1, NCFBLK_2
CHARACTER(LEN=256) ::C_CORE(100)

write(*,*) ' Three options are provided: &
1. Add the CSFs in the second list to the end of the first list without checking duplicates. &
2. Only the CSFs in the second list have n bigger than nmax of the first list are added. &
3. Only the CSFs in the second list have different core electrons with the first list are added.'
read(*,*)add
if(add==1)then
        continue
elseif(add==2)then
        write(*,*)'Input the max n quntum number of list1'
        read(*,*)num
elseif(add==3)then
        write(*,*)'Input the number of relativistic core subshells except for closed shells'
        read(*,*)num
else
        stop 'Wrong input'
end if
write(*,*) 'Full name of list1'
read(*,'(a)') list1
write(*,*) 'Full name of list2'
read(*,'(a)') list2
write(*,*) 'Number of blocks'
read(*,'(i2)') nblock

open(unit=21,file=trim(list1),status='old')
open(unit=20,file=trim(list2),status='old')
open(unit=22,file='rcsf.out',status='unknown')

READ (21, '(1A15)', IOSTAT=IOS) RECORD_1
READ (20, '(1A15)', IOSTAT=IOS) RECORD_2
WRITE (22, '(1A15)') RECORD_1

READ (21, '(A)') S_closed_1
READ (20, '(A)') S_closed_2
IF(S_closed_1 /= S_closed_2) then
	STOP "Diffeent close shells"
END IF
I = LEN_TRIM(S_closed_1)
WRITE (22,'(A)') S_closed_1(1:I)

READ (21, '(1A7)', IOSTAT=IOS) RECORD_1
READ (20, '(1A7)', IOSTAT=IOS) RECORD_2
WRITE (22,'(1A7)') RECORD_1

READ (21, '(A)') S_orbitals_1
READ (20, '(A)') S_orbitals_2
I = LEN_TRIM(S_orbitals_1)
IF(LEN_TRIM(S_orbitals_1) >= LEN_TRIM(S_orbitals_2)) THEN
	WRITE (22, '(A)') TRIM(S_orbitals_1)
ELSE
	WRITE (22, '(A)') TRIM(S_orbitals_2)
END IF

READ (21, '(1A7)', IOSTAT=IOS) RECORD_1
READ (20, '(1A7)', IOSTAT=IOS) RECORD_2
WRITE (22,'(1A7)') RECORD_1
WRITE (*, *) "  Block    First List   Complete List"

DO BLOCK=1, NBLOCK
        CALL LODCSL_FIRST(NCFBLK_1,ADD,NUM,NCORE)
        CALL LODCSL_SECOND(NCFBLK_2,ADD,NUM,NCORE)
	write(*,'(3X,I2,3X,I14,I16)') BLOCK, NCFBLK_1, NCFBLK_1+NCFBLK_2
	IF(BLOCK == NBLOCK) EXIT
	WRITE(22,'(A2)') ' *'
END DO

STOP
END PROGRAM RCSFADD

SUBROUTINE LODCSL_FIRST(NCFBLK_1,ADD,NUM,NCORE)
integer :: IOS, NCF, N, NEW_CORE, ADD, NUM
INTEGER, INTENT(OUT)::  NCFBLK_1, NCORE
CHARACTER(LEN=256):: RECORD_C_shell, RECORD_C_quant, RECORD_C_coupl
CHARACTER(LEN=2):: N_C
CHARACTER(LEN=256):: C_CORE(100)

NCFBLK_1 = 0
NCF = 0
NCORE = 1
DO
	READ (21, '(A)', IOSTAT=IOS) RECORD_C_shell
	IF (IOS == 0) THEN
		IF (RECORD_C_shell(1:2) == ' *') THEN
			NCFBLK_1 = NCF
			RETURN
		ENDIF
		NCF = NCF + 1
                IF(ADD==2)THEN
		        N_C=RECORD_C_shell(len_trim(RECORD_C_shell)-7:len_trim(RECORD_C_shell)-6)
		        read(N_C,'(i2)')N
                        IF(N>NUM)THEN
                                STOP 'Maxn in list1 is larger than input'
                        ENDIF
                ENDIF
                IF(ADD==3)THEN
		        NEW_CORE=1
        		do I=1,NCORE
	        		if(RECORD_C_shell(1:NUM*9)==C_CORE(I)(1:NUM*9))then
		        		NEW_CORE=NEW_CORE-1
			        endif
        		enddo
	        	if(NEW_CORE==1)then
		        	C_CORE(NCORE) = RECORD_C_shell(1:NUM*9)
			        NCORE=NCORE+1
		        endif
                ENDIF
		WRITE(22,'(A)') TRIM(RECORD_C_shell)
		READ (21, '(A)',  IOSTAT=IOS) RECORD_C_quant
		WRITE(22,'(A)') TRIM(RECORD_C_quant)
		READ (21, '(A)',  IOSTAT=IOS) RECORD_C_coupl
		WRITE(22,'(A)') TRIM(RECORD_C_coupl)
	ELSE
		EXIT
	ENDIF
END DO
NCFBLK_1 = NCF
RETURN
END SUBROUTINE LODCSL_FIRST

SUBROUTINE LODCSL_SECOND(NCFBLK_2,ADD,NUM,NCORE)
integer :: IOS, NCF, I, N, FOUND, ADD, NCORE
INTEGER, INTENT(OUT) ::  NCFBLK_2
CHARACTER(LEN=256) :: RECORD_C_shell, RECORD_C_quant, RECORD_C_coupl
CHARACTER(LEN=2)   :: N_C
CHARACTER(LEN=256):: C_CORE(100)

NCF = 0
FOUND = 0
DO
	READ (20, '(A)', IOSTAT=IOS) RECORD_C_shell
	IF (IOS == 0) THEN
		IF (RECORD_C_shell(1:2) == ' *') THEN
			NCFBLK_2 = NCF
			RETURN
		ENDIF
		READ (20, '(A)',  IOSTAT=IOS) RECORD_C_quant
		READ (20, '(A)',  IOSTAT=IOS) RECORD_C_coupl
                IF(ADD==1)THEN
		        NCF = NCF + 1  
		        WRITE(22,'(A)') TRIM(RECORD_C_shell)
		        WRITE(22,'(A)') TRIM(RECORD_C_quant)
		        WRITE(22,'(A)') TRIM(RECORD_C_coupl)
                ELSEIF(ADD==2)THEN
		        N_C=RECORD_C_shell(len_trim(RECORD_C_shell)-7:len_trim(RECORD_C_shell)-6)
		        read(N_C,'(i2)')N
		        IF(N > NUM) THEN
			        NCF = NCF + 1  
			        WRITE(22,'(A)') TRIM(RECORD_C_shell)
			        WRITE(22,'(A)') TRIM(RECORD_C_quant)
			        WRITE(22,'(A)') TRIM(RECORD_C_coupl)
                        ENDIF
                ELSEIF(ADD==3)THEN
                        FOUND=0
		        DO I=1,NCORE
			        IF(RECORD_C_shell(1:NUM*9)==C_CORE(I)(1:NUM*9))THEN
				        FOUND=1
			        ENDIF
		        END DO
		        IF(FOUND==0)THEN
			        NCF = NCF + 1  
			        WRITE(22,'(A)') TRIM(RECORD_C_shell)
			        WRITE(22,'(A)') TRIM(RECORD_C_quant)
			        WRITE(22,'(A)') TRIM(RECORD_C_coupl)
		        ENDIF
		ENDIF
	ELSE
		EXIT
                
	ENDIF
END DO
NCFBLK_2 = NCF
RETURN
END SUBROUTINE LODCSL_SECOND

