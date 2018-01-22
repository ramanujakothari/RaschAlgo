	module wig_sym
	
	use fgsl
	use parameters
	
	implicit none
	
	double precision, allocatable, dimension(:,:):: wig_arr
	
	contains
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				initialize the wigner 3j computation
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine init_wig()
	
	max_index = max_len(l_min, l_max)
	allocate(wig_arr(1:max_index, 2))
	call wignstore(l_min, l_max, max_index, wig_arr)
	
	end subroutine init_wig
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Maximum length of the Wigner storing array
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function max_len(l_min, l_max)
	
	implicit none
	
	integer, intent(in) :: l_max, l_min
	integer::l1, l2, l3, m1, m2, m3, counter, counter1
	integer:: max_len, index_loc, no_trans
	real, dimension(2,3) :: wigin
	integer*8 :: tot_num

	max_len = 0
	tot_num = 0
	
	do l1=l_min,l_max+ecare
		do l2=l_min,l_max+ecare
			do l3=abs(l1-l2), l1+l2 !ensure triangular condition
				do m1 = -l1, l1
					do m2 = -l2, l2
						m3 = -m1-m2
						if(abs(m3) <= l3) then
							call par2wig(1.0*l1, 1.0*l2, 1.0*l3, 1.0*m1, 1.0*m2, 1.0*m3, wigin)
							call wignstandindex(wigin, index_loc, no_trans)
							tot_num = tot_num + 1
							if(index_loc > max_len) then
								max_len = index_loc
							end if
						end if
					end do
				end do
			end do
		end do
	end do
	
	if(want_nsymb .eqv. .true.) then
		print*, 'Symbols Required by Code:', tot_num 
	end if
	
	!print*, 'Maximum value of Index', max_len
	
	end function
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Indexing function
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	this indexing scheme has been taken from Rasch & Yu SIAM 25, 1416
	function index_func(wigin)
	
	implicit none
	
	real, intent(in), dimension(2,3):: wigin
	integer:: L, X, T, B, S
	integer:: index_func
	
	S = -wigin(1,1) + wigin(1,2) + wigin(1,3)
	L = wigin(1,1) - wigin(1,2) + wigin(1,3)
	X = wigin(1,1) - wigin(2,1)
	B = wigin(1,2) - wigin(2,2)
	T = wigin(1,3) + wigin(2,3)
	
	index_func = reg2index(L, X, T, B, S)
	
	end function index_func

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Converts Regge enteries into index
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function reg2index(L, X, T, B, S)
	
	implicit none
	
	integer :: L, X, T, B, S, reg2index
	integer*8 :: term1, term2, term3, term4

!	the strange _8 has been concatenated to convert into 8 bit integer
	term1 = (L**5_8 + 1_8*10*L**4 + 1_8*35*L**3 + 1_8*50*L**2 + 1_8*24*L)/120
	term2 = (6_8*X + 11_8*X**2 + 6_8*X**3 + 1_8*X**4)/24
	term3 = (2_8*T + 3_8*T**2 + 1_8*T**3)/6
	term4 = (1_8*B**2 + 1_8*B)/2
	reg2index = term1 + term2 + term3 + term4 + S + 1
	
	if (reg2index<0 .or. term1<0 .or. term2<0 .or. term3<0 .or. term4<0) then
		print*, 'Integer overflow in reg2index'
		stop
	endif
	
	end function reg2index
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Wigner to Regge Symbol Conversion
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine wig2reg(wigin, regmat)
	
	implicit none
	
	real, intent(in), dimension(2,3):: wigin
	integer, intent(out), dimension(3,3) :: regmat
	
	regmat(1,1) = -wigin(1,1)+wigin(1,2)+wigin(1,3)
	regmat(1,2) = wigin(1,1)-wigin(1,2)+wigin(1,3)
	regmat(1,3) = wigin(1,1)+wigin(1,2)-wigin(1,3)
	regmat(2,1) = wigin(1,1)-wigin(2,1)
	regmat(2,2) = wigin(1,2)-wigin(2,2)
	regmat(2,3) = wigin(1,3)-wigin(2,3)
	regmat(3,1) = wigin(1,1)+wigin(2,1)
	regmat(3,2) = wigin(1,2)+wigin(2,2)
	regmat(3,3) = wigin(1,3)+wigin(2,3)
	
	end subroutine wig2reg

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Subroutine for swapping rows or columns
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine swap(a,b)
	
	implicit none
	
	integer :: i
	integer, dimension(3) :: a, b, temp
	
	do i=1,3
		temp(i) = a(i)
		a(i) = b(i)
		b(i) = temp(i)
	end do
	
	end subroutine
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Stopping Condition for Interchange of rows &/or columns
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function stopcond(l, x, t, b, s)
	
	implicit none
	
	integer :: l, x, t, b, s
	logical :: stopcond
	
	stopcond = .true.
	if ((l >= x) .and. (x >= t) .and. (t >= b) .and. (b >= s)) then
		stopcond = .false.
	end if
	
	end function
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Print the matrix
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine printmat(a,m,n)
!	for most of the time this will be used to print the integer matrix
!	but sometimes it might be 1/2 integer matrix as well
	
	implicit none
	
	real, dimension(m,n) :: a
	integer:: i,j,m,n
	
	print*, ''
	
	do i=1,m
		print*, (a(i,j),j=1,n)
	end do
	
	print*, ''
	
	end subroutine printmat

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!					Location of maximum element of matrix
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine locmax(mat, rowloc, colloc)
	
	implicit none
	
	integer, intent(in), dimension(3,3) :: mat
	integer, intent(out) :: colloc, rowloc
	integer :: i, j, maxi
	
	maxi = 0
	
	do i=1,3
		do j=1,3
			if(maxi <= mat(i,j)) then
				maxi = mat(i,j)
				rowloc = i
				colloc = j
			end if
		end do
	end do
	
	end subroutine
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	                Wigner 3j Symbol Calculation
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function wigner3j(wigin)
	
	implicit none
	
	real, intent(in), dimension(2,3):: wigin
	integer(fgsl_int), dimension(2,3):: wig
!	here I'm doing the type conversion from the real to GSL's fgsl_int
	double precision :: wigner3j
	
	wig = 2.0*wigin
	wigner3j = fgsl_sf_coupling_3j(wig(1,1), wig(1,2), wig(1,3), &
		& wig(2,1), wig(2,2), wig(2,3))
	
	end function
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!						Parameter to Wigner Matrix
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine par2wig(l1, l2, l3, m1, m2, m3, wigmat)
	
	implicit none
	
	real, intent(in):: l1, l2, l3, m1, m2, m3
	real, intent(out), dimension(2,3) :: wigmat

	wigmat(1,1) = 1.0*l1
	wigmat(1,2) = 1.0*l2
	wigmat(1,3) = 1.0*l3
	wigmat(2,1) = 1.0*m1
	wigmat(2,2) = 1.0*m2
	wigmat(2,3) = 1.0*m3
	
	end subroutine

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Calculate Index after standardisation
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	1. find the largest element in matrix
!	2. swap the column containing that element with column 2
!	3. now swap the row containing maximum element with row 1
!	4. this makes a12 as the maximum element of the array
!	5. Check whether L>=X>=T>=B>=S is satisfied
!	6. If not swap col1 and col3
!	7. Check which row has maximum element
!	8. Swap it with row 2
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine wignstandindex(wigin, index_loc, no_trans)
	
	implicit none
	
	real, intent(in), dimension(2,3) :: wigin
	integer, intent(out) :: index_loc, no_trans
	integer, dimension(3,3):: rmat
	integer L, X, T, B, S, temp, locmaxrow, locmaxcol
	integer, dimension(3) :: row_maxcol, col_maxcol
	
	no_trans = 0
	
	call wig2reg(wigin, rmat)
!	perform row, column transformations or transpose iff required
	if(stopcond(rmat(1,2),rmat(2,1),rmat(3,3),rmat(2,2),rmat(1,1))) then
! 	determine which column has the maximum element present,
!	locmaxcol is location of column which contains the maximum element of matrix
		call locmax(rmat, locmaxrow, locmaxcol)
!	determine whether minimum and maximum elements are in a row or column
		col_maxcol = rmat(:,locmaxcol)
		row_maxcol = rmat(locmaxrow,:)

!	transpose operation
		if(minval(row_maxcol) > minval(col_maxcol))then
!	this corresponds to the case when maximum and minimum values appear in a column
			rmat = transpose(rmat)
!	reassign wigin to the new one for the next transformation
			temp = locmaxrow
			locmaxrow = locmaxcol
			locmaxcol = temp
		end if
!		row swapping
		if(locmaxrow /= 1) then
			call swap(rmat(1,:), rmat(locmaxrow,:))
			no_trans = no_trans + 1
		endif
!		column swapping
		if(locmaxcol /= 2) then
			call swap(rmat(:,2), rmat(:,locmaxcol))
			no_trans = no_trans + 1
		endif
		
		if(rmat(1,1) > rmat(1,3)) then
			call swap(rmat(:,1), rmat(:,3))
			no_trans = no_trans + 1
		end if
		
		if(rmat(2,2) > rmat(3,2)) then
			call swap(rmat(2,:), rmat(3,:))
			no_trans = no_trans + 1
		end if
		
		if(rmat(2,2) == rmat(3,2) .and. rmat(2,3) > rmat(3,3)) then
			call swap(rmat(2,:), rmat(3,:))
			no_trans = no_trans + 1
		end if
	end if

	if(stopcond(rmat(1,2),rmat(2,1),rmat(3,3),rmat(2,2),rmat(1,1)) .eqv. .True.) then
		print*, 'Problem in Algo Implementation'
	end if
!	calculate the corresponding index using L, X, T, B, S
	L = rmat(1,2)
	S = rmat(1,1)
	T = rmat(3,3)
	B = rmat(2,2)
	X = rmat(2,1)
	
	index_loc = reg2index(L, X, T, B, S)
	
	end subroutine wignstandindex

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Calculation of phase factor
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!	1. 	if l1+l2+l3 is even then the phase will always be 1
!	2. 	if l1+l2+l3 is odd, 
!	2a.	if the sum of no of transformations of the stored wigner symbol 
!		& the one considered is even then phase=1
!	2b. if the sum is odd then phase is -1
	function phase_fac(wigin, trans1, trans2)
	
	implicit none
	
	real, intent(in), dimension(2,3) :: wigin
	integer, intent(in) :: trans1, trans2
	integer:: phase_fac, l1, l2, l3
	
	l1 = int(wigin(1,1))
	l2 = int(wigin(1,2))
	l3 = int(wigin(1,3))
	
	if(mod(l1+l2+l3,2) == 0) then
		phase_fac = 1
	else
		if(mod(trans1 + trans2,2) == 0)then
			phase_fac = 1
		else
			phase_fac = -1
		end if
	end if
	
	end function phase_fac

!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				This is the most important subroutine
!				Calculates and stores Wigner Symbols & passes the array
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine wignstore(l_min, l_max, max_ind, wigarr)
	
	implicit none
	
	integer, intent(in) :: max_ind, l_max, l_min
	double precision, intent(out), dimension(:,:), allocatable :: wigarr
	integer :: len_max, counter, l1, l2, l3, m1, m2, m3, index_loc, no_trans
	real, dimension(2,3) :: wigin
	real :: eps
	
	allocate(wigarr(1:max_ind,1:2))
	
!	print*, 'Computing & Storing Wigner Symbols...'
	counter = 0
	eps = 1.0E-20
	wigarr = 0.0

!	calculation of Wigner3j symbols as per the algorithm given in Rasch & Yu SIAM 25, 1416
	do l1 = l_min, l_max + ecare
		do l2 = l_min, l_max + ecare
			do l3 = abs(l1-l2), l1+l2 !ensure triangular condition
				do m1 = -l1, l1
					do m2 = -l2, l2
						m3 = -m1-m2
						if(abs(m3) <= l3) then
							call par2wig(1.0*l1, 1.0*l2, 1.0*l3, 1.0*m1, 1.0*m2, 1.0*m3, wigin)
							call wignstandindex(wigin, index_loc, no_trans)
							
							if(index_loc > max_index) then
								print*, 'Not enough Wigner Symbols Stored'
								stop
							end if
							
							if((abs(wigarr(index_loc,1)) == 0.0))then
!	this is an improved way, I've just added eps to the thing
								wigarr(index_loc,1) = wigner3j(wigin) + eps
								wigarr(index_loc,2) = no_trans
								counter = counter + 1
							end if
						end if
					end do
				end do
			end do
		end do
	end do
	
	print*, 'Stored Symbols:  ', counter
	
	end subroutine wignstore
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Retrieval part of the Wigner Symbols
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine wignretr(l1, l2, l3, m1, m2, m3, wigarr, max_ind, ans)
	
	implicit none

	integer, intent(in) :: l1, l2, l3, m1, m2, m3
	integer, intent(in) :: max_ind
	double precision, intent(in), dimension(max_ind,2) :: wigarr
	real, dimension(2,3) :: wigin
	double precision :: ans
	integer :: index_loc, no_trans, phase
	
	wigin(1,1) = 1.0*l1
	wigin(1,2) = 1.0*l2
	wigin(1,3) = 1.0*l3
	wigin(2,1) = 1.0*m1
	wigin(2,2) = 1.0*m2
	wigin(2,3) = 1.0*m3

	call wignstandindex(wigin, index_loc, no_trans)
	phase = phase_fac(wigin, no_trans, int(wigarr(index_loc,2)))
	ans = wigarr(index_loc,1)*phase
	
	end subroutine wignretr
	
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!				Convert seconds to hours
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	subroutine elaptime(time, hour, mins, seconds)
	
	implicit none
	
	real, intent(in) :: time
	integer :: hour, mins
	real :: seconds
	hour = int(time/3600.0)
	mins = int((time - hour*3600.0)/60.0)
	seconds = time - hour*3600.0 - mins*60.0

	end subroutine elaptime
!	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	end module wig_sym
