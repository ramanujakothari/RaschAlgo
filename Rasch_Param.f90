!	this is the parameter file that contains all the global variables to
!	to be used with the code
	module parameters
	
	use healpix_types
	
	integer(i4b), parameter :: l_min = 2
	integer(i4b), parameter :: l_max = 16
	integer(i4b), parameter :: ecare = 2
!	extra care variable to ensure enough wigner symbols are taken

!	whether or not one needs to print the number of symbols	
	logical :: want_nsymb = .false.
	end module parameters
