Module kinds_module
integer,parameter:: kr8=
integer,parameter:: kr4=
integer,parameter:: ki4=
integer,parameter:: kL4=
integer,parameter:: ki8=



! ieee 754 standard defined single and double precision reals
integer,parameter :: 

    ! Kind types for 64-, 32-, 16-, and 8-bit signed integers
    integer, parameter :: INT64 = selected_int_kind(18)
    integer, parameter :: INT32 = selected_int_kind(9)
    integer, parameter :: INT16 = selected_int_kind(4)
    integer, parameter :: INT08 = selected_int_kind(2)

    ! Kind types for IEEE 754/IEC 60559 single- and double-precision reals
    integer, parameter :: IEEE32 = selected_real_kind(  6,  37 )
    integer, parameter :: IEEE64 = selected_real_kind( 15, 307 )


end module kinds_module
