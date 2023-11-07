# this file contains a bunch of shortcuts for manipulating ICNTL.

export set_error_stream!, set_diagnostics_stream!, set_info_stream!, set_print_level!,
suppress_printing!, toggle_printing!, suppress_display!, toggle_display!,
sparse_matrix!, dense_matrix!,
sparse_rhs!, dense_rhs!,
toggle_null_pivot!

const ICNTL_DEFAULT = (6, 0, 6, 2, 0, 7, 7, 77, 1, 0, 0, 1, 0, 20, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -32, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 333, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0)

"""
    default_icntl!(mumps)

reset ICNTL to its default
"""
function default_icntl!(mumps::Mumps)
    mumps.icntl = ICNTL_DEFAULT
    return nothing
end

set_error_stream!(mumps::Mumps,i) = set_icntl!(mumps,1,i; displaylevel=0)
set_diagnostics_stream!(mumps::Mumps,i) = set_icntl!(mumps,2,i; displaylevel=0)
set_info_stream!(mumps::Mumps,i) = set_icntl!(mumps,3,i; displaylevel=0)
set_print_level!(mumps::Mumps,i) = set_icntl!(mumps,4,i; displaylevel=0)

suppress_printing!(mumps::Mumps) = set_print_level!(mumps,1)
toggle_printing!(mumps::Mumps) = set_print_level!(mumps,mod1(mumps.icntl[4]+1,2))
suppress_display! = suppress_printing!
toggle_display! = toggle_printing!

sparse_matrix!(mumps::Mumps) = set_icntl!(mumps,5,0; displaylevel=0)
dense_matrix!(mumps::Mumps) = set_icntl!(mumps,5,1; displaylevel=0)

LinearAlgebra.transpose!(mumps::Mumps) = set_icntl!(mumps,9,mod(mumps.icntl[9]+1,2); displaylevel=0)

sparse_rhs!(mumps::Mumps) = set_icntl!(mumps, 20, 1; displaylevel=0)
dense_rhs!(mumps::Mumps) = set_icntl!(mumps, 20, 0; displaylevel=0)

toggle_null_pivot!(mumps::Mumps) = set_icntl!(mumps,24,mod(mumps.icntl[24]+1,2); displaylevel=0)
