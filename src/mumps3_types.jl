# this file mirrors the "mumps_c_types.h" header file of MUMPS 5.2.2

const MUMPS_INT = Cint
const MUMPS_INT8 = Int64

const SMUMPS_COMPLEX = Cfloat
const SMUMPS_REAL = Cfloat

const DMUMPS_COMPLEX = Cdouble
const DMUMPS_REAL = Cdouble

mutable struct Mumps_complex r::Cfloat; i::Cfloat end
mutable struct Mumps_double_complex r::Cdouble; i::Cdouble end
const CMUMPS_COMPLEX = Mumps_complex
const CMUMPS_REAL = Cfloat

const ZMUMPS_COMPLEX = Mumps_double_complex
const ZMUMPS_REAL = Cdouble

const mumps_ftnlen = MUMPS_INT

const MUMPS_ARITH_s = 1
const MUMPS_ARITH_d = 2
const MUMPS_ARITH_c = 4
const MUMPS_ARITH_z = 8

const MUMPS_ARITH_REAL = ( MUMPS_ARITH_s | MUMPS_ARITH_d )
const MUMPS_ARITH_CMPLX = ( MUMPS_ARITH_c | MUMPS_ARITH_z )
const MUMPS_ARITH_SINGLE = ( MUMPS_ARITH_s | MUMPS_ARITH_c )
const MUMPS_ARITH_DBL = ( MUMPS_ARITH_d | MUMPS_ARITH_z )
