all: rrc_new rrc_cum rrc_abmon rrc_timevarying_method_GENEnew

rrc_new: rrc_new.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f
	gfortran -lm -o rrc_new rrc_new.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f

rrc_cum: rrc_cum.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f
	gfortran -lm -o rrc_cum rrc_cum.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f

rrc_abmon: rrc_abmon.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f
	gfortran -lm -o rrc_abmon rrc_abmon.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f

rrc_timevarying_method_GENEnew: rrc_timevarying_method_GENEnew.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f
	gfortran -lm -o rrc_timevarying_method_GENEnew rrc_timevarying_method_GENEnew.f quicksort.f90 lapack.f chidf.f fzero.f qld.f ffsqp.f blas.f
