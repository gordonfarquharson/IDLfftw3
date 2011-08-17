/* 
   IDLfftw3

   Copyright (C) 2005 Paco López Dekker <paco.dekker@gmail.com>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
      
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>
#include <idl_export.h>

#define MAXDIM 4

#define IDL_FFTW_NONE 0
#define IDL_FFTW_C2C  1
#define IDL_FFTW_R2C  2
#define IDL_FFTW_C2R  3

#define PLAN_INFO_LENGTH 9

IDL_VPTR idlfftw_plan(int argc, IDL_VPTR argv[], char *argk)
{
	IDL_VPTR datain, dataout = NULL, IDL_plan;
	int ndim, dim_index, dims[MAXDIM], dimsout[MAXDIM];
	char infomessage[200];
	fftwf_complex *c_in, *c_out;
	float *f_in, *f_out;
	/* A plan is an opaque pointer to a structure. I am going to
	   store the pointer as an unsigned long, this may not scale to
	   64-bit machines. This may be a mistake though */
	fftwf_plan plan = 0;
	IDL_ULONG *planinfo;

	/* Stuff used for keyword processing */
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;	/* Must be first entry in this structure */
		IDL_LONG c2r;
		IDL_LONG dimension;
		IDL_LONG fftw_backward;
		IDL_LONG lazy;
		IDL_LONG fftw_test;
		IDL_LONG verbose;
	} KW_RESULT;
	static IDL_KW_PAR kw_pars[] = {
		{"C2R", IDL_TYP_LONG, 1, IDL_KW_ZERO,
		 0, IDL_KW_OFFSETOF(c2r)},
		{"DIMENSION", IDL_TYP_LONG, 1, IDL_KW_ZERO,
		 0, IDL_KW_OFFSETOF(dimension)},
		{"INVERSE", IDL_TYP_LONG, 1, IDL_KW_ZERO,
		 0, IDL_KW_OFFSETOF(fftw_backward)},
		{"LAZY", IDL_TYP_LONG, 1, IDL_KW_ZERO,
		 0, IDL_KW_OFFSETOF(lazy)},
		{"TEST", IDL_TYP_LONG, 1, IDL_KW_ZERO,
		 0, IDL_KW_OFFSETOF(fftw_test)},
		{"VERBOSE", IDL_TYP_LONG, 1, IDL_KW_ZERO,
		 0, IDL_KW_OFFSETOF(verbose)},
		{NULL}
	};
	KW_RESULT kw;

	/* Execution flags */
	unsigned verbose = 0,
	    fftw_dir = FFTW_FORWARD,
	    fftw_dimension = 0,
	    fftw_type = IDL_FFTW_NONE,
	    fftw_c2r = 0, fftw_test = 0, planner_flag = FFTW_MEASURE;

	/* Keyword processing */
	argc = IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);
	if (kw.fftw_backward)
		fftw_dir = FFTW_BACKWARD;
	if (kw.dimension)
		fftw_dimension = kw.dimension;
	if (kw.lazy)
		planner_flag = planner_flag | FFTW_UNALIGNED;
	if (kw.fftw_test)
		fftw_test = 1;
	if (kw.verbose) {
		verbose = 1;
		IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO, "Verbose.");
	}
	if (kw.c2r)
		fftw_c2r = 1;
	IDL_KW_FREE;

	planinfo =
	    (IDL_ULONG *) IDL_MakeTempVector(IDL_TYP_ULONG,
					     (IDL_MEMINT) PLAN_INFO_LENGTH,
					     IDL_ARR_INI_ZERO, &IDL_plan);

	datain = argv[0];
	IDL_ENSURE_ARRAY(datain);

	/* See what kind of variable we have */
	/* Lets first see how many dimensions */
	ndim = datain->value.arr->n_dim;
	if (ndim > MAXDIM) {
		sprintf(infomessage,
			"Input array shoud have at most %i dimensions.\n",
			(int)MAXDIM);
		IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, infomessage);
	}
	if (fftw_dimension > ndim)
		IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
			    "Trying to plan a fftw over a non-existing dimension.\n");
	if ((fftw_dimension != 0) && (fftw_dimension != 1)
	    && (fftw_dimension != ndim))
		IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
			    "At this time the dimension keyword only works for the first or last dimension.\n");
	if (argc > 1) {
		dataout = argv[1];
		IDL_ENSURE_ARRAY(dataout);
	}
	for (dim_index = 0; dim_index < ndim; dim_index++) {
		/*  arrays are stored in row format or whatever in IDL */
		dims[dim_index] = datain->value.arr->dim[ndim - 1 - dim_index];
		if (argc > 1)
			dimsout[dim_index] =
			    dataout->value.arr->dim[ndim - 1 - dim_index];
	}
	/*  Determine what kind of transform we have to do */
	switch (datain->type) {
	case IDL_TYP_FLOAT:
		fftw_type = IDL_FFTW_R2C;
		break;
	case IDL_TYP_COMPLEX:
		fftw_type = IDL_FFTW_C2C;
		if (fftw_c2r)
			fftw_type = IDL_FFTW_C2R;
		break;
	default:
		fftw_type = IDL_FFTW_NONE;
	}
	switch (fftw_type) {
	case IDL_FFTW_R2C:
		if (verbose)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
				    "We have an array of floats.\n");
		if (argc == 1)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
				    "Author's note: idlfftw will not deal with in place real to complex fftw.\n");
		f_in = (float *)datain->value.arr->data;
		if (dataout->type != IDL_TYP_COMPLEX)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
				    "Second argument should be a complex array.\n");
		/* Check size of output array */
		{
			int critical_dimension = ndim - 1;
			int wrong_size = 0;
			if (fftw_dimension)
				critical_dimension = ndim - fftw_dimension;
			for (dim_index = 0; dim_index < ndim - 1; dim_index++) {
				if (dim_index == critical_dimension)
					wrong_size = wrong_size
					    || ((floor(dims[dim_index] / 2) + 1)
						!= dimsout[dim_index]);
				else
					wrong_size = wrong_size
					    || (dims[dim_index] !=
						dimsout[dim_index]);
			}
			/* dims[ndim-1] = (floor(dims[ndim-1]/2)+1); */
			if (wrong_size)
				IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
					    "Dimensions of output array don't match.\n");
		}
		c_out = (fftwf_complex *) dataout->value.arr->data;
		if (fftw_dimension == 0) {
			plan =
			    fftwf_plan_dft_r2c(ndim, dims, f_in, c_out,
					       planner_flag);
		} else {
			int howmany =
			    (int)(datain->value.arr->n_elts /
				  dims[ndim - fftw_dimension]);
			if (verbose) {
				sprintf(infomessage, "Will do %i FFTs.\n",
					howmany);
				IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
					    infomessage);
			}
			int istride = 1, idist = 1, ostride, odist;
			if (fftw_dimension == 1) {
				idist = dims[ndim - 1];
				odist = dimsout[ndim - 1];
				ostride = istride = 1;

			} else {
				odist = idist = 1;
				istride = dims[ndim - 1];
				ostride = dimsout[ndim - 1];
			}
			plan =
			    fftwf_plan_many_dft_r2c(1,
						    &(dims
						      [ndim - fftw_dimension]),
						    howmany, f_in, NULL,
						    istride, idist, c_out, NULL,
						    ostride, odist,
						    planner_flag);
		}
		break;
	case IDL_FFTW_C2R:
		if (verbose)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
				    "We have an array of complex.\n");
		if (argc == 1)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
				    "Author's note: idlfftw will not deal with in place complex to real fftw.\n");
		c_in = (fftwf_complex *) datain->value.arr->data;
		if (dataout->type != IDL_TYP_FLOAT)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
				    "Second argument should be a float array.\n");
		/* Check size of output array */
		{
			int critical_dimension = ndim - 1;
			int wrong_size = 0;
			if (fftw_dimension)
				critical_dimension = ndim - fftw_dimension;
			for (dim_index = 0; dim_index < ndim - 1; dim_index++) {
				if (dim_index == critical_dimension)
					wrong_size = wrong_size
					    ||
					    ((floor(dimsout[dim_index] / 2) +
					      1) != dims[dim_index]);
				else
					wrong_size = wrong_size
					    || (dims[dim_index] !=
						dimsout[dim_index]);
			}
			/* dims[ndim-1] = (floor(dims[ndim-1]/2)+1); */
			if (wrong_size)
				IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
					    "Dimensions of output array don't match.\n");
		}
		f_out = (float *)dataout->value.arr->data;
		if (fftw_dimension == 0) {
			plan =
			    fftwf_plan_dft_c2r(ndim, dimsout, c_in, f_out,
					       planner_flag);
		} else {
			int howmany =
			    (int)(datain->value.arr->n_elts /
				  dims[ndim - fftw_dimension]);
			int istride = 1, idist = 1, odist, ostride;
			if (fftw_dimension == 1) {
				idist = dims[ndim - 1];
				odist = dimsout[ndim - 1];
				ostride = istride = 1;
			} else {
				odist = idist = 1;
				istride = dims[ndim - 1];
				ostride = dimsout[ndim - 1];
			}
			plan =
			    fftwf_plan_many_dft_c2r(1,
						    &(dimsout
						      [ndim - fftw_dimension]),
						    howmany, c_in, NULL,
						    istride, idist, f_out, NULL,
						    ostride, odist,
						    planner_flag);
		}
		break;
	case IDL_FFTW_C2C:
		if (verbose)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
				    "We have an array of complex.\n");
		c_out = c_in = (fftwf_complex *) datain->value.arr->data;
		if (argc > 1) {
			/* Check that the output array is valid */
			if (dataout->type != IDL_TYP_COMPLEX)
				IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
					    "Second argument should also be a complex array.\n");
			if (datain->value.arr->n_elts !=
			    dataout->value.arr->n_elts)
				IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
					    "Second argument should have same size as first argument.\n");
			c_out = (fftwf_complex *) dataout->value.arr->data;
		}
		if (fftw_dimension == 0) {
			plan =
			    fftwf_plan_dft(ndim, dims, c_in, c_out, fftw_dir,
					   planner_flag);
		} else {
			int howmany;
			int istride, idist;
			howmany =
			    (int)((datain->value.arr->n_elts) /
				  dims[ndim - fftw_dimension]);
			if (verbose) {
				sprintf(infomessage, "Will do %i FFTs.\n",
					howmany);
				IDL_Message(IDL_M_GENERIC, IDL_MSG_INFO,
					    infomessage);
			}
			if (fftw_dimension == 1) {
				idist = dims[ndim - 1];
				istride = 1;
			} else {
				idist = 1;
				istride = dims[ndim - 1];
			}
			plan =
			    fftwf_plan_many_dft(1,
						&(dims[ndim - fftw_dimension]),
						howmany, c_in, NULL, istride,
						idist, c_out, NULL, istride,
						idist, fftw_dir, planner_flag);

		}
		break;
	default:
		IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
			    "Not a valid type: Input array has to be float or complex.\n");
	}

	if (fftw_test)
		fftwf_execute(plan);

	/* Now it is time to pupulate the planinfo array */
	{
		planinfo[0] = (IDL_ULONG) plan;	/*  This is the plan */
		if (argc == 1) {
			planinfo[1] = 1;
			planinfo[2] = (IDL_ULONG) datain->value.arr->data;	/* A pointer to the first array; */
			planinfo[3] = (IDL_ULONG) datain->value.arr->arr_len;	/* Physical size o the array in bytes.
										   This is useful to prevent segfaults, I hope */
		} else {
			planinfo[1] = 2;
			planinfo[2] = (IDL_ULONG) datain->value.arr->data;	/* A pointer to the first array; */
			planinfo[3] = (IDL_ULONG) datain->value.arr->arr_len;
			planinfo[4] = (IDL_ULONG) dataout->value.arr->data;	/* A pointer to the first array; */
			planinfo[5] = (IDL_ULONG) dataout->value.arr->arr_len;
		}
		planinfo[6] = 13579;	/* these two are useless, but I will use them to check that it is a plan */
		planinfo[7] = 2468;
		planinfo[8] = fftw_type;	/*This is so that we can execute the right function if playing with the guru mode */
	}
	return (IDL_plan);

}

int idlfftw_checkplan(IDL_VPTR IDL_plan)
{
	int badplan = 0;
	IDL_ULONG *planinfo;

	IDL_ENSURE_ARRAY(IDL_plan);

	if (IDL_plan->type != IDL_TYP_ULONG)
		badplan = 1;
	if (IDL_plan->value.arr->n_elts != PLAN_INFO_LENGTH)
		badplan = 2;
	if (badplan != 1) {
		planinfo = (IDL_ULONG *) IDL_plan->value.arr->data;
		if ((planinfo[6] != 13579) || (planinfo[7] != 2468)
		    || (planinfo[0] == 0))
			badplan = 3;
	}
	return (badplan);
}

void idlfftw_delplan(int argc, IDL_VPTR argv[])
{
	IDL_VPTR IDL_plan;
	fftwf_plan plan;
	IDL_ULONG *planinfo;
	IDL_plan = argv[0];
	int error;

	/* Let's check that it may be a plan. This may still easily go
	   wrong and cause a segmentation fault or something evil like
	   that. */
	if ((error = idlfftw_checkplan(IDL_plan))) {
		char error_msg[128];
		sprintf(error_msg, "Not a valid plan so let's not try this (%d).\n",
			error);
		IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, error_msg);
	}
	planinfo = (IDL_ULONG *) IDL_plan->value.arr->data;
	plan = (fftwf_plan) planinfo[0];
	fftwf_destroy_plan(plan);
	/* Now a very basic security measure to identify it as a
	   deleted plan */
	planinfo[0] = 0;
}

void idlfftw(int argc, IDL_VPTR argv[], char *argk)
{
	IDL_VPTR IDL_plan, IDL_in, IDL_out = NULL;
	fftwf_plan plan;
	IDL_ULONG *planinfo, fftw_type;
	int bad_data = 0, gurumode = 0;
	int error;

	/*  Stuff used for keyword processing */
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;	/* Must be first entry in this structure */
		IDL_LONG guru;
	} KW_RESULT;
	static IDL_KW_PAR kw_pars[] = {
		{"GURU", IDL_TYP_LONG, 1, IDL_KW_ZERO,
		 0, IDL_KW_OFFSETOF(guru)},
		{NULL}
	};
	KW_RESULT kw;

	argc = IDL_KWProcessByOffset(argc, argv, argk, kw_pars, NULL, 1, &kw);
	if (kw.guru)
		gurumode = 1;
	IDL_KW_FREE;
	if (argc < 2)
		IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
			    "IDLfftw requires at least a plan and an array.\n");
	/* Retrieve the plan */
	IDL_plan = argv[0];
	/* Let's check that it may be a plan. This may still easily go
	   wrong and cause a segmentation fault or something evil like
	   that. */
	if ((error = idlfftw_checkplan(IDL_plan))) {
		char error_msg[128];
		sprintf(error_msg, "Not a valid plan so let's not try this (%d).\n",
			error);
		IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP, error_msg);
	}
	planinfo = (IDL_ULONG *) IDL_plan->value.arr->data;
	if (gurumode == 0) {
		/* This is how I think it should be done, just use the
		   same arrays as where given to the planner */
		IDL_in = argv[1];
		/* Check that the arrays are where expected and have
		   the right size */
		{
			IDL_ENSURE_ARRAY(IDL_in);
			if (planinfo[2] != (IDL_ULONG) IDL_in->value.arr->data)
				bad_data = 1;	/* A pointer to the first array; */
			if (planinfo[3] !=
			    (IDL_ULONG) IDL_in->value.arr->arr_len)
				bad_data = 1;	/* Physical size o the array in bytes. */
		}

		if (argc > 2) {
			IDL_out = argv[2];
			IDL_ENSURE_ARRAY(IDL_out);
			if (planinfo[1] != 2)
				bad_data = 1;
			if (planinfo[4] != (IDL_ULONG) IDL_out->value.arr->data)
				bad_data = 1;	/* A pointer to the first array; */
			if (planinfo[5] !=
			    (IDL_ULONG) IDL_out->value.arr->arr_len)
				bad_data = 1;	/* Physical size o the array in bytes. */
		} else {
			if (planinfo[1] != 1)
				bad_data = 1;
		}
		if (bad_data)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
				    "Arrays do not correspond to plan: this would end in a segmentation fault or worse.\n");
		plan = (fftwf_plan) planinfo[0];
		fftwf_execute(plan);
	} else {
		/* In guru mode we can reuse a plan with different
		   arrays. Depending on how IDL aligns the arrays in
		   memory this may work well or be a huge disaster. In
		   principle, if the plan is created with the lazy
		   flag set it should be safe to do this (then you are
		   a lazy guru ;-) */
		IDL_in = argv[1];
		/* Check that the arrays are have the right size */
		{
			IDL_ENSURE_ARRAY(IDL_in);
			if (planinfo[3] !=
			    (IDL_ULONG) IDL_in->value.arr->arr_len)
				bad_data = 1;	/* Physical size o the array in bytes. */
		}
		if (argc > 2) {
			IDL_out = argv[2];
			IDL_ENSURE_ARRAY(IDL_out);
			if (planinfo[1] != 2)
				bad_data = 1;
			if (planinfo[5] !=
			    (IDL_ULONG) IDL_out->value.arr->arr_len)
				bad_data = 1;	/* Physical size o the array in bytes. */
		} else {
			if (planinfo[1] != 1)
				bad_data = 1;
		}
		if (bad_data)
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
				    "Arrays do not correspond to plan: this would end in a segmentation fault or worse.\n");
		plan = (fftwf_plan) planinfo[0];
		fftw_type = planinfo[8];
		switch (fftw_type) {
		case IDL_FFTW_C2C:
			if (argc > 2)
				fftwf_execute_dft(plan,
						  (fftwf_complex *)
						  IDL_in->value.arr->data,
						  (fftwf_complex *)
						  IDL_out->value.arr->data);
			else
				fftwf_execute_dft(plan,
						  (fftwf_complex *)
						  IDL_in->value.arr->data,
						  (fftwf_complex *)
						  IDL_in->value.arr->data);
			break;
		case IDL_FFTW_R2C:
			fftwf_execute_dft_r2c(plan,
					      (float *)IDL_in->value.arr->data,
					      (fftwf_complex *) IDL_out->
					      value.arr->data);
			break;
		case IDL_FFTW_C2R:
			fftwf_execute_dft_c2r(plan,
					      (fftwf_complex *) IDL_in->
					      value.arr->data,
					      (float *)IDL_out->value.
					      arr->data);
			break;
		default:
			IDL_Message(IDL_M_GENERIC, IDL_MSG_LONGJMP,
				    "Invalid fftw request (this should never happen).\n");
		}
	}
}

int IDL_Load(void)
{
	/* These tables contain information on the functions and
	 * procedures that make up the DLM. The information
	 * contained in these tables must be identical to that
	 * contained in the .dlm file.
	 */
	static IDL_SYSFUN_DEF2 function_addr[] = {
		{ { idlfftw_plan }, 
		  "IDLFFTW_PLAN", 1, 2, IDL_SYSFUN_DEF_F_KEYWORDS, 0 },
	};
	static IDL_SYSFUN_DEF2 procedure_addr[] = {
		{ { (IDL_SYSRTN_GENERIC) idlfftw_delplan },
		  "IDLFFTW_DELPLAN", 1, 1, 0, 0 },
		{ { (IDL_SYSRTN_GENERIC) idlfftw },
		  "IDLFFTW", 1, 3, IDL_SYSFUN_DEF_F_KEYWORDS, 0},
	};

	/* Register our routine. The routines must be specified
	 * exactly the same as in testmodule.dlm. */
	return IDL_SysRtnAdd(function_addr, TRUE, 
			     IDL_CARRAY_ELTS(function_addr)) &&
	       IDL_SysRtnAdd(procedure_addr, FALSE,
		             IDL_CARRAY_ELTS(procedure_addr));
}
