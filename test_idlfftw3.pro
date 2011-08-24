PRO test_idlfftw3

    t0 = systime(/SECONDS)

    ;; Test 2-d complex array with separate input and output variables
    mat = complexarr(64, 32)
    mat2 = complexarr(64, 32)
    plan = idlfftw_plan(mat, mat2, /INVERSE, /LAZY, RANK=2)
    mat = make_array(64, 32, /COMPLEX, VALUE=complex(1.,0.))
    idlfftw, plan, mat, mat2, /GURU
    idlfftw_delplan, plan

    ;; Test in-place 2-d FFT of a complex array
    mat = complexarr(64, 32)
    plan = idlfftw_plan(mat, /INVERSE, /LAZY, RANK=2)
    mat = make_array(64, 32, /COMPLEX, VALUE=complex(1.,0.))
    idlfftw, plan,mat, /GURU
    idlfftw_delplan, plan

    ;; Test multiple 2-d FFTs of a 4-d matrix
    I = [4,4]
    Na = 64
    Nr = 64
    Mr = 32
    Ma = 32
    inv_corr_zp = complexarr(I[0]*Nr, I[1]*Na, Mr, Ma)
    ;;inv_corr_zp_it = complexarr(I[0]*Nr, I[1]*Na, Mr, Ma)
    fftw_plan = idlfftw_plan(inv_corr_zp, /INVERSE, /LAZY, RANK=2)
    inv_corr_zp[0:Mr-1,0:Mr-1,0:Mr-1:2,*] = complex(1.,0.)
    idlfftw, fftw_plan, inv_corr_zp, /GURU
    ;;inv_corr_zp_it = 0
    ;;first_fft = 0

    IF NOT keyword_set(silent) THEN print, "Total time: ", $
                                           systime(/SECONDS) - t0

END
