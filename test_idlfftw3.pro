PRO test_idlfftw3

    a = [1.,0.,1.,0.]
    b = complexarr(4)

    plan = idlfftw_plan(a, b, /VERBOSE)

    help, plan
    print, plan

    idlfftw, plan, a, b

    print, a
    print, b
    print, fft(a)

END
