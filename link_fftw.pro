PRO link_fftw

    linkimage, 'idlfftw_plan', 'libIDLfftw.so', /FUN, /KEYWORD
    linkimage, 'idlfftw_delplan', 'libIDLfftw.so'
    linkimage, 'idlfftw', 'libIDLfftw.so', /KEYWORD

END
