    INTEGER :: FFTW_R2HC
    PARAMETER (FFTW_R2HC=0)
    INTEGER :: FFTW_HC2R
    PARAMETER (FFTW_HC2R=1)
    INTEGER :: FFTW_DHT
    PARAMETER (FFTW_DHT=2)
    INTEGER :: FFTW_REDFT00
    PARAMETER (FFTW_REDFT00=3)
    INTEGER :: FFTW_REDFT01
    PARAMETER (FFTW_REDFT01=4)
    INTEGER :: FFTW_REDFT10
    PARAMETER (FFTW_REDFT10=5)
    INTEGER :: FFTW_REDFT11
    PARAMETER (FFTW_REDFT11=6)
    INTEGER :: FFTW_RODFT00
    PARAMETER (FFTW_RODFT00=7)
    INTEGER :: FFTW_RODFT01
    PARAMETER (FFTW_RODFT01=8)
    INTEGER :: FFTW_RODFT10
    PARAMETER (FFTW_RODFT10=9)
    INTEGER :: FFTW_RODFT11
    PARAMETER (FFTW_RODFT11=10)
    INTEGER :: FFTW_FORWARD
    PARAMETER (FFTW_FORWARD=-1)
    INTEGER :: FFTW_BACKWARD
    PARAMETER (FFTW_BACKWARD=+1)
    INTEGER :: FFTW_MEASURE
    PARAMETER (FFTW_MEASURE=0)
    INTEGER :: FFTW_DESTROY_INPUT
    PARAMETER (FFTW_DESTROY_INPUT=1)
    INTEGER :: FFTW_UNALIGNED
    PARAMETER (FFTW_UNALIGNED=2)
    INTEGER :: FFTW_CONSERVE_MEMORY
    PARAMETER (FFTW_CONSERVE_MEMORY=4)
    INTEGER :: FFTW_EXHAUSTIVE
    PARAMETER (FFTW_EXHAUSTIVE=8)
    INTEGER :: FFTW_PRESERVE_INPUT
    PARAMETER (FFTW_PRESERVE_INPUT=16)
    INTEGER :: FFTW_PATIENT
    PARAMETER (FFTW_PATIENT=32)
    INTEGER :: FFTW_ESTIMATE
    PARAMETER (FFTW_ESTIMATE=64)
    INTEGER :: FFTW_ESTIMATE_PATIENT
    PARAMETER (FFTW_ESTIMATE_PATIENT=128)
    INTEGER :: FFTW_BELIEVE_PCOST
    PARAMETER (FFTW_BELIEVE_PCOST=256)
    INTEGER :: FFTW_NO_DFT_R2HC
    PARAMETER (FFTW_NO_DFT_R2HC=512)
    INTEGER :: FFTW_NO_NONTHREADED
    PARAMETER (FFTW_NO_NONTHREADED=1024)
    INTEGER :: FFTW_NO_BUFFERING
    PARAMETER (FFTW_NO_BUFFERING=2048)
    INTEGER :: FFTW_NO_INDIRECT_OP
    PARAMETER (FFTW_NO_INDIRECT_OP=4096)
    INTEGER :: FFTW_ALLOW_LARGE_GENERIC
    PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
    INTEGER :: FFTW_NO_RANK_SPLITS
    PARAMETER (FFTW_NO_RANK_SPLITS=16384)
    INTEGER :: FFTW_NO_VRANK_SPLITS
    PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
    INTEGER :: FFTW_NO_VRECURSE
    PARAMETER (FFTW_NO_VRECURSE=65536)
    INTEGER :: FFTW_NO_SIMD
    PARAMETER (FFTW_NO_SIMD=131072)
    INTEGER :: FFTW_NO_SLOW
    PARAMETER (FFTW_NO_SLOW=262144)
    INTEGER :: FFTW_NO_FIXED_RADIX_LARGE_N
    PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
    INTEGER :: FFTW_ALLOW_PRUNING
    PARAMETER (FFTW_ALLOW_PRUNING=1048576)
    INTEGER :: FFTW_WISDOM_ONLY
    PARAMETER (FFTW_WISDOM_ONLY=2097152)
