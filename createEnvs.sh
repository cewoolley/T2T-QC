# Create environments w/ required tools. ncrf needs an older version of Python so gets its own venv.

mamba create -n t2t-assess -c bioconda -c conda-forge \
        python=3.10 \
        seqtk=1.3 \
        gawk=5.1.0 \
        bioawk=1.0 \
        ncrf=1.01.02 \
        mashmap=3.1.3 \
        coreutils=9.5 \
        blis \
        -y

mamba create -n ncrf-py2 -c bioconda -c conda-forge \
        python=2.7 \
        ncrf \
        -y
