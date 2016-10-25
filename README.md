# KrATER (K-mer Analysis Tool ER)

I. Requirements

    1. Python libraries:
        Matplotlib
        Numpy
        Scipy
        
    2. Jellyfish 2

II. Installation

    1. Install requirements
        - run following command for global installation
            sudo pip install matplotlib numpy scipy 
        - or following command for local installation if you don't have root permissions
            pip install --user matplotlib numpy scipy

    2. Install KrATER
        Variant 1: install using pip
        
                TODO: add to pip
            
        Variant 2: install from source code with root permissions
        
            git clone https://github.com/mahajrod/krater
            cd KrATER
            python setup.py build
            sudo python setup.py install
        
        
        Variant 3: install from source code without installation
        
            Get KRATER
                git clone https://github.com/mahajrod/krater
        
            Add following strings to ~/.profile and ~/.bashrc (create files if absent). Don't forget to replace <KRATER_DIR> with actual path
            
                PYTHONPATH=${PYTHONPATH}:<KRATER_DIR>
                export PYTHONPATH
    
            Run in terminal
                source ~/.profile
            
    
III. RUN

    1. Input file/files is/are fastq:
    
        Typical usage:
            <path to KRATER dir>/scripts/kmer/draw_kmer_distribution_from_fastq.py -m 23  -t ${THREAD_NUMBER} -b -s 30G -e png -o ${OUTPUT_PREFIX} -i {COMMA_SEPARATED_LIST_OF_FASTQ_FILES} -w ${MINIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE}  -g ${MAXIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE}
        
        Parameter description:
            
        
    2. Input file is histogram file produced from jellyfish
    
        <path to KRATER dir>/draw_kmer_distribution_from_histo.py  -i ${JELLYFISH_HISTO_FILE}  -o ${OUTPUT_PREFIX} -w ${MINIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE} -g ${MAXIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE} -e png
        
        Parameter_description:
        