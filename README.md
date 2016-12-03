# KrATER (K-mer Analysis Tool ER)

I. Requirements

    1. Python libraries:
        Matplotlib
        Numpy
        Scipy
        
    2. Jellyfish 2

II. Installation

    1. Install requirements (not necessary if you going to install KrATER via pip)
        - run following command for global installation
        
            sudo pip install matplotlib numpy scipy 
            
        - or following command for local installation if you don't have root permissions
        
            pip install --user matplotlib numpy scipy
            
    2. Install KrATER
        Variant 1: install using pip
        
           - run following command for global installation
           
                sudo pip install krater
           
           - run following commands for local installation
                
                pip install --user krater
                
              Then add ~/.local/bin in your PATH variable
              
                cat "PATH=\${PATH}:~/.local/bin/" >> ~/.profile
                cat "PATH=\${PATH}:~/.local/bin/" >> ~/.bashrc
                
              Load updated PATH variable
              
                source ~/.profile
            
        Variant 2: install from source code with root permissions
        
            git clone https://github.com/mahajrod/krater
            cd krater
            python setup.py build
            sudo python setup.py install
        
        
        Variant 3: install from source code without root permissions
        
            Get KRATER
                git clone https://github.com/mahajrod/krater
        
            Add following strings to ~/.profile and ~/.bashrc (create files if absent). Don't forget to replace <KRATER_DIR> with actual path
            
                PYTHONPATH=${PYTHONPATH}:<KRATER_DIR>
                PATH=${PATH}:<KRATER_DIR>
                export PYTHONPATH
                export PATH
                
            Run in terminal
                source ~/.profile
            
    
III. RUN

    1. If input file/files is/are fastq:
    
        Typical usage:
            draw_kmer_distribution_from_fastq.py -m 23  -t ${THREAD_NUMBER} -b -s 30G -e png -o ${OUTPUT_PREFIX} -i {COMMA_SEPARATED_LIST_OF_FASTQ_FILES} -w ${MINIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE}  -g ${MAXIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE}
        
        Parameter description:
        
    2. If input file is jellyfish database
        
        Typical usage:
            draw_kmer_distribution_from_jellyfish_database.py -i ${JELLYFISH_DATABASE}  -o ${OUTPUT_PREFIX} -w ${MINIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE} -g ${MAXIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE} -e png
        
    2. If input file is histogram file produced by Jellyfish:
    
        Typical usage
            draw_kmer_distribution_from_histo.py  -i ${JELLYFISH_HISTO_FILE}  -o ${OUTPUT_PREFIX} -w ${MINIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE} -g ${MAXIMUM_COVERAGE_LIMIT_FOR_NON_LOG_PICTURE} -e png
        
        Parameter_description:
        