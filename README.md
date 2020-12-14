[![DOI](https://zenodo.org/badge/71775920.svg)](https://zenodo.org/badge/latestdoi/71775920)

# KrATER (K-mer Analysis Tool Easy to Run)

I. Requirements

    1. Python libraries:
        Matplotlib
        Numpy
        Scipy
        
        RouToolPa (https://github.com/mahajrod/RouToolPa)
        
    2. Jellyfish 2

II. Installation

    0. Fast way
        - install jellyfish 2.x.x  from http://www.genome.umd.edu/jellyfish.html
        
        - sudo pip install routoolpa krater
        
        or
        
        - pip install --user routoolpa krater 

    1. Install requirements 
        - install jellyfish 2.x.x  from http://www.genome.umd.edu/jellyfish.html
    
        - install python libraries(not necessary if you going to install KrATER via pip)
            - run following command for global installation
        
                sudo pip install matplotlib numpy scipy 
            
            - or following command for local installation if you don't have root permissions
        
                pip install --user matplotlib numpy scipy
        
        - install RouToolPa
            From pip with root permission:
                sudo pip install routoolpa
            
            From pip without root permissions:
                pip install --user routoolpa
        
            Get RouToolPa
                git clone https://github.com/mahajrod/routoolpa
        
            Add following strings to ~/.profile and ~/.bashrc (create files if absent). Don't forget to replace <ROUTOOLPA_DIR> with actual path
            
                PYTHONPATH=${PYTHONPATH}:<ROUTOOLPA_DIR>
                PATH=${PATH}:<ROUTOOLPA_DIR>
                export PYTHONPATH
                export PATH
                
            Run in terminal
                source ~/.profile
        
        
    2. Install KrATER
        Use variant 3 to install newest version of KrATER
        Variant 1: install using pip (old version)
        
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
        
