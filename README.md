# Core-lncLocation
A computational tool based multi-source heterogeneous features to predict lncRNA subcellular localization.

1. LncLocation.py

	Usage: 	python LncLocation.py input_file output_profix
    
    		Input_file must be fasta format
    
    		Users can use the output_profix key value to specify the output folder of the predicted results file

	Note:	1.You must install the necessary programs to support LncLocation, including the R programming language and python.

		2.Your Python version number must be 3 or more.

		3.The input_file must be in fasta format and you can specify the input path to the file, for example "./test.fa".
		
		4.The output__profix is the result folder and you can name it. 

		5.The test.fa is an example, and you can use the following instructions to run the test program: "python LncLocation.py ./test.fa result".

2. Requirements

	The environment you need includes:        R (version >= 4.0.2)       python (version >= 3.7)

			R packages:    LncFinder (version >= 1.1.4)    seqinr (version>=	3.6-1）

			python packages:    sklearn (version >= 0.22.1)    pandas (version >= 0.25.3)    numpy (version >= 1.18.1)    
				               rpy2 (version >= 3.3.5)    Keras (version >= 2.4.3)    tensorflow (version >= 2.3.0)

3. Installing

	You can configure it from the Linux command line using the following statement.
	
	1) Installation of R language and R language package：
		
		R language:
		
			sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common
			
			sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
      
			sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
			
			sudo apt install r-base
			
			sudo apt install build-essential
			
		R packages:
		
			R
			
			install.packages("LncFinder")
      
			install.packages("seqinr")
		
	2) Installation of python packages:
		
			pip install sklearn==0.22.1

			pip install pandas==0.25.3

			pip install numpy==1.18.1

			pip install keras==2.4.3

			pip install tensorflow==2.3.0

			sudo apt-get install python-rpy2

	3) Installation of viennaRNA package:
		
			wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.14.tar.gz

			tar -zxvf ViennaRNA-2.4.14.tar.gz

			cd ViennaRNA-2.4.14

			./configure

			make

			sudo make install

4. Sample output

        NR_002599.1|@2@|Homo|sapiens|small|nucleolar|RNA|host|gene|6|(SNHG6),|long|non-coding|RNA:
          cytoplasm:0.9601980499899129	exosome:0.00909585233825203	nucleus:0.006730082757155168	ribosome:0.023976014914679853
          predicted results:cytoplasm
        NR_002212.4|@2@|Homo|sapiens|nudix|hydrolase|4|pseudogene|1|(NUDT4P1),|non-coding|RNA:
          cytoplasm:0.00017943483633488412	exosome:0.04874539208880222	nucleus:4.6539513170652725e-06	ribosome:0.9510705191235461
          predicted results:ribosome
        NR_024037.1|@1@|Homo|sapiens|rhabdomyosarcoma|2|associated|transcript|(non-protein|coding)|(RMST),|long|non-coding|RNA:
          cytoplasm:0.0034342366478661065	exosome:0.07221840896907318	nucleus:0.5519930647840583	ribosome:0.372354289599003
          predicted results:nucleus
        NR_001549.1|@4@|Homo|sapiens|testis-specific|transcript,|Y-linked|19|(non-protein|coding)|(TTTY19),|long|non-coding|RNA:
          cytoplasm:2.9193545634207932e-06	exosome:0.9992812304227594	nucleus:2.404225867754603e-05	ribosome:0.0006918079639996525
          predicted results:exosome


please check reference:
[1] Feng Shiyao1,2, Liang Yanchun1,2, Du Wei1, Zhang Yu1, Li Ying1* [2020]
 LncLocation: efficient subcellular location prediction of long non-coding RNAs based multi-source heterogeneous features fusion.
	
--



