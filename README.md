# *RamanomeSpec*
!!! Under re-construction.
----
## :eyes: Project Goal 
This project is focused on the development of AI-powered automated analytics for **single-cell Raman Spectroscopy (SCRS)** platform

----
## :exclamation: COPYING.txt
The file contains GNU general public license for  permissions of this strong copyleft license. All rights reserved. No commerical-relevant distributions are allowed.

----
## :innocent: Setup.py
The setup.py file is the build script for the package. The setup function from setuptools will build the package for upload to PyPI. Setuptools includes info about the package, version numbber, and which other packages are required for users.

----
## :baby: CheckVersion.py
To check the versions of all imported packages for the setup.py section

----
## :computer: Python Modules
----
### :old_key: :exclamation: Input Data Format
Data is separated by ``` '\t' ``` with shape of ```n samples x d features```. The 1st row represents wavenumbers. 

n samples should contain x background spectra in the end. or it can be merged from separate file.

![Data](https://github.com/ZJLEOWANG3/RamanomeSpec/blob/b229a0899022aed7831da27871425b9af7df5e91/media/data.format.png)

----
### :beginner: Prepare your SCRS Data
- **Load.py** : load all modules for other modules

- **Raman_read.py** : read txt files to pandas.dataframe for downstream analysis

- **Raman_preprocess.py** : Raman spectroscopy preprocessing module
    - [x] Background subtraction
    - [x] Baseline correction
    - [x] Smooth
    - **Optional**
        - [ ] Binning: to align your data based on differnet stepsize
        - [ ] Normalization: to normalize your data for later feature and AI relevant analysis

----
### :chart_with_upwards_trend: Get your peak & Assign its molecular
- **Raman_find_polymer.py** : find peaks and polymers of interests
    - :dart: **Targeted Analysis**
        - find single polymer
        - find multiple polymer
    - :kite: **Non-Targeted Analysis**
        - find all peaks using Raman_find_polymer.find_peak with a window of np.arange(1,30) using 
        [wavelet transform-based methodology](https://academic.oup.com/bioinformatics/article/22/17/2059/274284?login=true)
        
        - [ ] Under-development analysis sections

- **Raman_molecule.py** : side module for (wavenumber,molecule) pair
    - Load the molecular dictionary library
    - Given wavenumber, return molecular name

----
### :art: Advanced Analysis
- **Raman_cluster.py** : Intra-strain clustering modules for SCRS
    - Pick clustering algorithm using various metrics based on hard voting
    - Visualize strains w/ PCA biplot and related average, std, min, max spectra

- **Raman_chmap.py** : Clustering heatmap for SCRS to show differences
    - Clustering heatmaps for averaged Raman spectra given one baseline standard

- **Raman_stat.py** : Statistical modules for SCRS 
    - ANOSIM for comparison of inter-strain and intra-strain
    - Tukey test for ML accuracy
    - Jaccard and consistency computation
    - [ ] Underdevelopment
        - [ ] Chaos coefficients
        - [ ] **Network Analysis from another separate self-developed repository**
        - [ ] Correlation analysis for 2 trees

- **Raman_tree.py** : Tree-based analysis modles
    - Phenotyping tree establishment
    - 16S tree establishment w/ same label color
    - Mantel analysis for pairwise distance tree correlation

- [ ] **Raman_Ranking.py** : Under development; used to conduct feature ranking for SCRS
- [ ] **AI** : Under development; establish common FCNN, CNN, Transformer, Transformer variants for SCRS

----
### :books: Link to R markdown using specific eco-evolution R packages
- **Raman_BlombergK.Rmd**
    - To calculate Blomberg K statistics for the comparison of two trees (16S and phenotyping trees)
- **Raman_BLT.Rmd**
    - To perform branch length transformation using various methods, including delta, lambda, OU, rate change, two rate, and exponential rate transformations
- **Test_Tree.Rmd**
    - To test how the function of rescale of package geiger process the phylogenetic tree data by various branch length transformation methods

----
### :gift: Side Modules
- **Raman_plot.py** : SCRS data visualization
    - plot_line
    - add_line
    - [ ] Underdevelopment by incoporating other modules 
    
- **Raman_figmerge.py** : figure combination package
    - Save figures
    - Combine figures of intra-strain clustering based on different criterion 
        - one strain under different growth stages; 
        - one representative strain from each genus in stationary phase 1 with best accuracy; 
        - all strains under Bacillus genus in stationary phase 1 with best accuracy

- **Raman_color.py** : To align colors for above analysis
    - Generate customized colors
    - Generate label and color for microbial taxonomy by Raman analysis

----
### :crystal_ball: Examplary Workflow using the abovementioned modules

<!---
----
## DATA
- Processed : processed Raman datasets from raw datasets
    - RamanData_combined_TXT_ZIJIAN2
    - 
- molecule\_dict.json : current Raman polymer librai; dict key format : A ; B, where A is abbreviation, B is full name
- bacterial\_label.txt : bacteria label number ~ genus species
- average_accuracy.zijian.json : ML training accuracy of selected model
- 36-strain-changed2.nwk : 16S tree data for 36 strains
- TEMP : save temporary data

- OUTPUT : save important necessary datasets like the data in figures
    - row_colors.npy : numpy file for row_colors to visualize heatmap
    - genus_to_color.npy : numpy file for consistency of genus colors
    - HEATMAP : dir for heamap data
        - Multiple csv files with various filtration condition
    - Cluster : dir for cluster data
        - phase_taxonomy
            - data1 : PCA_cluster.csv
            - data2 : total_spec.csv
            - data3 : sep_spec.csv
            - verbose : verbose.txt
            - consistency of two methods : Consistency.txt
    - Merge2.npy : dictionary for script ClusterMerge2.py to customized figure merge; dict["phase_dict", "tax_dict"][phase/tax] = list of figure path
    - Dictionary : save the three key dictionary files for chmap, tree, clusters
        - chmap.npy : the original file for chmap
        - tree.npy : the original file for tree 
        - cluster.npy : the original file for clusters
    - Tree\_test : save the data of K-statiscs and p-values
        - Normed: the tree dataset is normalized
            - p_value.csv
            - K_value.csv
        - Unnormed : the tree dataset is not normalized
            - p_value.csv
            - k_value.csv
    - Stat : dir to save statistical results
        - Verbose_stat.csv
        - consistency describe: cons_des.csv
    - ANOSIM : dir to save ANOSIM input and output
        - df_nxd.npy : np_nxd+label
        - ANOSIM_verbose.txt : verbose for ANOSIM
    - EVR : dir to save the PCA EVR data
        - EVR.csv : plot the curve
    
----
## FIGURE
This section is to save figures. Figures are listed below.
- genus\_color.png : the genus color code for tree
- HEATMAP : dir to save clustering heatmaps 
- Cluster : dir to save Sep and Com
    - Sep : store sepearated figures 
        - phase_taxonomy
            - fig1 : PCA_cluster.png
            - fig2 : total_spec.png
            - fig3 : sep_spec.png
    - Com : store combined figures
        - phase_taxonomy
            - combined.png
    - Com2 : store combined figures2 (2 columns)
        - growth_stage dir
            - Exp.png, S1.png, S2.png, S3.png
        - Genus dir
            - Taxonomy,png
                       
- Tree : dir to save 16S & HC phenotyping tree
- Tree\_test : dir to save the test results
    - Branch length transform method names for dir
    - Statistics barplot for OU and trend : TreePair.png
- Stat : dir to save statistical results
        - Verbose_stat.png
- consistency.png: consistency of two clustering methods
----
## EXAMPLE
This section is to show some examples for the project results
- pipeline.py : standard pipeline for all the figures except for ML training
- ML\_model.py : to visualize the training accuracy results
- Tree.py : to establish and visualize a tree either for phenotyping unrooted data or 16S rooted data
- Cluster.py : to generate the clustering figure and data
- ClusterMerge.py : to merge figures of seperate cluster figure
- ClusterMerge2.py : to merge figures based on different conditions like Growth Stage or Genus type. It takes argv (1 - process and save figure path list; )
- TreePair.py : to visualize the results of R code for 16S-SCRS tree comparison
- Statistics.py : to perform some customized statistics
    - Verbose for the PCA-Kmeans Clustering results
    - Verbose for ANOSIM 
- consistency.py : to compute the consistency of two algorithms
- PCA_EVR.py : to compute the PCA components impacts on consistency or Jaccard
- pipeline.sh : bash script to perform my codes
)
--->
