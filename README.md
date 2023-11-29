# Chemometric data analysis tutorials

This repository contains a series of tutorials on multivariate analysis of metabolic profiling datasets: 
 - Import, scaling & normalisation.ipynb: Introduction to normalisation, scaling and data transformations.  
 - Multivariate Analysis - PCA.ipynb: Multivariate chemometric analysis using Principal Component Analysis.
 - Multivariate Analysis - Supervised Analysis with PLS-DA.ipynb: Discrimnination of 2 classes with PLS-DA.
 
To run these tutorials download the contents of this repository and run the Jupyter Notebooks. Alternatively, these can be run on the browser via Google Colaboratory by clickinking on the links at the IPTC - Data Analysis Course [repository](https://github.com/IPTC-DataAnalysisCourse).  

All the data required to run these notebooks is provided on the 'data' folder. The dataset used in this tutorial comes from the following publication:
- Lovestone, S., Francis, P., Kloszewska, I., Mecocci, P., Simmons, A., Soininen, H., Spenger, C., Tsolaki, M., Vellas, B., Wahlund, L.-O., Ward, M. and (2009), AddNeuroMed—The European Collaboration for the Discovery of Novel Biomarkers for Alzheimer's Disease. Annals of the New York Academy of Sciences, 1180: 36-46. [https://doi.org/10.1111/j.1749-6632.2009.05064.x](https://doi.org/10.1111/j.1749-6632.2009.05064.x).

It is a set of urine samples of 577 individuals processes by using Liquid Chromatography - Mass Spectrometery instrument.  
There are two main biological sources of variation in the dataset:
- Gender/Sex (0: Biological female (n=294)) , 1: Biological male (n=283))
- Age 


## Installation instructions:
To run these tutorials localy, first download and install the latest Anaconda Python 3.10.x distribution available. 
The following packages also need to be installed (installed specific version if stated for compatibility):
- plotly: Available via conda (open an Anaconda prompt and type "conda install plotly")
- SciPy: 1.11.3
- NumPy: 1.26.0
- scikit-learn: 1.2.1
- Pandas: 2.1.1
- Matplotlib: 3.6.2
- seaborn
- kneed: 0.8.5 

All other dependencies required to run the tutorials (pyChemometrics toolbox) are bundled with the repository, and no specific installation is required. 

After installing plotly, either clone this repository or download as .zip file. Access the notebooks contents via the Jupyter notebook environment. 

For more information on how to use Jupyter notebooks see [Using the Jupyter notebook](https://docs.anaconda.com/ae-notebooks/user-guide/basic-tasks/apps/jupyter/)

