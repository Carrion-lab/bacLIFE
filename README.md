# microLife: an automated genome mining tool for identification of lifestyle associated genes

microLife is a streamlined computational workflow that annotates bacterial genomes and performs large-scale comparative genomics to predict bacterial lifestyles and to pinpoint candidate genes, denominated  **lifestyle-associated genes (LAGs)**, and biosynthetic gene clusters associated with each lifestyle detected. This whole process is divided into different modules:

- **Clustering module**
	Predicts, clusters and annotates the genes of every input genome
- **Lifestyle prediction**
	Employs a machine learning model to forecast bacterial lifestyle or other specified metadata
- **Analitical module (Shiny app)**
	Results from the previous modules are embedded in a user-friendly interface for comprehensive and interactive comparative genomics. An example of the app with a demo dataset (genomes present in `data`) is available at http://178.128.251.24:3838/microLife_linux


![workflow](https://user-images.githubusercontent.com/69348873/231155358-7fbebb3c-f6f6-406a-989b-9d273b83aa1e.png)




