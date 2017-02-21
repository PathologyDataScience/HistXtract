## HistXtract
HistXtract is a pipeline for extracting nuclear morphometry features from whole-slide images. This algorithm delineates individual cell nuclei in images of hematoxylin & eosin stained sections and calculates a quantitative feature profile describing the shape, texture and color of each nucleus. This information can be used with machine-learning algorithms to identify 

For more information, see our papers in [PLoS One](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081049) and [Laboratory Investigation](http://www.nature.com/labinvest/journal/v95/n4/abs/labinvest2014153a.html):

J Kong, LAD Cooper, F Wang, J Gao, G Teodoro, L Scarpace, T Mikkelsen, MJ Schniederjan, CS Moreno, JH Saltz, DJ Brat, *Machine-Based Morphologic Analysis of Glioblastoma Using Whole-Slide Pathology Images Uncovers Clinically Relevant Molecular Correlates*, PLoS One (2013) 8 (11), http://dx.doi.org/10.1371/journal.pone.0081049

LAD Cooper, J Kong, DA Gutman, WD Dunn, M Nalisnik, DJ Brat, *Novel genotype-phenotype associations in human cancers enabled by advanced molecular platforms and computational analysis of whole slide images*, Laboratory Investigation (2015) 95, 366–376; doi:10.1038/labinvest.2014.153; published online 19 January 2015

## Dependencies
HistXtract uses functions from [MatDigitalPathology](https://github.com/CancerDataScience/MatDigitalPathology) and the [OpenSlide](http://openslide.org/) project to analyze whole-slide images.
