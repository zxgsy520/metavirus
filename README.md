# metavirus
Metavirus analysis process

Dependent software
---------------------
### find_virus.py (Prediction of viral and phage sequences from metagenomic)
* [VirSorter2](https://github.com/jiarong/VirSorter2)
* [DeepVirFinder](https://github.com/jessieren/DeepVirFinder)
* [VIBRANT](https://github.com/AnantharamanLab/VIBRANT)
* [VirFinder](https://github.com/jessieren/VirFinder)
* [VirKraken](https://github.com/Strong-Lab/VirKraken)
* [checkv](https://jgi.doe.gov/data-and-tools/software-tools/checkv/)
* [Unicycler](https://github.com/rrwick/Unicycler) 修改Unicycler中的spades_func.py脚本（使用宏病毒模型）

#### Taxonomic assignment
* [vConTACT2](https://bitbucket.org/MAVERICLab/vcontact2/wiki/Home) 安装：conda create --prefix=/software/meta/vcontact2/v0.9.19 -c bioconda vcontact2=0.9.19  mcl blast diamond numpy=1.19.5 pandas=0.25.3
* [CAT](https://github.com/dutilh/CAT)

### Required R packages
* [complexheatmap](https://bioconductor.org/packages/3.14/bioc/html/ComplexHeatmap.html)
* 
