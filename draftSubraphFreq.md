# Predictions of selection against instability: does subgraph frequency relate to quasi sign-stability in food webs
Jonathan Borrelli  
Stony Brook University  
Stony Brook NY, 13850  

## Abstract
Food webs structure can be defined by the particular frequencies of three node subgraphs. Of the thirteen possible configurations of three species food webs, some are consistently over and underrepresented in larger (whole community) food webs. This is a robust pattern that spans multiple environments. Any potential explanation must also be able to apply without respect to the particulars of marine, freshwater, or terrestrial environments. I argue that the elimination of unstable subgraphs during the development of the food web can explain the observed pattern. A clear prediction of this hypothesis is that there should be differences in the probability of stability in different subgraphs, and that this probability should be related to their frequency in food webs. Using 50 food webs collected from a variety of databases I determined the frequency of each of the thirteen possible subgraphs with respect to randomized webs. Then by numerical simulation I determined the quasi sign stability (QSS) of each subgraph. My results clearly show the relationship between QSS and over/underrepresentation of the different subgraphs. **(172 words)**

## Introduction  

Larger networks are composed of many smaller subnetworks that are assembled together. Any large network with N nodes can be decomposed into smaller networks of size 1 to N-1. It is likely that the dynamics that occur within larger networks will at least be influenced by the dynamics of their respective components.  

There are thirteen possible (connected, directed) configurations of three nodes, five of which require only single direction links and 8 which combine single and bi-directional links. **Milo et al. (2002)** showed that the frequency that these thirteen configurations of three nodes are observed change in different types of networks (e.g., food webs, mutualistic networks, gene transcription networks, neural networks, etc.). For example, in food webs the tritrophic chain tends to be overrepresented compared to random. Those subgraphs which tend to be overrepresented are commonly termed _motifs_.

Further analysis of the composition of food webs with respect to the thirteen possible three node subgraphs shows that there tends to be distinct patterns of representation. Representation refers to the frequency with which a given subgraph is observed in a network with respect to a null model of randomized interactions. **Camacho et al. (2002, 2007)** determined that these patterns are observed regardless of the type of environment the food web is in (e.g., freshwater, marine, terrestrial), and that the best food web models were able to reproduce these patterns. 


 

## Methods  
### Data  
***\*NOTE TO SELF: ADD MORE INFO ABOUT EACH OF THE DATASETS*\***

I used 50 food webs collected from a variety of sources. Three food webs were downloaded from the [Dryad Digital Repository](http://datadryad.org/resource/doi:10.5061/dryad.c213h) (Roopnarine & Hertog 2012a, 2012b). Another seven were available from [Ecological Archives](http://esapubs.org/Archive/search.php) (Hechinger et al. 2011; Mouritsen et al. 2011; Thieltges et al. 2011; Zander et al. 2011; Preston et al. 2012). Fourteen webs were provided by Jennifer Dunne of the [PEaCE Lab](http://peacelab.net/) (Baird & Ulanowicz 1989; Warren 1989; Polis 1991; Hall & Raffaelli 1991; Martinez 1991; Christensen & Pauly 1992; Havens 1992; Goldwasser & Roughgarden 1993; Opitz 1996; Waide & Reagan 1996; Yodzis 1998, 2000; Christian & Luczkovich 1999; Martinez et al. 1999; Memmott et al. 2000; Link 2002) that were analyzed in (Dunne et al. 2002, 2004). The remaining 26 food webs were downloaded from the [Interaction Web Database](http://www.nceas.ucsb.edu/interactionweb/html/thomps_towns.html) (Jaarsma et al. 1998; Townsend et al. 1998; Thompson & Townsend 1999, 2000, 2003, 2005; Thompson & Edwards 2001).

### Subgraph Frequency
The frequency of each subgraph was found using the `triad.census` function from the `igraph` package (**Csardi and Nepusz 2006**) in R _version 3.0_ (**R Core Team 2014**). Counts of each of the thirteen subgraphs were determined for each of the 50 food webs described above. Each frequency was then compared against a null distribution. 

To create the null distribution each of the fifty adjacency matrices (each food web) were permuted, maintaining both the number of predators of a species and the number of prey a species has (maintiaining row and column sums). 1000 permuted matrices were generated for each food web, and the frequency of each subgraph type was found. Z-scores were computed using the formula:  

$$
z_i = \frac{(X_i - \overline{X_i})}{\sigma^2}
$$  

where $$$X_i$$$ is the frequency of the $$$i^{th}$$$ subgraph in each empirical food web, $$$\overline{X_i}$$$ is the mean  frequency of the $$$i^{th}$$$ subgraph in the permuted matrices, and $$$\sigma_i^2$$$ is the variance. 

### Subgraph Stability
Using numerical simulations in R _version 3.0_ (**R Core Team 2014**) I determined the probability of a given subgraph will be stable (quasi sign-stability, QSS).  

## Results  
I found that quasi sign-stability of a given subgraph is positively correlated with the frequency it is observed in empirical networks (relative to random). 

![image](\users\jonathanborrelli\Desktop\Github\Subgraph-Stability\Plots\stability-freqPLOT.png)

Subgraphs with a higher quasi sign-stability occur more frequently than expected by random chance, while those with lower quasi sign-stability tend to occur less frequently than expected by chance. 

## Discussion  

## References  
  Csardi G, Nepusz T. 2006. The igraph software package for complex network research, InterJournal,
  Complex Systems 1695. [http://igraph.org](http://igraph.org)
