# miRNA_predict
Project developed for the conclusion project for graduation at UFRGS - Computer Science.

Combining Performance and Diversity Measures for Optimizing Classification Ensembles via a Genetic Algorithm in the miRNA-Target Prediction Problem

MicroRNAs, also called miRNAs, are a large family of non-coding RNAs of approximately 22 nucleotides (nt) in length, which act as post-transcriptional gene silencers via translational repression or degradation of targets mRNAs, and have an important role in metabolism and genesis of different genetic diseases, such as cancers. The miRNA target prediction problem is considered a difficult challenge in the molecular biology area. There are millions of possible miRNA-mRNA possible combinations, and to experimentally find the functional combinations takes a large quantity of effort, therefore time and investment. 

The scientific community is actively researching computational approaches to overcome this cost with Machine Learning and their predictive models to better understand the interactions between miRNA-mRNA, and how they influence metabolic and disease processes. The purpose of this work is to study the effect of combining performance and diversity measures in a Genetic Algorithm's (GA) fitness function that learns the best combination of classifiers in an heterogeneous ensemble classifier in the miRNA-Target prediction problem.

Through experimentation, we've concluded that the challenge presented by the unbalanced and relatively small available datasets overshadows the possible benefits that the diversity measure could bring to the GA fitness function. Although the ensemble optimization combining performance and diversity measures has achieved better solutions than performance-based optimization in some cases, on average, the former solution did not surpass the latter. This doesn't allow us to conclude if the combination of performance and diversity measures results in better ensembles or not in our problem.
