
# Project 

## Goal

Train a Machine learning Algorithm to predict Enzyme Comission Number (EC -number ) from Protein features.

 - Sequence 	> Structure 		> Function
 - Sequence 	> (features -> ML)      > Function

## EC - Numbers
| Class                | Reaction catalyzed                                                                                                        | Typical reaction                                | Enzyme example                          |
|:-------------|:-------------------------------|:-------------|:-------------|
| EC 1 Oxidoreductases | Oxidation/reduction reactions; transfer of H and O atoms or electrons from one substance to another                       | AH + B → A + BH (reduced) A + O → AO (oxidized) | Dehydrogenase, oxidase                  |
| EC 2 Transferases    | Transfer of a functional groupf from one substance to another. The group may be methyl-, acyl-, amino- or phosphate group | AB + C → A + BC                                 | Transaminase, kinase                    |
| EC 3 Hydrolases      | Formation of two products from a substrate by hydrolysis                                                                  | AB + H2O → AOH + BH                             | Lipase, amylase, peptidase, phosphatase |
| EC 4 Lyases          | Non-hydrolytic addition or removal of groups from substrates. C-C, C-N, C-O or C-S bonds may be cleaved                   | RCOCOOH → RCOH + CO2 or [X-A+B-Y] → [A=B + X-Y] | Decarboxylase                           |
| EC 5 Isomerases      | Intramolecule rearrangement, i.e. isomerization changes within a single molecule                                          | ABC → BCA                                       | Isomerase, mutase                       |
| EC 6 Ligases         | Join together two molecules by synthesis of new C-O, C-S, C-N or C-C bonds with simultaneous breakdown of ATP             | X + Y + ATP → XY + ADP + Pi                     | Synthetase                              |
| EC 7 Translocases    | Catalyse the movement of ions or molecules across membranes or their separation within membranes                          |                                                 | Transporter                             |
(from <https://www.abbexa.com/enzyme-commission-number> *altered*)

## Challenges

-   SAR -paradoxon: SAR: \_Structure--Activity *Relationship* ... refers to the fact that it is not the case that all similar molecules have similar activities.



## Training Data:


### Creation of training data:

To create a larger training and validation dataset we used the training and validation data that had been used in the creation of _ECPred_ a Java application for ECN prediction (see: https://github.com/cansyl/ECPred/ ). The benefits in using this dataset was that it is already curated for machine learning, is very vast. The drawback was that the repository only contains UniProt accession numbers and not the sequences of proteins. For this reason we had to download the sequences and build the dataset ourselfs afterwards.

To keep the scope of the project reasonlable we focused only on the prediction of the main enzyme identifier number EC{1,2,3,4,5}
"Hydrolases","Isomerases","Ligases","Oxidoreductases","Transferases"


### To get accession number lists:


```{bash, eval= FALSE}
git clone --depth 1 --filter=blob:none --sparse https://github.com/cansyl/ECPred
cd ECPred
git sparse-checkout set ECPred\ Datasets
cd ./ECPred\ Datasets
rm -r -f !\(EC_Main*\)  # only keep main EC number classification
```

## The ECPred 

```{bash}
cd ./ECPred/ECPred\ Datasets
tree
```

The script "createTrainingdata.R" includes the code to query the data from Uniprot. 

## Creation of the evaluation dataset:

TTo evaluate the performance of the different ML models we used the human Porteom 
To create the evaluation data we queried the Uniprot Database for reviewed human enzyme sequences with know EC numbers, using their REST-Api.



