
# Project Biological data Analysis MCB COIMBRA 2022

## Goal

Compare different Machine learning approches for the prediction of the Enzyme Comission Number (EC) from Protein features such as aminoacid content.

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



## Creation of the training data:

To create a larger training and validation dataset we took inspiration in the training and validation data that had been used in the creation of the _ECPred_ a Java application for ECN prediction (see: https://github.com/cansyl/ECPred/ ). The benefits in using this dataset was that it is already curated for machine learning, is very vast. The drawback was that the repository only contains UniProt accession numbers and not the sequences of proteins. For this reason we had to download the sequences and build the dataset ourselfs afterwards.

To keep the scope of the project reasonlable we focused only on the prediction of the main enzyme identifier number EC{1,2,3,4,5}
"Hydrolases","Isomerases","Ligases","Oxidoreductases","Transferases"


### Accession number lists:

Files cointaining the identifiers where downloaded form github:

```{bash, eval= FALSE}
git clone --depth 1 --filter=blob:none --sparse https://github.com/cansyl/ECPred
cd ECPred
git sparse-checkout set ECPred\ Datasets
cd ./ECPred\ Datasets
rm -r -f !\(EC_Main*\)  # only keep main EC number classification
```

### Data retrival:


The script "createTrainingdata.R" includes the code to query the data from Uniprot from files containing Uniprot IDs taken from ECPret . 
It calls the helper functions sourced from utils : "getAASample.R" , "createAAmatrix.R"
In getAASample function a subset of identifiers gets sampled, then a costum query function "getUniProt_custom.R" is used, since the getUniProt function from protr will break if an identifier is bad, the custom version implements basic error handling using a tryCatch statement. The getAASample function also performs some basic data preprocession by removing gaps, and checking sequences for unknown AA identifiers.

### Creation of the evaluation dataset:

To evaluate the performance of the different ML models we used the human Porteom with assigned known EC numbers.
To create the evaluation data we queried the Uniprot Database for reviewed human enzyme sequences with know EC numbers, using their REST-Api.

```{r,eval=FALSE}
# create directories
for(i in c("uniprot","expasy")){
dir.create(path = paste0("data/raw/",i),recursive = TRUE)
}
# Download fasta files using wget:
for( i in 1:7){
  url<-paste0("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%20AND%20%28ec%3A",i,"%29%29%20AND%20%28reviewed%3Atrue%29")
  dest<-paste0("data/raw/uniprot/EC_",i,"_HomoSapiens")
  download.file(url,method = 'wget',destfile = dest)
  }
```

```{r,eval=FALSE}
EC_list<-empty.dump()
for(i in 1:7){
 EC_list[[i]]<-seqinr::read.fasta(
   paste0("data/raw/uniprot/EC_",i,"_HomoSapiens"),
                        seqtype = "AA",
                        seqonly = FALSE)
}
EC_matrix_list<-lapply(EC_list,SeqFastaAA_list_to_matrix)
# Add EC number:
names(EC_matrix_list)<-c(paste("EC",1:7))
## I call cbind.data.frame to keep the cells as int (cbind converts to character)
for(i in 1:7){
EC_matrix_list[[i]]<- cbind.data.frame( EC_matrix_list[[i]],"EC"=names(EC_matrix_list)[i])
}
# Bind rows to one big df:
comp_matrix<-do.call(rbind,EC_matrix_list)
saveRDS(comp_matrix,"data/comp_matrix_HSapiens.rds")
```

## Implementation of different Machine-learning approaches:



