## this is the final clustering system and the pipeline assembly
There are two phrases
1. Mapping and clustering the files 
2. Assemble the specific positions you would like to see

## In order to map and cluster fly genomes
1. srrSubmit.py with the specific args will run both srrSubmit.sh and srrSubmitM.sh depending on the input value of the files
```srrSubmit.py <experimentalNames.txt> ```
2. this should produce all the files summary and breakpoint after an hour or two 

## Assembly of the selected genomes
1. run pipelineAssembly - once per potential breakpoints
1.1. TODO: work on multiple potential breakpoints of a selected genome and not just one at a time 
2. this should produce the contigs files which will be blast into flybase


## OPTIONAL