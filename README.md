# DBgDel: Database-Enhanced Gene Deletion Framework for Growth-Coupled Production in Genome-Scale Metabolic Models

## About DBgDel

DBgDel is a powerful framework designed to determine gene deletion strategies. It integrates an extracted initial remaining gene pool from the database with advanced algorithms to narrow the algorithmic search space.

DBgDel comprises two steps: (1) **STEP 1**, extracting the gene set from a gene deletion database that corresponds to the essential core components for growth-coupled production for a given constraint-based model;
and (2) **STEP 2**, using an extended version of the gDel_minRN algorithm that incorporates the gene set extracted from STEP 1 as the initial gene pool to narrow the algorithmic search space.

DBgDel Framework Overview|
:-------------------------:|
| <img width="1000" alt="image" src="https://github.com/MetNetComp/DBgDel/blob/main/ov.png">

## Necessary Environments

To process the gene deletion strategy data downloaded from MetNetComp, you need an environment with the following:

- Python
- OS
- PANDAS

To run DBgDel, you need an environment with the following:

- MATLAB
- CPLEX
- COBRA Toolbox

## Download the Gene Deletion Strategy Data

The maximal gene deletion strategy data (specifying a model and a target metabolite) can be downloaded from MetNetComp via the following URL:

```
https://metnetcomp.github.io/database1/csv/<model>&<target metabolite>&core.csv
```
For example, to download the maximal gene deletion strategy of succ_e in the e_coli_core model, change the last part of the above link to ’/e_coli_core&succ_e.csv’.

The `Example_gene_deletion_data` file provides CSV files of 41 available maximal gene deletion strategies of the metabolites in e_coli_core model.

## Calculate the Initial Remaining Gene Set

`Calculate_Initial_Remaining_Gene_Set.py` provides a Python script to read all downloaded gene deletion strategy CSV files in a folder, calculate the initial remaining gene set, and save the result to a new CSV file.

Explanation:

1. Import Libraries:

- `os` for directory operations.
- `pandas` for reading and writing CSV files.

2. Function `read_csv_files(folder_path)`:

- Lists all CSV files in the specified folder.
- Reads each CSV file into a DataFrame.
- Converts the single column of data into a set and adds it to a list.

3. Function `calculate_intersection(data_sets)`:

- Takes a list of sets and returns their intersection.
- If the list is empty, returns an empty set.

4. Function `save_to_csv(data_set, output_file)`:

- Converts the intersection set to a DataFrame.
- Saves the DataFrame to a CSV file.

5. Function `main(folder_path, output_file)`:

- Calls the functions to read CSV files, calculate, and save the result.
- Prints a confirmation message with the path to the output file.

Execution:

- The script's main block sets the folder path containing the CSV files and the output file path.
- Calls the following to execute the entire process:

```
main()
```
- To use this script, the default folder and output path are set to the script's directory. You can also use the code in lines 60-64 to replace "path" and "path/Remaining_gene.csv" with the actual paths where your CSV files are located and where you want to save the result.

For more detailed information, please refer to the comments within the source code.

## Output Details

- The output CSV file `Remaining_gene.csv` will have a single column titled Remaining_gene, containing the genes in the initial remaining gene set.

## Running the Test Code

To run the test code for DBgDel, use the following command:

```
test()
```
The test() function performs the following steps:

1. Initializes the COBRA Toolbox environment using initCobraToolbox.
2. Loads:
  (1) A MATLAB matfile `e_coli_core.mat`, which contains a core metabolic model of E. coli.
  (2) A CSV file `e_coli_core_ex.csv`, which lists the initial remaining genes of e_coli_core.
Employs DBgDel to obtain the gene deletion strategy for growth coupling of succinate.

## Example Code
DBgDel can be used to calculate gene deletion strategies for various metabolites. Here are some examples:

- Example 1: Gene deletion strategy for pantothenate in iMM904.
- Example 2: Gene deletion strategy for succinate in iMM904.
- Example 3: Gene deletion strategy for biotin in iML1515.
- Example 4: Gene deletion strategy for riboflavin in iML1515.

The calculation results for these examples are available in the following files:
`Example_results/pantothenate_Strategy.mat`
`Example_results/succinate_Strategy.mat`
`Example_results/biotin_Strategy.mat`
`Example_results/riboflavin_Strategy.mat`

## Ablation Test and Baselines
DBgDel ablation variants are based on different initial remaining gene pool settings, including:

- Randomly Chosen (RC) genes.
- Growth Essential (GE) genes.
- $G_{\text{remain}}$ genes.

For example, the RC gene set is given by Function `RC_gene_selector(model, RC_size)`:

- Calls this function to randomly select a specified number of genes (`RS_size`) from a metabolic model (`model`).
- The resulting RC gene list is saved to `Ablation_Test_and_Baselines/RS_genes_model_name.csv`.

The other initial remaining gene pool settings are also available for testing in file path `Ablation_Test_and_Baselines`:
- RC genes: `Ablation_Test_and_Baselines/RC_genes_model_name.csv`.
- GE genes: `Ablation_Test_and_Baselines/GE_genes_model_name.csv`.
- $G_{\text{remain}}$ genes: `Ablation_Test_and_Baselines/G_remain_model_name.csv`.

Here are some examples using RC genes as the initial remaining gene pool:
- RC_Example 1: Gene deletion strategy for pantothenate in iMM904, using RC genes as the initial remaining gene pool.
- RC_Example 2: Gene deletion strategy for succinate in iMM904, using RC genes as the initial remaining gene pool.
- RC_Example 3: Gene deletion strategy for biotin in iML1515, using RC genes as the initial remaining gene pool.
- RC_Example 4: Gene deletion strategy for riboflavin in iML1515, using RC genes as the initial remaining gene pool.

The baseline methods include GDLS [1], optGene [2], gMCSE [3], and gDel_minRN [4]. 
GDLS and optGene are available in the MATLAB COBRA Toolbox [5]; please refer to [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html) for detailed instructions.

gMCSE and gDel_minRN are provided as open-source tools implemented in MATLAB; please refer to the original papers for detailed instructions and environments setup, and refer to [gMCSE](https://www2.mpi-magdeburg.mpg.de/projects/cna/etcdownloads.html) and [gDel_minRN](https://github.com/MetNetComp/gDel-minRN) for the resources.

## Output Details
The output contains a matrix `gvalue`, where:

- The first column lists the genes.
- The second column is a 0/1 vector indicating which genes should be deleted:
  - 0: Gene to be deleted.
  - 1: Gene to remain.
    
For more detailed information, please refer to the comments within the source code.

## References
[1] Lun D S, Rockwell G, Guido N J, et al. Large‐scale identification of genetic design strategies using local search[J]. molecular systems biology, 2009, 5(1): 296.

[2] Rocha I, Maia P, Rocha M, et al. OptGene: a framework for in silico metabolic engineering[C]//10th International Conference on Chemical and Biological Engineering. Portugal: University of Minho, 2008: 218-219.

[3] von Kamp A, Klamt S. Growth-coupled overproduction is feasible for almost all metabolites in five major production organisms[J]. Nature communications, 2017, 8(1): 15956.

[4] Tamura T, Muto-Fujita A, Tohsato Y, et al. Gene deletion algorithms for minimum reaction network design by mixed-integer linear programming for metabolite production in constraint-based models: gDel_minRN[J]. Journal of Computational Biology, 2023, 30(5): 553-568.

[5] Heirendt L, Arreckx S, Pfau T, et al. Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v. 3.0[J]. Nature protocols, 2019, 14(3): 639-702.
