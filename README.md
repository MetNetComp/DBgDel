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

The maximal gene deletion strategy data (specifying a model and a target metabolite) can be downloaded from [MetNetComp](https://metnetcomp.github.io/database1/indexFiles/index.html) via the following URL:

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

## Output Details
The output contains a matrix `gvalue`, where:

- The first column lists the genes.
- The second column is a 0/1 vector indicating which genes should be deleted:
  - 0: Gene to be deleted.
  - 1: Gene to remain.
    
For more detailed information, please refer to the comments within the source code.

## DBgDel Ablation Test
DBgDel ablation variants are based on different initial remaining gene pool settings, including:

- Randomly Chosen (RC) genes.
- Growth Essential (GE) genes.
- $G_{\text{remain}}$ genes.

For example, the RC gene set is given by Function `RC_gene_selector(model, RC_size)`:

- Calls this function to randomly select a specified number of genes (`RC_size`) from a metabolic model (`model`).
- The resulting RC gene list is saved to `initial_remaining_gene/RC_genes_model_name.csv`.

The other initial remaining gene pool settings are also available for testing in file path `initial_remaining_gene`:
- RC genes: `initial_remaining_gene/RC_genes_model_name.csv`.
- GE genes: `initial_remaining_gene/GE_genes_model_name.csv`.
- $G_{\text{remain}}$ genes: `initial_remaining_gene/G_remain_model_name.csv`.

Here are some examples using RC genes as the initial remaining gene pool:
- RC_Example 1: Gene deletion strategy for pantothenate in iMM904, using RC genes as the initial remaining gene pool.
- RC_Example 2: Gene deletion strategy for succinate in iMM904, using RC genes as the initial remaining gene pool.
- RC_Example 3: Gene deletion strategy for biotin in iML1515, using RC genes as the initial remaining gene pool.
- RC_Example 4: Gene deletion strategy for riboflavin in iML1515, using RC genes as the initial remaining gene pool.

The output file format of the above ablation variants remains the same as the original DBgDel.

## Baseline Test (Updated)

Baseline methods for comparison include **gMCSE** [1], **gDel_minRN** [2], **GDLS** [3], and **optGene** [4]. Among them, gMCSE is a minimal cut set (MCS)-based method, while the others are based on elementary flux vectors (EFVs).
Below, we provide detailed experimental setups and source code to ensure fair reproduction and facilitate further research.
All the codes are available in file path `Baselines`.

### 1. gMCSE
gMCSE is provided as an open-source tool implemented in MATLAB; please refer to the original papers for detailed instructions, and refer to [gMCSE](https://www2.mpi-magdeburg.mpg.de/projects/cna/etcdownloads.html) for the resources.

**(a) Environments Setup**

To implement large-scale baseline experiments on gMCSE, we used the API function `CNAgeneMCSEnumerator2` from CellNetAnalyzer [5] (ver. 2023.1). For more details and access to the latest version of CellNetAnalyzer, please refer to [CellNetAnalyzer](https://www2.mpi-magdeburg.mpg.de/projects/cna/etcdownloads.html).
For small-scale tests on a few targets, we recommend using the GUI provided by CellNetAnalyzer, as it enables more straightforward parameter configuration.

gMCSE requires MATLAB computational tools such as EFMTool, which rely on Java to handle large-scale calculations. In our cases, MATLAB runs on a Java Virtual Machine (JVM), and the default memory limit may be insufficient for memory-intensive tasks like MCS computation. This can result in errors, slowdowns, interruptions, or incorrect termination of computations. To prevent this, we provide the function `setJavaHeapSize(sizeGB),` which is strongly recommended to use before experiments for manually adjusting the Java heap size based on the available system RAM.


**(b) gMCSE Parameter Setup**

The gMCSE method requires several key parameters (`D`, `d`, `T`, and `t`) to set up the Desired Region and Target Region during calculations. Incorrect setup of these parameters can lead to false results.The recommended setup for gMCSE is as follows:
- Desired Region (D * r <= d): Ensure growth is greater than or equal to the GR_threshold. (The metabolism must support growth.)
- Target Region (T * r <= t): Ensure growth is greater than or equal to the GR_threshold, and production is less than or equal to the PR_threshold. (Metabolism must no longer be able to grow while no product is produced.)

For compatibility with gMCSE, it is necessary to rewrite these constraints as matrix-vector multiplications:
- D * r <= d (1xN matrix and 1x1 vector): -growth <= -GR_threshold
- T * r <= t (2xN matrix and 2x1 vector): -growth <= -GR_threshold, production <= PR_threshold

To simplify the setup of gMCSE, we provide two functions that help easily define the constraints as described above:
- `initializeDesiredRegion(cnap_numr, growth_idx, GR_threshold)`: Initializes `D` and `d` based on the given growth reaction in the model and the GR_threshold.
- `initializeTargetRegion(cnap_numr, growth_idx, production_idx, GR_threshold, PR_threshold)`: Initializes `T` and `t` based on the given growth reaction, target production reaction in the model, GR_threshold, and PR_threshold.

**(c) Result Reports on Computational Experiments**

We provide scripts to reproduce the tests on the reported results for gMCSE:
- `Baselines/test_gMCSE_e_coli_core.m`: computational experiments on **e_coli_core** model (refer to **TABLE VIII**).
- `Baselines/test_gMCSE_iMM904.m`: computational experiments on **iMM904** model (refer to **TABLE X**).
- `Baselines/test_gMCSE_iML1515.m`: computational experiments on **iML1515** model (refer to **TABLE XII**).

Note_1: As mentioned above in **(a) Environment Setup**, we provide an additional function `cleanUpBeforeNewModel()` to reset MATLAB memory and parallel pool usage before loading and running a new model. 
We strongly recommend testing one model at a time and calling this function before each new round of computation.

Note_2: The final report is based on the results of the `GRPRchecker()` function, which checks whether the resulting gene deletions achieve GCP for each target in the original model considering the GPR rules setup in this study.

The results for each target are recorded, and the experiment summary is reported in the file `Baselines/model_name_results_check.csv`
- ProductionIdx: Refers to the target ID in each model.
- TimeCost: Elapsed time for each target, along with the average elapsed time.
- Status: The output from gMCSE, with the following codes (only status codes 0 and 3 give feasible gene deletions):
  - 0: Successful
  - 1: Timeout with no solution
  - 2: Infeasible 
  - 3: Timeout with some solutions
- DeletedGenes: The set of deleted genes.
- GR: The GR check result based on the knockout of DeletedGenes.
- PR: The PR check result based on the knockout of DeletedGenes.
- Success:
  - 0: Failed GCP
  - 1: Successful GCP

**(d) Example Test**

We provide a script to make a quick small-scale example test on gMCSE:
- `Baselines/test_gMCSE.m`: test on 10 target metabolites on **iMM904** model.

### 2. gDel_minRN

**(a) Computational Experiments Setup**

gDel_minRN is provided as an open-source tool implemented in MATLAB; please refer to the original papers for detailed instructions and environment setup, and refer to [gDel_minRN](https://github.com/MetNetComp/gDel-minRN) for the resource download.

**(b) Result Reports on Computational Experiments**

We provide the script to reproduce the tests on the reported results for gDel_minRN:
- `Baselines/test_gDel_minRN.m`: computational experiments on all three models (**e_coli_core, iMM904, and iML1515**).


### 3. GDLS and optGene

**(a) Computational Experiments Setup**

GDLS and optGene are available in the MATLAB COBRA Toolbox [6]; please refer to [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html) for detailed instructions and resource download.

**(b) Result Reports on Computational Experiments**

- `Baselines/test_GDLS_optGene.m`: computational experiments on all three models (**e_coli_core, iMM904, and iML1515**).

## References
[1] von Kamp A, Klamt S. Growth-coupled overproduction is feasible for almost all metabolites in five major production organisms[J]. Nature communications, 2017, 8(1): 15956.

[2] Tamura T, Muto-Fujita A, Tohsato Y, et al. Gene deletion algorithms for minimum reaction network design by mixed-integer linear programming for metabolite production in constraint-based models: gDel_minRN[J]. Journal of Computational Biology, 2023, 30(5): 553-568.

[3] Lun D S, Rockwell G, Guido N J, et al. Large‐scale identification of genetic design strategies using local search[J]. molecular systems biology, 2009, 5(1): 296.

[4] Rocha I, Maia P, Rocha M, et al. OptGene: a framework for in silico metabolic engineering[C]//10th International Conference on Chemical and Biological Engineering. Portugal: University of Minho, 2008: 218-219.

[5] Klamt S, Saez-Rodriguez J, Gilles E D. Structural and functional analysis of cellular networks with CellNetAnalyzer[J]. BMC systems biology, 2007, 1: 1-13.

[6] Heirendt L, Arreckx S, Pfau T, et al. Creation and analysis of biochemical constraint-based models using the COBRA Toolbox v. 3.0[J]. Nature protocols, 2019, 14(3): 639-702.
