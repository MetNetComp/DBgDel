# DBgDel: Database-Enhanced Gene Deletion Framework for Growth-Coupled Production in Genome-Scale Metabolic Models

## About DBgDel

DBgDel is a powerful framework designed to determine gene deletion strategies. It integrates an extracted initial remaining gene pool from the database with advanced algorithms to narrow the algorithmic search space.

DBgDel comprises two steps: (1) **STEP 1**, extracting the gene set from a gene deletion database that corresponds to the essential core components for growth-coupled production for a given constraint-based model;
and (2) **STEP 2**, using an extended version of the gDel_minRN algorithm that incorporates the gene set extracted from STEP 1 as the initial gene pool to narrow the algorithmic search space.

DBgDel Framework Overview|
:-------------------------:|
| <img width="1000" alt="image" src="https://github.com/MetNetComp/DBgDel/blob/main/ov.png">

## Necessary Environments

To run DBgDel, you need an environment with the following:

- MATLAB
- CPLEX
- COBRA Toolbox

## Running the Test Code

To run the test code for DBgDel, use the following command:

```matlab
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

- Example 1: Gene deletion strategy for biotin in iML1515.
- Example 2: Gene deletion strategy for riboflavin in iML1515.
- Example 3: Gene deletion strategy for pantothenate in iMM904.
- Example 4: Gene deletion strategy for succinate in iMM904.
The calculation results for these examples are available in the following files:
`Example_results/biotin_Strategy.mat`
`Example_results/riboflavin_Strategy.mat`
`Example_results/pantothenate_Strategy.mat`
`Example_results/succinate_Strategy.mat`

## Output Details
The output contains a matrix `gvalue`, where:

- The first column lists the genes.
- The second column is a 0/1 vector indicating which genes should be deleted:
  - 0: Gene to be deleted.
  - 1: Gene to remain.
    
For more detailed information, please refer to the comments within the source code.
