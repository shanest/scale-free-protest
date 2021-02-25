# scale-free-protest

Companion repository to the paper ``How Social Networks Affect the Repression-Dissent Puzzle''.

Runs network simulations of the spread of protest with varying network structure and repression mechanisms.

Data for some empirical results can also be found at: http://dx.doi.org/10.17605/OSF.IO/6U9RW

Setup:
* Make a new venv/conda environment, with python >= 3.6
* pip install -r requirements.txt

Basic usage:

```
python simulations.py --exp {file}.yml --num_procs 1 --out_dir OUTPUT_DIR
```

The experiments from the paper can be found as specified in YAML files in the `exps` directory.


## Analysis scripts

The simulations will output a large amount of raw data about the dynamics of the protest across trials.  There are several scripts in `analysis` for analyzing this data, briefly described here:

1. `01_ProcessSimulationData_v1.R`
	- Merges datasets
	- Input:
		- exp3a-graph_type=powerlaw_cluster_graph-num_nodes=1000-threshold_type=uniform-m=3-repression_type=node_removal.csv
		- exp2a-threshold_type=uniform-repression_type=node_removal-graph_type=watts_strogatz_graph-k=15-num_nodes=1000.csv
		- exp1a-1-num_nodes=1000-graph_type=scale_free_graph-threshold_type=uniform-repression_type=node_removal.csv
		- exp1a-2-num_nodes=1000-graph_type=scale_free_graph-repression_type=node_removal-threshold_type=uniform.csv
		- exp1a-3-repression_type=node_removal-num_nodes=1000-graph_type=scale_free_graph-threshold_type=uniform.csv
		- exp1a-4-repression_type=node_removal-num_nodes=1000-graph_type=scale_free_graph-threshold_type=uniform.csv
	- Output:
		- processedData/01_experimentsCombined.csv

2. `02_MakeFigures_v1.R`
	- Make plots to explore the data.  Main goal is to facet results by network type, DV, one facet plot per variable of interest.
	- Input:
		- processedData/01_experimentsCombined.csv
	- Output:
		- A ton of figures

3. `03_RandomForestRegression_v1.R`

4. `04_Simulations_Regression_InvestigatedScaleFreeAndRepression.R`
	- The purpose of this script is to see if regression inferences change based on regression on different subsets of the scale-free dataset.  I have found in its parent script that repression_rate can sometimes have heterogenous effects, when low enough and sample is small enough.
	- Input:
		- processedData/01_experimentsCombined_eigenCore.csv
	- Output:
		- Lots of figures

5. `05_Simulations_RepressionFigures_v1.R`
	- The purpose of this script is to make the figures detailed on Page 132 of my Moleskine journal.
	- Input:
		- processedData/01_experimentsCombined_eigenCore.csv
	- Output:
		- Lots of faceted figures

6. `06_tStatSampleSizeChange_v1.R`
	- The purpose of this script is to see how a t-statistic changes as increase sample size.

7. `07_FatalityRate_v1.R`
	- The purpose of this script is to generate an estimate of the fatality rate at protests.  It uses ACLED and NAVCO.
	- Input (in OSF project linked above):
		- ACLED_protests_1900-01-01-2019-08-07.csv
		- navco3.csv
	- Output:
		- Just an average value, no saved output.

8. `08_RegressionResultsOverTime_v1.R`
	- The purpose of this script is to see if regression results change as we get new data about an event.  The idea is that it replicates how our knowledge changes as new data comes in with history.  
	- Input:
		- Data/MassMobilization_Binghampton/mm_public_120615.csv
		- Data/MassMobilizationsAutocracyDataset/EventLevel/events.csv
	- Output: