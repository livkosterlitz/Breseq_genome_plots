These programs can be used to generate a genome-wide plot for multiple samples [example plot](https://github.com/livkosterlitz/Breseq_genome_plots/blob/main/BasicRun/Breseq_plots/out.pdf). Step 1: running the `breseq_HTML_to_CSV.py` program to convert the _breseq_ output into CSV files. Step 2: running the `breseq_plotting.R` program to create the genome-wide plot. The example provided in the repository is from [Jordt. et. al. 2020](https://doi.org/10.1038/s41559-020-1170-1). 


## Step 1 `breseq_HTML_to_CSV.py`

#### Inputs
The `breseq_HTML_to_CSV.py` program has one required input and two optional arguments.

Required:
* -f/--file_pairs: User-specified sample name and paths to be processed (e.g., sample sample.html). User can pass multiple pairs (e.g., sample1 sample1.html sample2 sample2.html)

Optional:
* -p/--plotting_csv: If passed the program will create a combined CSV file for all of the samples processed and a CSV file in the correct format for a downstream R program for genome plotting. 
* -c/--treatment_csv: A CSV with two columns, treatment and strains (i.e., samples), which will set groupings of samples for the downstream R program for genome plotting [example CSV output](https://github.com/livkosterlitz/Breseq_genome_plots/blob/main/BasicRun/Breseq_CSV_output/treatments.csv)  

#### Command to run example
```bash
python /Users/oliviakosterlitz/Programs/breseq_HTML_to_CSV.py -p T -c treatments.csv -f EC_29_anc EC_29_anc/index.html EC_29_1 EC_29_1/index.html EC_29_2 EC_29_2/index.html EC_29_3 EC_29_3/index.html EC_29_4 EC_29_4/index.html EC_29_5 EC_29_5/index.html EC_29_6 EC_29_6/index.html
```

#### Outputs
The program puts file output in the working directory. 

Ouput CSV files [example CSV output](https://github.com/livkosterlitz/Breseq_genome_plots/blob/main/BasicRun/Breseq_CSV_output/EC_29_1_predicted_mutations.csv)
> For each sample, the program converts the _breseq_ index.html into CSV files. One for each section: predicted mutations, unassigned missing coverage evidence. and unassigned new junction evidence. 

Additional outputs if plotting option turned on:
1. combined_predicted_mutations.csv [example CSV output](https://github.com/livkosterlitz/Breseq_genome_plots/blob/main/BasicRun/Breseq_CSV_output/combined_predicted_mutations.csv)
> A CSV file with all of the predicted mutations for each sample. 
2. plot_variants.csv [example CSV output](https://github.com/livkosterlitz/Breseq_genome_plots/blob/main/BasicRun/Breseq_CSV_output/plot_variants.csv)
> A CSV file in the correct input format for the `breseq_plotting.R` program. 

## Step 2 `breseq_plotting.R`

#### Inputs
The `breseq_plotting.R` program two required arguments.

* -f/--file: The data file, plot_variants.csv, created by the `breseq_HTML_to_CSV.py` program
* -o/--plotting_options: The CSV file passing the customizations for the plot. All settings must be set. The example template serves as a good starting point. [example plotting options input](https://github.com/livkosterlitz/Breseq_genome_plots/blob/main/BasicRun/Breseq_plots/plotting_options.csv)

#### Command to run example
```bash
Rscript /Users/oliviakosterlitz/Programs/breseq_plotting.R -f plot_variants.csv -o plotting_options.csv
```
#### Outputs
A single plot in PDF format [example plot](https://github.com/livkosterlitz/Breseq_genome_plots/blob/main/BasicRun/Breseq_plots/out.pdf)
