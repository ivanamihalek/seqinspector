# Installation

No special installation required. Put raw data (fastq and bed files) in directory
called, well, ras_data. 


The scripts, in python folder, were tested to run under python3.8, 
in Linux environment. They can be executed from anywhere, but the toy_pipeline
directory should be of the format as given here, and open for read/write.
If you make the scripts executable you can run them by providing their name on the
command line, e.g. ./python/01_fastq_quality_check.py.

# Dependencies

Samtools, bcftools, cutadapt, as described in the report. Set the paths in 
config.py.

# Run

The pipeline should be executed in the order in which it is enumerated. The
numbers may not be consecutive, it does not reflect the completeness of the pipeline.
The scripts not prefixed by a number contain utility methods.