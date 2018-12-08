This software uses third party toolkits such as GATK (https://software.broadinstitute.org/gatk/) and ANNOVAR annotators to provide annotated variant files from Burrows-Wheeler aligned human NGS data sequenced with Illumina sequencers such as Illumina Miniseq, NextSeq and HiSeq on TrusightOne and TrusightCancer panels, although other panels can be used.

This is an open source program (licensed GNU GPL v.3.0) designed to run on UNIX systems (tested on Ubuntu). All third party toolkits must be installed manually. The program takes in Illumina NGS sequence reads in BAM format and variant files in VCF format.

The pipeline can aggregate all runs in a local environment to provide a lab-specific variant database. Custom annotations are provided from local representations of OMIM, HPO, HGNC data that can be easily downloaded and updated with an update tool provided you have API keys from the databases. HGNC gene symbols or whole gene panels (custom or from Genomics England PanelApp (https://panelapp.genomicsengland.co.uk/)) can be used to filter variant tables and create gene specific coverage tables for samples.

Currently the program expects input from either TrusightOne or TrusightCancer Illumina panels. These two pipelines vary slightly and thus have different entry points.

The benefits of using this pipeline manager are
1. ease-of-use after setup
2. running in single sample, multiple sample or batch mode
3. annotated variant files in excel format
4. automatic excel filtering based on set preconditions (highlight variants with an ExAC frequency <0.05)
5. gene specific custom annotations for variants (HPO phenotype terms, diseases related to the gene)
6. conversion of synonymous gene symbols to HGNC symbols automatically
7. gene or gene panel specific filtering based on sample specific order (e.g. prioritize variants corresponding to Intellectual Disability panel-specific genes)

Drawbacks include
1. terminal program
2. not being very easily customised (some knowledge of python required).
3. requires third party databases and programs
4. relies on updates as GATK and annotators evolve

The pipeline consists of
0. Loading gene panel data internally
1. Checking the validity of the command (all samples must have both the BAM and VCF file available) and creating a sample specific order
2. combining order specific VCF files for analysis heterogeneity analysis using VCF tools
3. copying the VCFs to a "database" (GATK CombineVariants of all previous VCFs)
4. third party annotation using preset settings with ANNOVAR
5. try converting all gene names to their HGNC declared symbols (synonym to HGNC)
6. custom annotations added to the ANNOVAR-annotated VCF files (custom annotations include OMIM, HPO terms, order-specific data for filtering.
7. Coverage analysis (if selected) with GATK DepthOfCoverage and DiagnoseTargets based on the order
8. converting table files to an excel file and running post-scripts to edit the excel document (adding autofilters, hyperlinks, converting text to be more readable, etc [all variants are conserved])


Online Mendelian Inheritance in Man, OMIM®. McKusick-Nathans Institute of Genetic Medicine, Johns Hopkins University (Baltimore, MD), 2018. World Wide Web URL: https://omim.org/

Genomics England PanelApp; https://panelapp.genomicsengland.co.uk

HUGO Gene Nomenclature Committee at the European Bioinformatics Institute
