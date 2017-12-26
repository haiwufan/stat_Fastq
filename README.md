# stat_Fastq
Statistic illumina Fastq file and Plot 

version:  1.0.2
Author:   haibo.li

Install:
  cd <source dir>
  git clone https://github.com/haiwufan/stat_Fastq.git
  cd stat_Fastq
  javac stat_Fastq.java
  
set Environment:
  export CLASSPATH=<source>/stat_Fastq
 
 usage:
  java stat_Fastq <33|64> <Output_prefix> <file1.fq|file1.fq.gz> [file2.fq|file2.fq.gz] [...]
  
  Rscript plot.R <sample_name> <R1_file_prefix> <R2_file_prefix>
  
  e.g:
    java stat_Fastq 33 Sample1 tmp.R1.fq tmp.R2.fq
    Rscript plot.R Sample1 tmp.R1 tmp.R2
  
 
  
