# gatk-demuxlet-workflow
A workflow performs GATK variant calling on one or multiple RNASeq bam(s) based on GATK best practice. 

If 'mode' == 'demuxlet', demuxlet (Kang. Nature 2017) will be used to perform de-multiplexing on a multiplexed single cell RNASeq BAM using VCF from last step.

![alt text](https://raw.githubusercontent.com/yh154/gatk-demuxlet-workflow/master/workflow.png)

