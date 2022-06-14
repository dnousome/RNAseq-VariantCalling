import os
import pandas as pd
import re
import sys
import glob
import datetime

configfile:"config.json"


BASEDIR=os.path.realpath(config['input_params']['BASE_OUTDIR'])

input_bamdir=os.path.join(BASEDIR,"input_files","bam")
bam_source=os.path.join(config['input_params']['BAM_SOURCE'])

output_bamdir=os.path.join(BASEDIR,config['output_params']['BAM'])

output_vcf=os.path.join(BASEDIR,"vcf")

SAMPLES, = glob_wildcards(config['input_params']['BAM_SOURCE'] + "/{samples}_star.genome.sorted.bam")

gvcfLst = expand(os.path.join(output_vcf, "{samples}.vcf") , samples=SAMPLES)


rule all:
    input:
        expand(os.path.join(output_bamdir,"rg_bams","{samples}.RG.bam"), samples=SAMPLES),
        expand(os.path.join(output_bamdir,"md_bams","{samples}.md.bam"), samples=SAMPLES),
        expand(os.path.join(output_bamdir,"md_bams","{samples}.md.bam.bai"), samples=SAMPLES),
        expand(os.path.join(output_bamdir,"splitncigar","{samples}_split.out.bam"), samples=SAMPLES),
        expand(os.path.join(output_bamdir,"recal1","{samples}_recal.table"), samples=SAMPLES),
        expand(os.path.join(output_bamdir,"bqsr1","{samples}_recal.pass1.bam"), samples=SAMPLES),
     #   expand(os.path.join(output_bamdir,"recal2","{samples}_recal.table"), samples=SAMPLES),
     #   expand(os.path.join(output_bamdir,"bqsr2","{samples}_recal.pass2.bam"), samples=SAMPLES),
        expand(os.path.join(output_vcf,"{samples}.vcf"), samples=SAMPLES),
        os.path.join(output_vcf,"CombinedGvcfs", "all.g.vcf")

##########Split to separate files?
rule AddRG:
     input:
        bam = os.path.join(bam_source, "{samples}_star.genome.sorted.bam"),
     output:
        RG = os.path.join(output_bamdir, "rg_bams", "{samples}.RG.bam"),
     params: "-RGLB lib1 -RGPL illumina -RGPU {samples} -RGSM {samples}"   
     envmodules:
        'picard/2.27.2'
     shell:"""
          java -jar ${{PICARDJARPATH}}/picard.jar AddOrReplaceReadGroups {params} -I {input.bam} -O {output.RG}
           """ 
           
rule mark_dups:
    input:
       bam = os.path.join(output_bamdir, "rg_bams","{samples}.RG.bam"),
    output:
       dbam = os.path.join(output_bamdir, "md_bams", "{samples}.md.bam"),
    resources:
       mem_mb = 10000
    envmodules:
        'sambamba/0.8.1'
    threads: 4
    shell: """
        sambamba markdup -t {threads} {input.bam} {output.dbam}
          """
         
#rule mark_dups:
#    input:
#        bam = os.path.join(output_bamdir, "rg_bams","{samples}.RG.bam"),
#    output:
#       dbam = os.path.join(output_bamdir, "md_bams", "{samples}.md.bam"),
     #  metric = os.path.join(output_bamdir, "md_bams", "{samples}.metrics.txt"),
#    resources:
#       mem_mb = 10000
#    envmodules:
#        'picard/2.27.2'       
    #shell: """
    #    java -jar ${{PICARDJARPATH}}/picard.jar MarkDuplicates -I {input.bam} -O {output.dbam} --METRICS_FILE {output.metric} --ASSUME_SORT_ORDER coordinate --OPTICAL_DUPLICATE_PIXEL_DISTANCE 100
    #      """

# Index bam file using samtools
rule index:
      input:
         bam = os.path.join(output_bamdir, "md_bams", "{samples}.md.bam")
      output:
         bai = os.path.join(output_bamdir, "md_bams", "{samples}.md.bam.bai")
      envmodules:
         'samtools/1.12'
      shell:"""
            samtools index {input.bam} {output.bai}
            """ 

#Splits N Cigar Reads from bam file
rule splitNcigar:
     input:
        bam = os.path.join(output_bamdir, "md_bams", "{samples}.md.bam")
     output:
        SBam = os.path.join(output_bamdir, "splitncigar", "{samples}_split.out.bam")
     resources:
        mem_mb= 60000
     params: 
        genome = config['references']['GENOME']
     envmodules:
        'GATK/4.2.0.0'
     shell:""" 
           gatk SplitNCigarReads \\
           -R {params.genome} \\
           -I {input.bam} \\
           -O {output.SBam} 
           """ 

# base recalibration
rule BQSR_Pass1:         
     input:
        bam = os.path.join(output_bamdir, "splitncigar", "{samples}_split.out.bam")
     output:
        Recall =os.path.join(output_bamdir, "recal1", "{samples}_recal.table")
     resources:
        mem_mb = 50000
     envmodules:
        'GATK/4.2.0.0'
     params: 
        genome = config['references']['GENOME'],
        known = config['references']['KNOWNRECAL']
     shell:"""
           gatk BaseRecalibrator \\
           -I {input.bam} \\
           -R {params.genome} \\
           {params.known} \\
           -O {output.Recall}
           """


rule ApplyBQSR1:
    input:
        bam = os.path.join(output_bamdir, "splitncigar", "{samples}_split.out.bam"),
        recal = os.path.join(output_bamdir, "recal1", "{samples}_recal.table")
    output:
        Rbam = os.path.join(output_bamdir, "bqsr1", "{samples}_recal.pass1.bam")
    resources:
        mem_mb = 50000
    envmodules:
        'GATK/4.2.0.0'
    params: 
        genome = config['references']['GENOME']
    shell:"""
          gatk ApplyBQSR \\
          -I {input.bam}  \\
          -R {params.genome} \\
          --bqsr-recal-file {input.recal} \\
          -O {output.Rbam}
          """



#Variant Calling 
rule gatk_HaplotypeCaller:
    input:
        bam =os.path.join(output_bamdir, "bqsr1", "{samples}_recal.pass1.bam")
    params:
        genome = config['references']['GENOME'],
    output:
        vcf =  os.path.join(output_vcf, "{samples}.vcf")
    resources:
        mem_mb = 50000
    envmodules:
        'GATK/4.2.0.0'
    shell:"""
           gatk HaplotypeCaller \\
           -R {params.genome} \\
           -I {input.bam} \\
           --dont-use-soft-clipped-bases \\
           -stand-call-conf 20.0  \\
           -O {output.vcf}
           """
           
 #Combine all the gVCFS for joint calling 
rule CombineGvfs:
    input:
        vcfs =  gvcfLst,
    output:
        combined = os.path.join(output_vcf,"CombinedGvcfs", "all.g.vcf")
    params:
        genome = config['references']['GENOME'],
        lst = " --variant " .join(gvcfLst)
    resources:
        mem_mb = 100000
    envmodules:
        'GATK/4.2.0.0'
    shell:"""
           export JAVA_TOOL_OPTIONS=-Xmx140g
           gatk CombineGVCFs \\
           -R {params.genome} \\
           -O {output.combined} \\
           --variant {params.lst}
           """         
