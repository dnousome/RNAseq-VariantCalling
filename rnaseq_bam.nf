#!/usr/bin/env nextflow
//Variant Calling for RNA-seq BAM files
nextflow.enable.dsl=2

date = new Date().format( 'yyyyMMdd' )

GENOME=file(params.genome)
GENOMEDICT=file(params.genomedict)
WGSREGION=file(params.wgsregion) 
MILLSINDEL=file(params.millsindel) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"// file(params.gold_indels1) //
SHAPEITINDEL=file(params.shapeitindel) //params.shapeitindel =  "/data/OpenOmics/references/genome-seek/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz" //file(params.gold_indels2) //
KGP=file(params.kgp) ///data/CCBR_Pipeliner/Exome-seek/hg38/GATK_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
DBSNP=file(params.dbsnp) //= "/data/OpenOmics/references/genome-seek/GATK_resource_bundle/hg38bundle/dbsnp_138.hg38.vcf.gz"
GNOMAD=file(params.gnomad) //= '/data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz' // /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz
PON=file(params.pon) 


results_dir=params.output
baminput= Channel.fromPath(params.bams, checkIfExists: true, type: 'file')
               .map { row -> tuple(
                        row.simpleName,row
                       )}

baiinput= Channel.fromPath(params.bai, checkIfExists: true, type: 'file')
               .map { row -> tuple(
                        row.simpleName,row
                       )}

baminput=baminput.join(baiinput)


intervalbed = Channel.fromPath(params.intervals,checkIfExists: true,type: 'file')
sample_sheet=Channel.fromPath(params.sample_sheet, checkIfExists: true)
                       .ifEmpty { exit 1, "sample sheet not found" }
                       .splitCsv(header:true, sep: "\t")
                       .map { row -> tuple(
                        row.Tumor,
                       row.Normal
                       )
                                  }
//Final Workflow
workflow {
    //baminput.view()
    
    
    bamwithsample=baminput.join(sample_sheet)
    .map{it.swap(3,0)}
    .join(baminput)
    .map{it.swap(3,0)}

    bambyinterval=bamwithsample.combine(intervalbed)
    
    
    //Paired Mutect2    
    mutect2(bambyinterval)
    pileup_paired_t(bambyinterval)
    pileup_paired_n(bambyinterval)
    

    //pileup_paired_t.out.view()
      
    pileup_paired_tout=pileup_paired_t.out.groupTuple()
    .map{samplename,pileups-> tuple( samplename,
    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tumor.pileup.table/)[0][1].toInteger() } ,
    )}

  
    pileup_paired_nout=pileup_paired_n.out.groupTuple()
    .map{samplename,pileups-> tuple( samplename,
    pileups.toSorted{ it -> (it.name =~ /${samplename}_(.*?).normal.pileup.table/)[0][1].toInteger() } ,
    )}

    
    pileup_paired_all=pileup_paired_tout.join(pileup_paired_nout)
    contamination_paired(pileup_paired_all)
    
    mut2out_lor=mutect2.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    f1r2.toSorted{ it -> (it.name =~ /${samplename}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } 
    )}

    learnreadorientationmodel(mut2out_lor)

    mut2out_mstats=mutect2.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    stats.toSorted{ it -> (it.name =~ /${samplename}_(.*?).mut2.vcf.gz.stats/)[0][1].toInteger() } 
    )}

    mergemut2stats(mut2out_mstats)

    allmut2tn=mutect2.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    vcfs.toSorted{ it -> (it.name =~ /${samplename}_(.*?).mut2.vcf.gz/)[0][1].toInteger() } 
    )}
    
    mut2tn_filter=allmut2tn
    .join(mergemut2stats.out)
    .join(learnreadorientationmodel.out)
    .join(contamination_paired.out)
    mutect2filter(mut2tn_filter)



    //Tumor Only Calling

    mutect2_tonly(bambyinterval)    
    
    //LOR     
    mut2tout_lor=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    f1r2.toSorted{ it -> (it.name =~ /${samplename}_(.*?).f1r2.tar.gz/)[0][1].toInteger() } 
    )}
    learnreadorientationmodel_tonly(mut2tout_lor)


    //Stats
    mut2tonly_mstats=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    stats.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tonly.mut2.vcf.gz.stats/)[0][1].toInteger() } 
    )}
    mergemut2stats_tonly(mut2out_mstats)


    //Contamination
    contamination_tumoronly(pileup_paired_tout)


    
    //Final TUMOR ONLY FILTER
    allmut2tonly=mutect2_tonly.out.groupTuple()
    .map { samplename,vcfs,f1r2,stats -> tuple( samplename,
    vcfs.toSorted{ it -> (it.name =~ /${samplename}_(.*?).tonly.mut2.vcf.gz/)[0][1].toInteger() } 
    )}
    
    mut2tonly_filter=allmut2tonly
    .join(mergemut2stats_tonly.out)
    .join(learnreadorientationmodel_tonly.out)
    .join(contamination_tumoronly.out)

    mutect2filter_tonly(mut2tonly_filter)


    //#To implement
        //CNMOPs from the BAM BQSRs
        //##VCF2MAF TO
    //inpvep=vcfinput.join(sample_sheet)
    //annotvep(inpvep)
    
    tn_vepin=mutect2filter.out
    .join(sample_sheet)

   annotvep_tn(tn_vepin)
   annotvep_tonly(mutect2filter_tonly.out)
    
}


process mutect2 {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz.stats")

    
    script:

    """
    gatk Mutect2 \
    --reference ${GENOME} \
    --intervals ${bed} \
    --input ${tumor} \
    --input ${normal} \
    --normal-sample ${normal.simpleName} \
    --tumor-sample ${tumor.simpleName} \
    --germline-resource ${GNOMAD} \
    --panel-of-normals ${PON} \
    --output ${tumor.simpleName}_${bed.simpleName}.mut2.vcf.gz \
    --f1r2-tar-gz ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates
    """
}
   
/*
//--germline-resource /data/CCBR_Pipeliner/Exome-seek/hg38/GNOMAD/somatic-hg38-af-only-gnomad.hg38.vcf.gz --panel-of-normals /data/CCBR_Pipeliner/Exome-seek/hg38/PON/hg38.noCOSMIC_ClinVar.pon.vcf.gz'"]
     
"--genome-reference /data/CCBR_Pipeliner/Exome-seek/hg38/genome/Homo_sapiens_assembly38.fasta",
     "--output-directory",args.outdir,
     "-dbsnp /data/CCBR_Pipeliner/Exome-seek/hg38/GATK_resource_bundle/dbsnp_138.hg38.vcf.gz --threads 64",
     "--run-mutect2",
     "--mutect2-arguments '
*/  
process mutect2filter {
    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups), path(normal_pileups),path(tumorcontamination),path(normalcontamination)
    output:
        tuple val(sample), path("${sample}.marked.vcf.gz"),path("${sample}.final.mut2.vcf.gz"),path("${sample}.marked.vcf.gz.filteringStats.tsv")
    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")


    """
    gatk GatherVcfs -I ${mut2in} -O ${sample}.concat.vcf.gz 
    gatk IndexFeatureFile -I ${sample}.concat.vcf.gz 
    gatk FilterMutectCalls \
        -R ${GENOME} \
        -V ${sample}.concat.vcf.gz \
        --ob-priors ${obs} \
        --contamination-table ${tumorcontamination} \
        --stats ${stats} \
        -O ${sample}.marked.vcf.gz


    gatk SelectVariants \
        -R ${GENOME} \
        --variant ${sample}.marked.vcf.gz \
        --exclude-filtered \
        --output ${sample}.final.mut2.vcf.gz
    """
}



process mutect2_tonly {
    input:
        tuple val(tumorname), path(tumor), path(tumorbai),val(normalname), path(normal), path(normalbai), path(bed)
    
    output:
        tuple val(tumorname),
        path("${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz"),
        path("${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz.stats")
    
    script:

    """
    gatk Mutect2 \
    --reference ${GENOME} \
    --intervals ${bed} \
    --input ${tumor} \
    --tumor-sample ${tumor.simpleName} \
    --germline-resource ${GNOMAD} \
    --panel-of-normals ${PON} \
    --output ${tumor.simpleName}_${bed.simpleName}.tonly.mut2.vcf.gz \
    --f1r2-tar-gz ${tumor.simpleName}_${bed.simpleName}.f1r2.tar.gz \
    --independent-mates    
    """
}



process mutect2filter_tonly {
    input:
        tuple val(sample), path(mutvcfs), path(stats), path(obs), path(pileups),path(tumorcontamination)
    output:
        tuple val(sample), path("${sample}.tonly.marked.vcf.gz"),path("${sample}.tonly.final.mut2.vcf.gz"),path("${sample}.tonly.marked.vcf.gz.filteringStats.tsv")
    script:
    //Include the stats and  concat ${mutvcfs} -Oz -o ${sample}.concat.vcf.gz
    mut2in = mutvcfs.join(" -I ")


    """
    gatk GatherVcfs -I ${mut2in} -O ${sample}.concat.vcf.gz 
    gatk IndexFeatureFile -I ${sample}.concat.vcf.gz 
    gatk FilterMutectCalls \
        -R ${GENOME} \
        -V ${sample}.concat.vcf.gz \
        --ob-priors ${obs} \
        --contamination-table ${tumorcontamination} \
        --stats ${stats} \
        -O ${sample}.tonly.marked.vcf.gz


    gatk SelectVariants \
        -R ${GENOME} \
        --variant ${sample}.tonly.marked.vcf.gz \
        --exclude-filtered \
        --output ${sample}.tonly.final.mut2.vcf.gz
    """
}




process annotvep_tn {
    module=['vcf2maf/1.6.21','VEP/106']
    
    publishDir("${results_dir}/mafs/", mode: "copy")

    input:
        tuple val(tumorsample), 
        path("${tumorsample}.marked.vcf.gz"), 
        path("${tumorsample}.final.mut2.vcf.gz"), 
        path("${tumorsample}.marked.vcf.gz.filteringStats.tsv"), 
        val(normalsample)

    output:
        path("${tumorsample}.maf")

    script:

    """
    
    zcat ${tumorsample}.final.mut2.vcf.gz > ${tumorsample}.final.mut2.vcf

    vcf2maf.pl \
    --vep-forks 16 --input-vcf ${tumorsample}.final.mut2.vcf \
    --output-maf ${tumorsample}.maf \
    --tumor-id ${tumorsample} \
    --normal-id ${normalsample} \
    --vep-path \${VEP_HOME}/bin \
    --vep-data \${VEP_CACHEDIR} \
    --ncbi-build GRCh38 --species homo_sapiens --ref-fasta ${GENOME}

    """
}


process annotvep_tonly {
    module=['vcf2maf/1.6.21','VEP/106']

    publishDir("${results_dir}/mafs/", mode: "copy")

    input:
        tuple val(tumorsample), 
        path("${tumorsample}.tonly.marked.vcf.gz"),
        path("${tumorsample}.tonly.final.mut2.vcf.gz"),
        path("${tumorsample}.tonly.marked.vcf.gz.filteringStats.tsv")

    output:
        path("${tumorsample}.tonly.maf")

    script:

    """
    
    zcat ${tumorsample}.tonly.final.mut2.vcf.gz  > ${tumorsample}.tonly.final.mut2.vcf

    vcf2maf.pl \
    --vep-forks 16 --input-vcf ${tumorsample}.tonly.final.mut2.vcf \
    --output-maf ${tumorsample}.tonly.maf \
    --tumor-id ${tumorsample} \
    --vep-path \${VEP_HOME}/bin \
    --vep-data \${VEP_CACHEDIR} \
    --ncbi-build GRCh38 --species homo_sapiens --ref-fasta ${GENOME}

    """
}