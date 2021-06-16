#! /usr/bin/env nextflow

nextflow.enable.dsl=2 

params.tmp_vcf
params.tmp_gen

process reshape_arrow {
    input: path(vcf)
    output: path("${vcf.simpleName}.reshaped.vcf.gz")
    script:
    """
    #! /usr/bin/env bash
    grep -v "#" ${vcf} | sed 's/,/;/g' > ${vcf.simpleName}.temp.reshaped.vcf
    bcftools view -h ${vcf} > ${vcf.simpleName}.temp.reshaped.header.vcf
    cat ${vcf.simpleName}.temp.reshaped.header.vcf ${vcf.simpleName}.temp.reshaped.vcf > ${vcf.simpleName}.temp.reshaped.combined.vcf
    rm ${vcf.simpleName}.temp.reshaped.header.vcf ${vcf.simpleName}.temp.reshaped.vcf
    bcftools annotate -h extra_header.vcf ${vcf.simpleName}.temp.reshaped.combined.vcf > ${vcf.simpleName}.temp.reshaped.vcf
    bcftools view -h ${vcf.simpleName}.temp.reshaped.vcf | sed 's/\tINFO/\tINFO\tFORMAT\tIND/g' > ${vcf.simpleName}.reshaped.vcf
    rm ${vcf.simpleName}.temp.reshaped.vcf 
    bcftools view -H ${vcf.simpleName}.temp.reshaped.combined.vcf | awk -F"\\t" -v OFS="\\t" '{gsub(/DP=/,".\\tGT:DP\\t1/1:",\$8);print \$0}' >> ${vcf.simpleName}.reshaped.vcf
    bcftools view ${vcf.simpleName}.reshaped.vcf -Oz > ${vcf.simpleName}.reshaped.vcf.gz
    rm ${vcf.simpleName}.reshaped.vcf 
    rm ${vcf.simpleName}.temp.reshaped.combined.vcf
    """
}

process merfin {
    input: path(vcf)
    output: path("merfin_filtered.vcf")
    script:
    """
    #! /usr/bin/env bash
    # figure out this section
    merfin -vmer ${vcf} > merfin_filtered.vcf
    """
}
workflow{ 
    vcf_ch = channel.fromPath(params.tmp_vcf, checkIfExists:true)
    gen_ch = chennel.fromPath(params.tmp_gen, checkIfExists:true)

    vcf_ch | reshape_arrow | merfin
}