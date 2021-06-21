#! /usr/bin/env nextflow

nextflow.enable.dsl=2 

params.tmp_vcf="*.vcf"
params.tmp_gen="*.fasta"
params.tmp_ill="*_{1,2}.fastq.gz"

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

process filter_freebayes_vcf {
    input: tuple path(genome), path(vcf)
    output: path("${vcf}.filtered_sorted.vcf.gz"), path("*")
    script:
    """
    #! /usr/bin/env bash
    PROC=\$((`nproc`-4))
    bcftools view -Ou -e'type="ref"' ${vcf}.vcf | bcftools norm -Ob -f ${genome} -o ${vcf}_norm.bcf --threads \$PROC
    bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=\$PROC ${vcf}_norm.bcf  > ${vcf}.filtered.vcf.gz
    bcftools sort -o ${vcf}.filtered_sorted.vcf.gz -m 10G -Oz ${vcf}.filtered.vcf.gz
    bcftools index ${vcf}.filtered_sorted.vcf.gz
    """
}

process jellyfish {
    input: tuple path(vcf), path(reads)
    output: val("\$peak")
    script:
    """
    #! /usr/bin/env bash
    # run jellyfish to get kmer data
    jellyfish count -C -m 21 -s 3000000000 -t 30 ${reads} -o ${vcf}reads.jf
    # create histogram from that data
    jellyfish histo -t 36 ${vcf}reads.jf > ${vcf}reads.hist
    # identify peak
    head -n 50 ${vcf}reads.hist
    peak=`more +5 ${vcf}reads.hist | sort -k 2n | tail -n 1 | awk '{print $1}'`
    """
}

process mk_meryldb {
    input: tuple path(genome), path(r1), path(r2)
    output: tuple path("$genome"), path("${vcf}all.meryl")
    script:
    """
    #! /usr/bin/env bash 
    ${meryl} count k=$kmer ${genome}  output ${genome}.meryl
    ${meryl} count k=$kmer ${readR1} output ${readR1}.meryl
    ${meryl} count k=$kmer ${readR2} output ${readR2}.meryl
    ${meryl} union-sum ${readR1}.meryl ${readR2}.meryl output ${vcf}all.meryl
    """
}

process merfin {
    input: tuple path(genome), path(genome_meryl), path(vcf), val(peak)
    output: path("merfin_filtered.vcf")
    script:
    """
    #! /usr/bin/env bash
    # figure out this section
    ${merfin} -vmer \
     -sequence ${genome} \
     -seqmers ${genome}.meryl \
     -readmers ${vcf}all.meryl \
     -peak ${peak} \
     -vcf ${vcf}.filtered_sorted.vcf.gz   \
     -output ${vcf}_out.dump.gz
    """
}

process build_consensus {
    input:
    output:
    script:
    """
    #! /usr/bin/env bash
    bcftools view -Oz ${vcf}_out.dump.gz.polish.vcf > ${vcf}_out.dump.gz.polish.vcf.gz               # bgzip merfin output
    bcftools index ${vcf}_out.dump.gz.polish.vcf.gz
    bcftools consensus ${vcf}_out.dump.gz.polish.vcf.gz -f ${genome} -H 1 > ${genome}_polished.fasta # -H 1 applies only first allele from GT at each positionscrips
    """
}


workflow{ 
    vcf_ch = channel.fromPath(params.tmp_vcf, checkIfExists:true)
    gen_ch = chennel.fromPath(params.tmp_gen, checkIfExists:true)
    ill_ch = channel.fromPath(params.tmp_ill, checkIfExists:true)

    // From Arrow
    peak_ch = vcf_ch | reshape_arrow | combine(ill_ch) | jellyfish 
    gen_ch | mk_meryldb | combine(reshape_arrow.out) | combine(jellyfish.out) | merfin
    
    // From Freebayes
    peak_ch = vcf_ch | filter_freebayes_vcf | combine(ill_ch) | jellyfish
    gen_ch | mk_meryldb | combine(filter_freebayes_vcf.out) | combine(jellyfish.out) | merfin
}