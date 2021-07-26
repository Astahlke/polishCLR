#! /usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
  log.info isuGIFHeader()
  log.info """
   Usage:
   The typical command for running the pipeline are as follows:
   nextflow run main.nf --primary_assembly "*fasta" --illumina_reads "*{1,2}.fastq.bz2" --pacbio_reads "*_subreads.bam" -resume

   Mandatory arguments:
   --primary_assembly             genome assembly fasta file to polish
   --illumina_reads               paired end illumina reads, to be used for Merqury QV scores, and freebayes polish primary assembly
   --pacbio_reads                 pacbio reads in bam format, to be used to arrow polish primary assembly

   Optional modifiers
   --k                            kmer to use in MerquryQV scoring [default:21]
   --same_specimen                if illumina and pacbio reads are from the same specimin [default: true].
   --falcon_unzip                 if primary assembly has already undergone falcon unzip [default: false]. If true, will Arrow polish once instead of twice.

   Optional configuration arguments
   --parallel_app                 Link to parallel executable [default: 'parallel']
   --bzcat_app                    Link to bzcat executable [default: 'bzcat']
   --pigz_app                     Link to pigz executable [default: 'pigz']
   --meryl_app                    Link to meryl executable [default: 'meryl']
   --merqury_sh                   Link to merqury script [default: '\$MERQURY/merqury.sh']
   --pbmm2_app                    Link to pbmm2 executable [default: 'pbmm2']
   --samtools_app                 Link to samtools executable [default: 'samtools']
   --gcpp_app                     Link to gcpp executable [default: 'gcpp']
   --bwamem2_app                  Link to bwamem2 executable [default: 'bwa-mem2']
   --freebayes_app                Link to freebayes executable [default: 'freebayes']
   --bcftools_app                 Link to bcftools executable [default: 'bcftools']

   Optional arguments:
   --outdir                       Output directory to place final output [default: 'PolishCLR_Results']
   --clusterOptions               Cluster options for slurm or sge profiles [default slurm: '-N 1 -n 40 -t 04:00:00'; default sge: ' ']
   --threads                      Number of CPUs to use during each job [default: 40]
   --queueSize                    Maximum number of jobs to be queued [default: 50]
   --account                      Some HPCs require you supply an account name for tracking usage.  You can supply that here.
   --help                         This usage statement.
  """
}
//   --bzip                         if illumina paired reads are in bz2 format [default: true]. If true, will convert to gz.

// Show help message
if (params.help || !params.primary_assembly || !params.illumina_reads || !params.pacbio_reads ) {
  helpMessage()
  exit 0
}

// 00 Preprocess: Unzip any bz2 files
process bz_to_gz {
    publishDir "${params.outdir}/00_Preprocess", mode: 'copy'
    input:tuple val(readname), path(illumina_reads)
    output: tuple val(readname), path("*.gz")
    script:
    template 'bz_to_gz.sh'
}

// 01 Merqury: Quality value of primary assembly
process meryl_count {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'symlink'
    input: tuple val(k), path(illumina_read)
    output: path("*.meryl")
    script:
    template 'meryl_count.sh'
}

process meryl_union {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'copy'
    input: path(illumina_meryls)
    output: path("illumina.meryl")
    script:
    template 'meryl_union.sh'
}

process MerquryQV_01 {
    publishDir "${params.outdir}/01_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// 1st Arrow Polish
process create_windows {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(assembly_fasta)
    output: tuple path("*.fai"), path("win.txt")
    script:
    template 'create_windows.sh'
}

process pbmm2_index {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(assembly_fasta)
    output: tuple val("$i"), path("$assembly_fasta"), path("*.mmi")
    script:
    template 'pbmm2_index.sh'
}

process pbmm2_align {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input:tuple val(i), path(assembly_fasta), path(assembly_mmi), path(pacbio_read)
    output: tuple val("$i"), path("*.bam"), path("*.bai")
    script:
    template 'pbmm2_align.sh'
}

process gcc_Arrow {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'symlink'
    input: tuple val(i), path(pacbio_bam), path(pacbio_bai),  path(assembly_fasta), path(assembly_fai), val(window)
    output: tuple val("$i"), path("*.fasta"), path("*.vcf")
    script:
    template 'gcc_arrow.sh'
}

process merge_consensus {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'copy'
    input: tuple val(i), path(windows_fasta)
    output: path("${i}_consensus.fasta")
    script:
    template 'merge_consensus.sh'
}

workflow ARROW_02 {
  take:
    asm_ch
    pac_ch
  main:
    win_ch = channel.of("2") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    newasm_ch = channel.of("2") | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow | 
      map { n -> [ n.get(0), n.get(1) ] } | groupTuple | merge_consensus
  
  emit:
    newasm_ch
}

process MerquryQV_03 {
    publishDir "${params.outdir}/03_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// 2nd Arrow run with merfin

process merfin {
    publishDir "${params.outdir}/0${i}_ArrowPolish", mode: 'copy'
    input: tuple val(i), path(vcf), path(asm)
    output: path("${i}_reshaped.reshaped.vcf.gz")
    script:
    template 'merfin.sh'
}

workflow ARROW_04 {
  take:
    asm_ch
    pac_ch
  main:
    win_ch = channel.of("4") | combine(asm_ch) | create_windows | 
      map { n -> n.get(1) } | splitText() {it.trim() }
    fai_ch = create_windows.out | map { n -> n.get(0) }

    newasm_ch = channel.of("4") | combine(asm_ch) | pbmm2_index | combine(pac_ch) | pbmm2_align |
      combine(asm_ch) | combine(fai_ch) | combine(win_ch) | gcc_Arrow | 
      map { n -> [ n.get(0), n.get(2) ] } | groupTuple |
      combine(asm_ch) |  merfin  //| 
//      groupTuple | merge_consensus
  
  emit:
    newasm_ch
}

process MerquryQV_05 {
    publishDir "${params.outdir}/05_MerquryQV", mode: 'copy'
    input: tuple path(illumina_db), path(assembly_fasta)
    output: path("*")
    script:
    template 'merquryqv.sh'
}

// // 1st FreeBayes Polish
// process align_shortreads {
//     publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
//     input: tuple path(assembly_fasta), path(illumina_one), path(illumina_two)
//     output: tuple path("*.bam"), path("*.bai")
//     script:
//     """
//     #! /usr/bin/env bash
//     PROC=\$((`nproc` /2+1))
//     mkdir tmp
//     ${bwamem2_app} index ${assembly_fasta}
//     ${bwamem2_app} mem -SP -t \${PROC} ${assembly_fasta} ${illumina_one} ${illumina_two} |
//       ${samtools_app} sort -T tmp -m 8G --threads 4 - > ${illumina_one.simpleName}_aln.bam
//     ${samtools_app} index -@ \${PROC} ${illumina_one.simpleName}_aln.bam
//     """
// }
// // -m 5G -@ 36
//
// process freebayes {
//     publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
//     input: tuple val(window), path(assembly_fasta), path(assembly_fai), path(illumina_bam), path(illumina_bai)
//     output: path("*.vcf")
//     script:
//     """
//     #! /usr/bin/env bash
//     ${freebayes_app} \
//       --region \"${window}\" \
//       --min-mapping-quality 0 \
//       --min-coverage 3 \
//       --min-supporting-allele-qsum 0 \
//       --ploidy 2 \
//       --min-alternate-fraction 0.2 \
//       --max-complex-gap 0 \
//       --bam ${illumina_bam} \
//       --vcf ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}.vcf \
//       --fasta-reference ${assembly_fasta}
//     """
// }
//
// process combineVCF {
//     publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
//     input: path(vcfs)
//     output: path("consensus_04.vcf")
//     script:
//     """
//     #! /usr/bin/env bash
//     cat ${vcfs.get(0)} | grep "^#" > consensus_04.vcf
//     cat ${vcfs} | grep -v "^#" >> consensus_04.vcf
//     """
// }
//
// process vcf_to_fasta {
//     publishDir "${params.outdir}/04_FreeBayesPolish", mode: 'symlink'
//     input: tuple path(vcf), path(genome_fasta)
//     output: path("consensus_04.fasta")
//     script:
//     """
//     #! /usr/bin/env bash
//     ${bcftools_app} view -Oz ${vcf} > ${vcf.simpleName}.gz
//     ${bcftools_app} index ${vcf.simpleName}.gz
//     ${bcftools_app} consensus ${vcf.simpleName}.gz -f ${genome_fasta} -H 1 > consensus_04.fasta
//     """
// }
//
// workflow Freebayes_06 {
//   take:
//     assembly_ch
//     illumina_ch
//   main:
//     assembly_ch | combine(illumina_ch.collect()) | align_shortreads
//     fai_ch = assembly_ch | create_windows | map { n -> n.get(0) }
//     create_windows.out |
//       map { n -> n.get(1) } |
//       splitText() {it.trim()} |
//       combine(assembly_ch) | combine(fai_ch) | combine(align_shortreads.out) |
//       freebayes |
//       collect |
//       combineVCF |
//       combine(assembly_ch) |
//       vcf_to_fasta
//     new_asm_ch = vcf_to_fasta.out         // <= New assembly
//   emit:
//     new_asm_ch
// }
//
//
// // 3rd Merqury QV value
// process MerquryQV_05 {
//     publishDir "${params.outdir}/05_MerquryQV", mode: 'symlink'
//     input: tuple path(illumina_db), path(assembly_fasta)
//     output: path("*")
//     script:
//     """
//     #! /usr/bin/env bash
//     $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
//     """
// }
//
// // 2nd Freebayes polish
// process align_shortreads_06 {
//     publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
//     input: tuple path(assembly_fasta), path(illumina_one), path(illumina_two)
//     output: tuple path("*.bam"), path("*.bai")
//     script:
//     """
//     #! /usr/bin/env bash
//     PROC=\$((`nproc` /2+1))
//     mkdir tmp
//     bwa-mem2 index ${assembly_fasta}
//     bwa-mem2 mem -SP -t \${PROC} ${assembly_fasta} ${illumina_one} ${illumina_two} |
//       samtools sort -T tmp -m 8G --threads 4 - > ${illumina_one.simpleName}_06_aln.bam
//     samtools index -@ \${PROC} ${illumina_one.simpleName}_06_aln.bam
//     """
// }
// // -m 5G -@ 36
//
// process create_windows_06 {
//     publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
//     input: path(assembly_fasta)
//     output: tuple path("*.fai"), path("win_06.txt")
//     shell:
//     """
//     #! /usr/bin/env bash
//     samtools faidx ${assembly_fasta}
//     cat ${assembly_fasta}.fai | awk '{print \$1 ":0-" \$2}' > win_06.txt
//     """
// }
//
// process freebayes_06 {
//     publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
//     input: tuple val(window), path(assembly_fasta), path(assembly_fai), path(illumina_bam), path(illumina_bai)
//     output: path("*.vcf")
//     script:
//     """
//     #! /usr/bin/env bash
//     freebayes \
//       --region \"${window}\" \
//       --min-mapping-quality 0 \
//       --min-coverage 3 \
//       --min-supporting-allele-qsum 0 \
//       --ploidy 2 \
//       --min-alternate-fraction 0.2 \
//       --max-complex-gap 0 \
//       --bam ${illumina_bam} \
//       --vcf ${illumina_bam.simpleName}_${window.replace(':','_').replace('|','_')}_06.vcf \
//       --fasta-reference ${assembly_fasta}
//     """
// }
//
// process combineVCF_06 {
//     publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
//     input: path(vcfs)
//     output: path("consensus_06.vcf")
//     script:
//     """
//     #! /usr/bin/env bash
//     cat ${vcfs.get(0)} | grep "^#" > consensus_06.vcf
//     cat ${vcfs} | grep -v "^#" >> consensus_06.vcf
//     """
// }
//
// process vcf_to_fasta_06 {
//     publishDir "${params.outdir}/06_FreeBayesPolish", mode: 'symlink'
//     input: tuple path(vcf), path(genome_fasta)
//     output: path("final_polished_assembly.fasta")
//     script:
//     """
//     #! /usr/bin/env bash
//     bcftools view -Oz ${vcf} > ${vcf}.gz
//     bcftools index ${vcf.simpleName}.vcf.gz
//     bcftools consensus ${vcf.simpleName}.vcf.gz -f ${genome_fasta} -H 1 > final_polished_assembly.fasta
//     """
// }
//
// // 3rd Merqury QV value
// process MerquryQV_07 {
//     publishDir "${params.outdir}/07_MerquryQV", mode: 'symlink'
//     input: tuple path(illumina_db), path(assembly_fasta)
//     output: path("*")
//     script:
//     """
//     #! /usr/bin/env bash
//     $MERQURY/merqury.sh $illumina_db $assembly_fasta ${assembly_fasta.simpleName}
//     """
// }

workflow {
    // Setup input channels, starting assembly (asm), Illumina reads (ill), and pacbio reads (pac)
    asm_ch = channel.fromPath(params.primary_assembly, checkIfExists:true)
    ill_ch = channel.fromFilePairs(params.illumina_reads, checkIfExists:true)
    pac_ch = channel.fromPath(params.pacbio_reads, checkIfExists:true)
    k_ch   = channel.of(params.k) // Either passed in or autodetect (there's a script for this)

    // Step 0: Preprocess illumina files from bz2 to gz files
    // Instead of a flag, auto detect, however it must be in the pattern, * will fail
    if(params.illumina_reads =~ /bz2$/){
      pill_ch = ill_ch | bz_to_gz | map { n -> n.get(1) } | flatten
    }else{
      pill_ch = ill_ch | map { n -> n.get(1) } | flatten
    }

    // Step 1: Check quality of assembly with Merqury
    merylDB_ch = k_ch | combine(pill_ch) | meryl_count | collect | meryl_union 
    merylDB_ch | combine(asm_ch) | MerquryQV_01
    
    // Step 2: Arrow Polish with PacBio reads
    asm_arrow_ch = ARROW_02(asm_ch, pac_ch)
    asm_arrow_ch | view
    //
    // Step 3: Check quality of new assembly with Merqury (turns out we can reuse the illumina database)
    merylDB_ch | combine(asm_arrow_ch) | MerquryQV_03

    // if the primary assembly came from falcon unzip, skip the 2nd arrow polish
    if(!params.falcon_unzip) {
      // Step 2b: Arrow Polish with PacBio reads
      asm_arrow2_ch = ARROW_04(asm_arrow_ch, pac_ch)
      // Step 3b: Check quality of new assembly with Merqury (turns out we can reuse the illumina database)
      merylDB_ch | combine(asm_arrow2_ch) | MerquryQV_05
    } else {
      asm_arrow2_ch = asm_arrow_ch
    }
    asm_arrow2_ch | view
    //
    // // Step 4: FreeBayes Polish with Illumina reads
    //
    // asm2b_ch | combine(pill_ch.collect()) | align_shortreads_04 | view
    // fai2_ch = asm2b_ch | create_windows_04 | map { n -> n.get(0) } | view
    // create_windows_04.out |
    //   map { n -> n.get(1) } |
    //   splitText() {it.trim()} |
    //   combine(asm2b_ch) | combine(fai2_ch) | combine(align_shortreads_04.out) |
    //   freebayes_04 |
    //   collect |
    //   combineVCF_04 |
    //   combine(asm2b_ch) |
    //   vcf_to_fasta_04
    // asm3_ch = vcf_to_fasta_04.out         // <= New assembly
    //
    // // Step 5: Check quality of assembly with Merqury
    // meryl_union_01.out | combine(asm3_ch) | MerquryQV_05
    //
    // // Step 6: FreeBayes Polish 2nd time
    // asm3_ch | combine(pill_ch.collect()) | align_shortreads_06 | view
    // fai3_ch = asm3_ch | create_windows_06 | map { n -> n.get(0) } | view
    // create_windows_06.out |
    //   map { n -> n.get(1) } |
    //   splitText() {it.trim()} |
    //   combine(asm3_ch) | combine(fai3_ch) | combine(align_shortreads_06.out) |
    //   freebayes_06 |
    //   collect |
    //   combineVCF_06 |
    //   combine(asm3_ch) |
    //   vcf_to_fasta_06
    // asm4_ch = vcf_to_fasta_06.out         // <= New assembly
    //
    // // Step 7: Check quality of assembly with Merqury
    // meryl_union_01.out | combine(asm4_ch) | MerquryQV_07
}

def isuGIFHeader() {
  // Log colors ANSI codes
  c_reset = params.monochrome_logs ? '' : "\033[0m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_black = params.monochrome_logs ? '' : "\033[1;90m";
  c_green = params.monochrome_logs ? '' : "\033[1;92m";
  c_yellow = params.monochrome_logs ? '' : "\033[1;93m";
  c_blue = params.monochrome_logs ? '' : "\033[1;94m";
  c_purple = params.monochrome_logs ? '' : "\033[1;95m";
  c_cyan = params.monochrome_logs ? '' : "\033[1;96m";
  c_white = params.monochrome_logs ? '' : "\033[1;97m";
  c_red = params.monochrome_logs ? '' :  "\033[1;91m";

  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_white}                                ${c_red   }\\\\------${c_yellow}---//       ${c_reset}
  ${c_white}  ___  ___        _   ___  ___  ${c_red   }  \\\\---${c_yellow}--//        ${c_reset}
  ${c_white}   |  (___  |  | / _   |   |_   ${c_red   }    \\-${c_yellow}//         ${c_reset}
  ${c_white}  _|_  ___) |__| \\_/  _|_  |    ${c_red  }    ${c_yellow}//${c_red  } \\        ${c_reset}
  ${c_white}                                ${c_red   }  ${c_yellow}//---${c_red  }--\\\\       ${c_reset}
  ${c_white}                                ${c_red   }${c_yellow}//------${c_red  }---\\\\       ${c_reset}
  ${c_cyan}  isugifNF/polishCLR  v${workflow.manifest.version}       ${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}
