#!/usr/bin/env nextflow
/*
 * Authors:
 * - Ye <yewangfaith@gmail.com>
 *
 */



// Define contigs here!
CONTIG_LIST = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14","15", "16", "17", "18", "19", "20", "21","22", "23", "24", "25", "26", "27", "28", "30", "31", "32", "33", "Z", "W", "MT"]
contigs = Channel.from(CONTIG_LIST)

/*
    Params
*/

params.debug = false
params.tmpdir = "tmp/"
params.email = ""
params.reference = "(required)"
params.email = "yewangfaith@gmail.com"


// Compressed Reference File
File reference = new File("${params.reference}")
if (params.reference != "(required)") {
   reference_handle = reference.getAbsolutePath();
   reference_handle_uncompressed = reference_handle.replace(".gz", "")
} else {
   reference_handle = "(required)"
}

// Debug
if (params.debug == true) {
    println """

        *** Using debug mode ***

    """
    params.fqs = "${workflow.projectDir}/test_data/sample_sheet.tsv"
    params.bamdir = "${params.out}/bam"
    File fq_file = new File(params.fqs);
    params.fq_file_prefix = ""

    // DEBUG Filter thresholds
    min_depth=0
    qual=10
    mq=10
    dv_dp=0.0


} else {
    // The SM sheet that is used is located in the root of the git repo
    params.bamdir = "(required)"
    params.fq_file_prefix = null;
    params.fqs = "sample_sheet.tsv"

    min_depth=10
    qual=30
    mq=40
    dv_dp=0.5

}

File fq_file = new File(params.fqs);

/*
    ==
    UX
    ==
*/

param_summary = '''


     ▄         ▄  ▄▄▄▄▄▄▄▄▄▄▄                         ▄▄        ▄  ▄▄▄▄▄▄▄▄▄▄▄
    ▐░▌       ▐░▌▐░░░░░░░░░░░▌                        ▐░░▌      ▐░▌▐░░░░░░░░░░░▌
    ▐░▌       ▐░▌ ▀▀▀▀█░█▀▀▀▀                        ▐░▌░▌     ▐░▌▐░█▀▀▀▀▀▀▀▀▀
    ▐░▌       ▐░▌     ▐░▌                             ▐░▌▐░▌    ▐░▌▐░▌
    ▐░▌   ▄   ▐░▌     ▐░▌           ▄▄▄▄▄▄▄▄▄▄▄      ▐░▌ ▐░▌   ▐░▌▐░█▄▄▄▄▄▄▄▄▄
    ▐░▌  ▐░▌  ▐░▌     ▐░▌          ▐░░░░░░░░░░░▌      ▐░▌  ▐░▌  ▐░▌▐░░░░░░░░░░░▌
    ▐░▌ ▐░▌░▌ ▐░▌     ▐░▌           ▀▀▀▀▀▀▀▀▀▀▀      ▐░▌   ▐░▌ ▐░▌▐░█▀▀▀▀▀▀▀▀▀
    ▐░▌▐░▌ ▐░▌▐░▌     ▐░▌                             ▐░▌    ▐░▌▐░▌▐░▌
    ▐░▌░▌   ▐░▐░▌ ▄▄▄▄█░█▄▄▄▄                        ▐░▌     ▐░▐░▌▐░▌
    ▐░░▌     ▐░░▌▐░░░░░░░░░░░▌                        ▐░▌      ▐░░▌▐░▌
     ▀▀       ▀▀  ▀▀▀▀▀▀▀▀▀▀▀                         ▀        ▀▀  ▀


''' + """

    parameters              description                    Set/Default
    ==========              ===========                    =======

    --debug                 Set to 'true' to test          ${params.debug}
    --cores                 Regular job cores              ${params.cores}
    --out                   Directory to output results    ${params.out}
    --fqs                   fastq file (see help)          ${params.fqs}
    --fq_file_prefix        fastq prefix                   ${params.fq_file_prefix}
    --reference             Reference Genome (w/ .gz)      ${params.reference}
    --annotation_reference  SnpEff annotation              ${params.annotation_reference}
    --bamdir                Location for bams              ${params.bamdir}
    --tmpdir                A temporary directory          ${params.tmpdir}
    --email                 Email to be sent results       ${params.email}

    HELP: http://andersenlab.org/dry-guide/pipeline-wi/

"""

println param_summary

if (params.reference == "(required)" || params.fqs == "(required)") {

    println """
    The Set/Default column shows what the value is currently set to
    or would be set to if it is not specified (it's default).
    """
    System.exit(1)
}

if (!reference.exists()) {
    println """

    Error: Reference does not exist

    """
    System.exit(1)
}

if (!fq_file.exists()) {
    println """

    Error: fastq sheet does not exist

    """
    System.exit(1)
}


// Read sample sheet
strainFile = new File(params.fqs)

if (params.fq_file_prefix != "") {
    fqs = Channel.from(fq_file.collect { it.tokenize( '\t' ) })
                 .map { SM, ID, LB, fq1, fq2, seq_folder -> [SM, ID, LB, file("${params.fq_file_prefix}/${fq1}"), file("${params.fq_file_prefix}/${fq2}"), seq_folder] }
} else {
    fqs = Channel.from(fq_file.collect { it.tokenize( '\t' ) })
                 .map { SM, ID, LB, fq1, fq2, seq_folder -> [SM, ID, LB, file("${fq1}"), file("${fq2}"), seq_folder] }
}


fqs.into {
    fqs_kmer
    fqs_align
}

/*
    =============
    Kmer counting
    =============
*/
/*
process kmer_counting {

    cpus params.cores

    tag { ID }

    input:
        set SM, ID, LB, fq1, fq2, seq_folder from fqs_kmer
    output:
        file("${ID}.kmer.tsv") into kmer_set

    """
        # fqs will have same number of lines
        export OFS="\t"
        fq_wc="`zcat ${fq1} | awk 'NR % 4 == 0' | wc -l`"
        
        zcat ${fq1} ${fq2} | \\
        fastq-kmers -k 6 | \\
        awk -v OFS="\t" -v ID=${ID} -v SM=${SM} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${ID}.kmer.tsv
    """
}


process merge_kmer {

    publishDir "${params.out}/phenotype", mode: "copy"

    input:
        file("kmer*.tsv") from kmer_set.collect()
    output:
        file("kmers.tsv")

    """
        cat <(echo "kmer\tfrequency\tSM\tID\twc") *.tsv > kmers.tsv
    """

}
*/


/*
    ===============
    Fastq alignment
    ===============

    The output looks strange below,
    but its designed to group like samples together - so leave it!

*/

process perform_alignment {

    cpus params.cores

    tag { ID }

    input:
        set SM, ID, LB, fq1, fq2, seq_folder from fqs_align
    output:
        set val(SM), file("${ID}.bam"), file("${ID}.bam.bai") into fq_bam_set


    """
        bwa mem -t ${task.cpus} -R '@RG\\tID:${ID}\\tLB:${LB}\\tSM:${SM}' ${reference_handle} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=${params.tmpdir} --out=${ID}.bam /dev/stdin
        sambamba index --nthreads=${task.cpus} ${ID}.bam

        if [[ ! \$(samtools view ${ID}.bam | head -n 10) ]]; then
            exit 1;
        fi

    """
}

/*
    ===================================
    Merge - Generate isotype-level BAMs
    ===================================
*/

process merge_bam {

    cpus params.cores

    tag { SM }

    input:
        set SM, bam, index from fq_bam_set.groupTuple()

    output:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") into SM_bam_set
        file("${SM}.picard.sam.markduplicates") into duplicates_set

    """
    count=`echo ${bam.join(" ")} | tr ' ' '\\n' | wc -l`

    if [ "\${count}" -eq "1" ]; then
        ln -s ${bam.join(" ")} ${SM}.merged.bam
        ln -s ${bam.join(" ")}.bai ${SM}.merged.bam.bai
    else
        sambamba merge --nthreads=${task.cpus} --show-progress ${SM}.merged.bam ${bam.join(" ")}
        sambamba index --nthreads=${task.cpus} ${SM}.merged.bam
    fi

    picard MarkDuplicates I=${SM}.merged.bam O=${SM}.bam M=${SM}.picard.sam.markduplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false
    sambamba index --nthreads=${task.cpus} ${SM}.bam
    """
}

SM_bam_set.into {
                  bam_publish;
                  bam_idxstats;
                  bam_stats;
                  bam_coverage;
                  bam_snp_individual;
                  bam_snp_union;
                  bam_telseq;
                  bam_isotype_stats;
                  bam_manta;
                  bam_tiddit;
                  bam_delly;
                  bam_delly_recall;
}

process bam_isotype_stats {

    cpus params.cores

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_isotype_stats

    output:
         file("${SM}.samtools.txt") into SM_samtools_stats_set
         file("${SM}.bamtools.txt") into SM_bamtools_stats_set
         file("${SM}_fastqc.zip") into SM_fastqc_stats_set
         file("${SM}.picard.*") into SM_picard_stats_set

    """
        samtools stats --threads=${task.cpus} ${SM}.bam > ${SM}.samtools.txt
        bamtools -in ${SM}.bam > ${SM}.bamtools.txt
        fastqc --threads ${task.cpus} ${SM}.bam
        picard CollectAlignmentSummaryMetrics R=${reference_handle} I=${SM}.bam O=${SM}.picard.alignment_metrics.txt
        picard CollectInsertSizeMetrics I=${SM}.bam O=${SM}.picard.insert_metrics.txt H=${SM}.picard.insert_histogram.txt
    """
}

process bam_publish {

    publishDir "${params.out}/BAM", mode: 'copy', pattern: '*.bam*'

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_publish
    output:
        set file("${SM}.bam"), file("${SM}.bam.bai")

    """
        echo "${SM} saved to publish folder you rockstar."
    """
}

process SM_idx_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_idxstats
    output:
        file("${SM}.bam_idxstats") into bam_idxstats_set
        file("${SM}.bam_idxstats") into bam_idxstats_multiqc

    """
        samtools idxstats ${SM}.bam | awk '{ print "${SM}\\t" \$0 }' > ${SM}.bam_idxstats
    """
}

process SM_combine_idx_stats {

    publishDir "${params.out}/alignment", mode: "copy"

    input:
        val bam_idxstats from bam_idxstats_set.toSortedList()

    output:
        file("isotype_bam_idxstats.tsv")

    """
        echo -e "SM\\treference\\treference_length\\tmapped_reads\\tunmapped_reads" > isotype_bam_idxstats.tsv
        cat ${bam_idxstats.join(" ")} >> isotype_bam_idxstats.tsv
    """
}

/*
    =================
    Isotype BAM stats
    =================
*/

process isotype_bam_stats {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_stats

    output:
        file 'bam_stat' into SM_bam_stat_files

    """
        samtools stats ${SM}.bam | \\
        grep ^SN | \\
        cut -f 2- | \\
        awk '{ print "${SM}\t" \$0 }' | \\
        sed 's/://g' > bam_stat
    """
}

process combine_isotype_bam_stats {

    publishDir "${params.out}/alignment", mode: "copy"

    input:
        val stat_files from SM_bam_stat_files.collect()

    output:
        file("isotype_bam_stats.tsv")

    """
        echo -e "fq_pair_id\\tvariable\\tvalue\\tcomment" > isotype_bam_stats.tsv
        cat ${stat_files.join(" ")} >> SM_stats.tsv
    """
}

/*
    ============
    Coverage BAM
    ============
*/
process coverage_SM {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_coverage

    output:
        val SM into isotype_coverage_sample
        file("${SM}.coverage.tsv") into isotype_coverage


    """
        bam coverage ${SM}.bam > ${SM}.coverage.tsv
    """
}

process coverage_SM_merge {

    publishDir "${params.out}/alignment", mode: 'copy'

    input:
        val sm_set from isotype_coverage.toSortedList()

    output:
        file("isotype_coverage.full.tsv") into mt_content
        file("isotype_coverage.tsv") into isotype_coverage_merged

    """
        echo -e 'bam\\tcontig\\tstart\\tend\\tproperty\\tvalue' > isotype_coverage.full.tsv
        cat ${sm_set.join(" ")} >> isotype_coverage.full.tsv

        # Generate condensed version
        cat <(echo -e 'strain\\tcoverage') <(cat isotype_coverage.full.tsv | grep 'genome' | grep 'depth_of_coverage' | cut -f 1,6) > isotype_coverage.tsv
    """
}

/*
    ==========
    MT content
    ==========
*/

process output_mt_content {

    executor 'local'

    publishDir "${params.out}/phenotype", mode: "copy"

    input:
        file("isotype_coverage.full.tsv") from mt_content

    output:
        file("MT_content.tsv")

    """
        cat <(echo -e 'isotype\\tmt_content') <(cat isotype_coverage.full.tsv | awk '/mt_nuclear_ratio/' | cut -f 1,6) > MT_content.tsv
    """
}

/*
    ======
    telseq
    ======
*/

process call_telseq {

    tag { SM }

    input:
        set val(SM), file("${SM}.bam"), file("${SM}.bam.bai") from bam_telseq
    output:
        file("telseq_out.txt") into telseq_results

    """
        telseq -z TTAGGC -H ${SM}.bam > telseq_out.txt
    """
}

process combine_telseq {

    executor 'local'

    publishDir "${params.out}/phenotype", mode: 'copy'

    input:
        file("ind_telseq?.txt") from telseq_results.toSortedList()

    output:
        file("telseq.tsv")

    '''
        telseq -h > telseq.tsv
        cat ind_telseq*.txt | egrep -v '\\[|BAMs' >> telseq.tsv
    '''
}

/*
    ========================
    Call Variants - BCFTools
    ========================    
*/

process split_genome_full {

    cpus 1
    memory '4 GB' 
    time '1h'

    output:
    file "scatter/*-scattered.intervals" into interval_ch

    script:
    
    """
    gatk --java-options "-Xmx4g -Xms4g" \\
     SplitIntervals \\
        -R ${reference_handle} \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \\
        --scatter-count 20 \\
        -ip 250 \\
        -O scatter
    """

  }

bam_snp_individual
    .spread(interval_ch)
    .set {bam_for_calling}

process gatk_call {

    memory '64 GB'

    tag { "${SM} - ${region}" }
    cpus params.cores

    input:
      set val(SM), file(smbam), file(smbambai), file(region) from bam_for_calling

    output:
      set val(SM), file("${SM}.${int_tag}.g.vcf"), file("${SM}.${int_tag}.g.vcf.idx") into individual_sites
      set val(SM), file(region), file("*sample_map.txt") into sample_map

    script:
      int_tag = region.baseName

    """
      gatk HaplotypeCaller --java-options "-Xmx32g -Xms32g" \\
          -R ${reference_handle} \\
          -I ${smbam} \\
          --emit-ref-confidence GVCF \\
          --genotyping-mode DISCOVERY \\
          --max-genotype-count 3000 \\
          --max-alternate-alleles 100 \\
          --annotation DepthPerAlleleBySample \\
          --annotation Coverage \\
          --annotation GenotypeSummaries \\
          --annotation TandemRepeat \\
          --annotation StrandBiasBySample \\
          --annotation ChromosomeCounts \\
          --annotation AS_QualByDepth \\
          --annotation AS_StrandOddsRatio \\
          --annotation AS_MappingQualityRankSumTest \\
          --annotation DepthPerSampleHC \\
          --annotation-group StandardAnnotation \\
          --annotation-group AS_StandardAnnotation \\
          --annotation-group StandardHCAnnotation \\
          -L ${region} \\
          -O ${SM}.${int_tag}.g.vcf   

      gatk IndexFeatureFile \\
          -F ${SM}.${int_tag}.g.vcf 

      intname=`ls *.intervals`
      gvcf=`find . -name '*.g.vcf' | cut -f2 -d "/"`

      echo "${SM} \${gvcf}" | awk -v OFS="\t" '\$1=\$1' > ${SM}.\$intname.sample_map.txt
      echo "samplemaps" > place_holder.txt
    """
}

sample_map
  .groupTuple()
  .into{merged_sample_maps;
       print_merged_sample_maps}

individual_sites
  .groupTuple()
  .join(merged_sample_maps)
  .into{merge_sm_int_gvcfs;
       print_sm_int_gvcfs}

process merge_sm_gvcfs {

    memory '64 GB'

    tag { "${SM}" }

    input:
      set val(SM), file(gvcf), file(index), file(intervals), file(sample_map) from merge_sm_int_gvcfs

    output:
      file("${SM}_merged-intervals.g.vcf") into merged_sample_gvcf
      file("${SM}_merged-intervals.g.vcf.idx") into merged_sample_gvcf_index

    """
      cat *sample_map.txt > ${SM}_merged_sample-map.txt
      cut -f2 ${SM}_merged_sample-map.txt > ${SM}_intervals.list

      gatk MergeVcfs --java-options "-Xmx16g -Xms16g" \\
        --INPUT ${SM}_intervals.list \\
        --OUTPUT ${SM}_merged-intervals.g.vcf
    """
}

merged_sample_gvcf
  .toSortedList()
  .into{merged_sample_gvcf_to_sample_map;
        merged_sample_gvcf_to_db;
        merged_sample_gvcf_other;}


merged_sample_gvcf_index
  .toSortedList()
  .into{merged_sample_gvcf_index_to_sample_map;
        merged_sample_gvcf_index_to_db;
        merged_sample_gvcf_index_to_other;}



process cohort_to_sample_map {

    memory '64 GB'

    publishDir "${params.out}/Isotype_gVCF", mode: 'copy'

    executor 'local'

    input:
      file(gvcfs) from merged_sample_gvcf_to_sample_map
      file(indices) from merged_sample_gvcf_index_to_sample_map

    output:
      file("cohort.sample_map") into cohort_map

    """
      ls *_merged-intervals.g.vcf > cohort.sample_map_tmp
      awk -F"_" 'BEGIN{OFS="\t"}; \$1=\$1' OFS="\t" cohort.sample_map_tmp | \\
      awk '{print \$1, \$1"_"\$2}' OFS="\t" > cohort.sample_map
    """
}

 process cohort_to_db {

    memory '64 GB'

    tag { "${chr}" }

    cpus 12

    publishDir "${params.out}/gVCF_db/${chr}", mode: 'copy'

    input:
      file(gvcfs) from merged_sample_gvcf_to_db
      file(indices) from merged_sample_gvcf_index_to_db
      file(samplemap) from cohort_map
      each chr from contigs

    output:
      set val(chr), file("${chr}_database.tar") into chromosomal_db

    """
      gatk --java-options "-Xmx32g -Xms32g" \\
       GenomicsDBImport \\
       --genomicsdb-workspace-path ${chr}_database \\
       --batch-size 16 \\
       -L ${chr} \\
       --sample-name-map ${samplemap} \\
       --reader-threads ${task.cpus} \\
       --genomicsdb-vcf-buffer-size 40000

      tar -cf ${chr}_database.tar ${chr}_database
    """
}

process genotype_cohort_gvcf_db {

    memory '64 GB'

    tag { "${chr}" }

    cpus 12

    publishDir "${params.out}/cohort_VCFs/${chr}", mode: 'copy'

    input:
      set val(chr), file(chr_db) from chromosomal_db

    output:
      set val(chr), file("${chr}_cohort.vcf"), file("${chr}_cohort.vcf.idx") into cohort_vcf

    """
      tar -xf ${chr_db}
      WORKSPACE=\$( basename ${chr_db} .tar)

      gatk --java-options "-Xmx48g -Xms48g" \\
        GenotypeGVCFs \\
        -R ${reference_handle} \\
        -V gendb://\$WORKSPACE \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        -G StandardHCAnnotation \\
        -L ${chr} \\
        --use-new-qual-calculator \\
        -O ${chr}_cohort.vcf
    """
}

process annotate_vcf {

    cpus 12

    tag { chrom }

    input:
        set val(chrom), file(vcf), file(vcfindex) from cohort_vcf

    output:
        set file("${chrom}.annotated.vcf.gz"), file("${chrom}.annotated.vcf.gz.tbi") into annotated_vcf
        file("snpeff_out.csv") into snpeff_multiqc


    """
      bcftools view -Oz -o ${vcf}.gz ${vcf}
      tabix -p vcf ${vcf}.gz

      # bcftools csq
      bcftools view --threads=${task.cpus-1} --regions ${chrom} ${vcf}.gz | \\
      snpEff eff -v chicken_new \\
                 -csvStats snpeff_out.csv \\
                 -no-downstream \\
                 -no-intergenic \\
                 -no-upstream \\
                 -nodownload \\ | \\
      bcftools view --threads=${task.cpus} -O z > ${chrom}.annotated.vcf.gz
      
      tabix -p vcf ${chrom}.annotated.vcf.gz

    """

}


contig_raw_vcf = CONTIG_LIST*.concat(".annotated.vcf.gz")

process concatenate_annotated_vcf {

    memory '64 GB'

    cpus 12

    publishDir "${params.out}/variation", mode: 'copy' 

    input:
        file(merge_vcf) from annotated_vcf.collect()

    output:
        set file("CH.annotated.vcf.gz"), file("CH.annotated.vcf.gz.tbi") into annotated_concatenated_vcf
        set val("soft"), file("CH.annotated.vcf.gz"), file("CH.annotated.vcf.gz.tbi") into annotated_sample_summary
        file("CH.annotated.stats.txt") into annotated_stats

    """
      bcftools concat --threads ${task.cpus-1} -O z ${contig_raw_vcf.join(" ")} > CH.annotated.vcf.gz
      tabix -p vcf CH.annotated.vcf.gz
      bcftools stats --verbose CH.annotated.vcf.gz > CH.annotated.stats.txt
    """
}

annotated_concatenated_vcf
  .into{ann_vcf_to_soft_filter;
        ann_vcf_to_hard_filter;
        ann_vcf_to_vqsr;
        ann_vcf_to_recal;
        ann_vcf_to_cnn_train;
        ann_vcf_to_cnn_apply;}


/*
================================================
~ > *                                      * < ~
~ ~ > *                                  * < ~ ~
~ ~ ~ > *   Apply Soft Filters on VCF  * < ~ ~ ~
~ ~ > *                                  * < ~ ~
~ > *                                      * < ~
================================================
*/

/*
============================================
~ ~ ~ > *   Split Indels and SNVs  * < ~ ~ ~
============================================
*/

process split_snv_indel {

    memory '64 GB'

    cpus 8

    publishDir "${params.out}/variation", mode: 'copy' 

    input:
      set file(merge_ann_vcf), file(merge_ann_index) from ann_vcf_to_soft_filter

    output:
      set file("merged_annotated_SNV.vcf"), file("merged_annotated_SNV.vcf.idx") into split_snv_vcf
      set file("merged_annotated_INDEL.vcf"), file("merged_annotated_INDEL.vcf.idx") into split_indel_vcf


  """
    gatk --java-options "-Xmx48g -Xms48g" \\
        SelectVariants \\
        -R ${reference_handle} \\
        -V ${merge_ann_vcf} \\
        --select-type-to-include SNP \\
        -O merged_annotated_SNV.vcf


    gatk --java-options "-Xmx48g -Xms48g" \\
        SelectVariants \\
        -R ${reference_handle} \\
        -V ${merge_ann_vcf} \\
        --select-type-to-include INDEL \\
        -O merged_annotated_INDEL.vcf
  """

}


/*
================================================
~ ~ ~ > *   Apply Indels Soft Filters  * < ~ ~ ~
================================================
*/

process apply_soft_filters_indel {

    publishDir "${params.out}/variation", mode: 'copy'

    cpus params.cores

    input:
      set file(indel_vcf), file(indel_index) from split_indel_vcf

    output:
      set file("soft_filtered_indels.vcf"), file("soft_filtered_indels.vcf.idx") into soft_filter_indels


    """
      bcftools norm --threads ${task.cpus} -m -any -Ov -o norm_indels.vcf ${indel_vcf} 

      gatk --java-options "-Xmx48g -Xms48g" \\
          VariantFiltration \\
          -R ${reference_handle} \\
          --variant norm_indels.vcf \\
          --genotype-filter-expression "DP < ${params.min_depth}" \\
          --genotype-filter-name "depth" \\
          --filter-expression "QD < ${params.quality_by_depth}" \\
          --filter-name "qd" \\
          --filter-expression "SOR > ${params.strand_odds_ratio}" \\
          --filter-name "sor" \\
          --filter-expression "QUAL < ${params.qual}" \\
          --filter-name "quality" \\
          --filter-expression "ReadPosRankSum < ${params.readbias}" \\
          --filter-name "readend" \\
          --filter-expression "FS > ${params.fisherstrand}" \\
          --filter-name "fisherstrand" \\
          -O soft_filtered_indels.vcf
    """
}

/*
=============================================
~ ~ ~ > *   Apply SNV Soft Filters  * < ~ ~ ~
=============================================
*/

process apply_soft_filters_snv {

    publishDir "${params.out}/variation", mode: 'copy'

    cpus params.cores

    input:
      set file(snp_vcf), file(snp_index) from split_snv_vcf

    output:
      set file("soft_filtered_snps.vcf"), file("soft_filtered_snps.vcf.idx") into soft_filter_snvs


    """
      gatk --java-options "-Xmx48g -Xms48g" \\
          VariantFiltration \\
          -R ${reference_handle} \\
          --variant ${snp_vcf} \\
          --genotype-filter-expression "DP < ${params.min_depth}" \\
          --genotype-filter-name "depth" \\
          --filter-expression "QUAL < ${params.qual}" \\
          --filter-name "quality" \\
          --filter-expression "ReadPosRankSum < ${params.readbias}" \\
          --filter-name "readend" \\
          --filter-expression "FS > ${params.fisherstrand}" \\
          --filter-name "fisherstrand" \\
          --filter-expression "QD < ${params.quality_by_depth}" \\
          --filter-name "qd" \\
          --filter-expression "SOR > ${params.strand_odds_ratio}" \\
          --filter-name "sor" \\
          -O soft_filtered_snps.vcf
    """
}


/*
==============================================
~ ~ ~ > *   Combine SNVs and Indels  * < ~ ~ ~
==============================================
*/

process combine_soft_filter_vcfs {

    publishDir "${params.out}/variation", mode: 'copy'

    cpus params.cores

    input:
      set file("soft_filtered_snps.vcf"), file("soft_filtered_snps.vcf.isx") from soft_filter_snvs
      set file("soft_filtered_indels.vcf"), file("soft_filtered_indels.vcf.idx") from soft_filter_indels

    output:
      set file("CH.soft-filter.vcf.gz"), file("CH.soft-filter.vcf.gz.tbi") into soft_filtered_cohort_vcf
      file("CH.soft-filter.stats.txt") into soft_filtered_stats_to_mqc


    """
      bcftools view -Oz -o soft_filtered_indels.vcf.gz soft_filtered_indels.vcf
      tabix -p vcf soft_filtered_indels.vcf.gz

      bcftools view -Oz -o soft_filtered_snps.vcf.gz soft_filtered_snps.vcf
      tabix -p vcf soft_filtered_snps.vcf.gz

      bcftools concat \\
      --threads ${task.cpus} \\
      --allow-overlaps \\
      soft_filtered_indels.vcf.gz \\
      soft_filtered_snps.vcf.gz | \\
      bcftools filter -Oz --threads ${task.cpus} --mode + --soft-filter high_missing --include "F_MISSING<=${params.missing}" > CH.soft-filter.vcf.gz

      tabix -p vcf CH.soft-filter.vcf.gz
      bcftools stats --verbose CH.soft-filter.vcf.gz > CH.soft-filter.stats.txt
    """
}

soft_filtered_cohort_vcf
  .into{soft_vcf_to_strain_list;
        soft_vcf_to_split_by_strain;
        soft_vcf_to_other}

/*
==========================================
~ ~ ~ > *   Split VCF by sample  * < ~ ~ ~
==========================================
*/

process generate_strain_list {

    executor 'local'

    input:
        set file(softvcf), file(softvcf_index) from soft_vcf_to_strain_list

    output:
        file('isotype_list.tsv') into isotype_list

    """
        bcftools query -l ${softvcf} > isotype_list.tsv
    """

}

isotype_list
  .splitText() { it.strip() } 
  .combine(soft_vcf_to_split_by_strain)
  .into{isotype_set_vcf; 
        isotype_set_tsv}


/*
===========================================
~ ~ ~ > *   Apply AD soft filter  * < ~ ~ ~
===========================================
*/

process apply_allele_depth_filter {

    publishDir "${params.out}/isotype/vcf", mode: 'copy'

    tag { isotype }

    memory '64 GB'

    input:
        set val(isotype), file(softvcf), file(softvcf_vcf) from isotype_set_vcf

    output:
        file("${isotype}.AD-filter.vcf.gz") into isotype_AD_soft_vcf
        file("${isotype}.AD-filter.vcf.gz.tbi") into isotype_AD_soft_vcf_index

    """
    bcftools view --samples ${isotype} ${softvcf} |\\
    bcftools filter -Ov --mode + --soft-filter dv_dp --include "((FORMAT/AD[*:1])/(FORMAT/DP) >= 0.5) || (FORMAT/GT == '0/0') || (TYPE == 'REF')" -Ov -o ${isotype}_temp.vcf

    gatk --java-options "-Xmx4g -Xms4g" \\
         VariantFiltration \\
         -R ${reference_handle} \\
         --variant ${isotype}_temp.vcf \\
         --genotype-filter-expression "FILTER != 'PASS'" \\
         --genotype-filter-name "dv_dp" \\
         -O temp_soft_filtered.vcf

    awk 'BEGIN{FS=OFS="\\t"} {gsub("dv_dp",\$7,\$10)} 1' temp_soft_filtered.vcf | \\
    bcftools view -Oz -o ${isotype}.AD-filter.vcf.gz

    tabix -p vcf ${isotype}.AD-filter.vcf.gz
    """

}

contigs_soft = Channel.from(CONTIG_LIST)

contigs_soft
  .into{contigs_vcf;
        contigs_index;
        }

contigs_vcf
  .spread(isotype_AD_soft_vcf)
  .groupTuple()
  .set{ to_merge_soft_sm_ad }

contigs_index
  .spread(isotype_AD_soft_vcf_index)
  .groupTuple()
  .into{ to_merge_soft_sm_ad_index;
        print_merged_index }

/*
==========================================
~ ~ ~ > *   Combine Cohort VCFs  * < ~ ~ ~
==========================================
*/

process merge_sm_soft_vcfs {

    tag { chrom }

    publishDir "${params.out}/variation/", mode: 'copy'

    memory '64 GB'
    cpus 4

    input:
        set val(chrom), file(softvcf) from to_merge_soft_sm_ad
        set val(chrom), file(softvcf_vcf) from to_merge_soft_sm_ad_index

    output:
        file("CH.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz") into cohort_soft_filter_vcf
        file("CH.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into cohort_soft_filter_vcf_index

    """
        bcftools merge \\
        -m none \\
        -r ${chrom} \\
        -Oz -o CH.${chrom}.COMPLETE.vcf.gz \\
        ${softvcf}

        bcftools query -l CH.${chrom}.COMPLETE.vcf.gz | sort > samples.txt

        bcftools view -S samples.txt CH.${chrom}.COMPLETE.vcf.gz -Oz -o CH.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz

        tabix -p vcf CH.${chrom}.COMPLETE-SOFT-FILTER.vcf.gz
    """
}

cohort_soft_filter_vcf
  .toSortedList()
  .set{ to_concat_soft_filter_vcf }

cohort_soft_filter_vcf_index
  .toSortedList()
  .set{ to_concat_soft_filter_vcf_index }

/*
====================================================
~ ~ ~ > *   Concatenate Cohort CHROM VCFs  * < ~ ~ ~
====================================================
*/

process concat_cohort_soft_vcfs {

    publishDir "${params.out}/variation/", mode: 'copy'

    memory '64 GB'
    cpus 20

    input:
        file(chromvcf) from to_concat_soft_filter_vcf
        file(chromvcf_index) from to_concat_soft_filter_vcf_index

    output:
        set file("CH.COMPLETE-SOFT-FILTER.vcf.gz"), file("CH.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into cohort_complete_soft_filter_vcf
        set val("soft"), file("CH.COMPLETE-SOFT-FILTER.vcf.gz"), file("CH.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into soft_vcf_summary
        set val("soft"), file("CH.COMPLETE-SOFT-FILTER.vcf.gz"), file("CH.COMPLETE-SOFT-FILTER.vcf.gz.tbi") into soft_sample_summary
        file("CH.COMPLETE-SOFT-FILTER.stats.txt") into complete_soft_filter_vcf_stats

    """
      bcftools concat \\
      --threads ${task.cpus} \\
      *.vcf.gz |\\
      bcftools filter -Oz --mode=+x --soft-filter="high_missing" --include 'F_MISSING  <= ${params.missing}' | \\
      bcftools filter -Oz --mode=+x --soft-filter="high_heterozygosity" --include '(COUNT(GT="het")/N_SAMPLES <= 0.10)' |\\
      bcftools view -Oz -o CH.COMPLETE-SOFT-FILTER.vcf.gz

      tabix -p vcf CH.COMPLETE-SOFT-FILTER.vcf.gz
      bcftools stats --verbose CH.COMPLETE-SOFT-FILTER.vcf.gz > CH.COMPLETE-SOFT-FILTER.stats.txt
    """
}

cohort_complete_soft_filter_vcf
  .into{soft_filtered_vcf_to_hard;
        soft_filtered_vcf_gtcheck;
        soft_filter_vcf_strain;
        soft_filter_vcf_isotype_list;
        soft_filter_vcf_mod_tracks;
        soft_filter_vcf_tsv
      }

/*
================================================
~ > *                                      * < ~
~ ~ > *                                  * < ~ ~
~ ~ ~ > *   Apply Hard Filters on VCF  * < ~ ~ ~
~ ~ > *                                  * < ~ ~
~ > *                                      * < ~
================================================
*/

process generate_hard_vcf {

    cpus 20

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        set file(softvcf), file(softvcfindex) from soft_filtered_vcf_to_hard

    output:
        set file("CH.HARD-FILTERED.vcf.gz"), file("CH.HARD-FILTERED.vcf.gz.tbi") into hard_vcf
        set val("hard"), file("CH.HARD-FILTERED.vcf.gz"), file("CH.HARD-FILTERED.vcf.gz.tbi") into hard_vcf_summary
        set val("hard"), file("CH.HARD-FILTERED.vcf.gz"), file("CH.HARD-FILTERED.vcf.gz.tbi") into hard_sample_summary
        file("CH.HARD-FILTERED.stats.txt") into hard_filter_stats


    """
        # Generate hard-filtered VCF
        function generate_hard_filter {
            bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${softvcf} |\\
            bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT != "PASS"' |\\
            bcftools filter -O u --include 'F_MISSING  <= ${params.missing}' |\\
            bcftools filter -O u --include '(COUNT(GT="het")/N_SAMPLES <= 0.10)' |\\
            bcftools view -O v --min-af 0.0000000000001 --max-af 0.999999999999 |\\
            vcffixup - | \\
            bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz
        }

        export -f generate_hard_filter

        parallel --verbose generate_hard_filter {} ::: I II III IV V X MtDNA

        bcftools concat -O z I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz > CH.HARD-FILTERED.vcf.gz
        
        tabix -p vcf CH.HARD-FILTERED.vcf.gz

        bcftools stats --verbose CH.HARD-FILTERED.vcf.gz > CH.HARD-FILTERED.stats.txt

        # Remove extra files
        rm I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz
    """
}

process multiqc_report {

    publishDir "${params.out}/report", mode: 'copy'

    input:
        file("CH.HARD-FILTERED.stats.txt") from hard_filter_stats
        file("CH.COMPLETE-SOFT-FILTER.stats.txt") from complete_soft_filter_vcf_stats
        file(bam_tools) from SM_bamtools_stats_set.toSortedList()
        file(samtools_stats) from SM_samtools_stats_set.toSortedList()
        file(duplicates) from duplicates_set.toSortedList()
        file(fastqc) from SM_fastqc_stats_set.toSortedList()
        file("bam*.idxstats") from bam_idxstats_multiqc.toSortedList()
        file("picard*.stats.txt") from SM_picard_stats_set.collect()
        file("snpeff_out.csv") from snpeff_multiqc

    output:
        file("multiqc_data/*.json") into multiqc_json_files
        file("multiqc.html")

    """
        multiqc -k json --filename multiqc.html .
    """

}



workflow.onComplete {

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]CH

    """

    println summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }

    // mail summary
    if (params.email) {
        ['mail', '-s', 'wi-nf', params.email].execute() << summary
    }


}