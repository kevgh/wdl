task GetBwaVersion {
  command {
    /usr/gitc/bwa 2>&1 | grep -e '^Version' | sed 's/Version: //'
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "1 GB"
  }
  output {
    String version = read_string(stdout())
  }
}

task SamToFastqAndBwaMem {
  File input_bam
  String bwa_commandline
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # only relevant for hg38 and on this file lists the reference sequences that are "alternative". 
  # references without this file should present an empty file
  File ref_alt

  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Int disk_size

  command <<<
    set -o pipefail
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    # if ref_alt has data in it, then we use hg38 reference
    if [ -s ${ref_alt} ]; then
      java -Xmx3000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=${input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        CLIPPING_ATTRIBUTE=XT \
        CLIPPING_ACTION=2 \
        NON_PF=true |
      /usr/gitc/${bwa_commandline} /dev/stdin - |
      /usr/gitc/k8 /usr/gitc/bwa-postalt.js ${ref_alt} |
      /usr/gitc/samtools view -1 - > ${output_bam_basename}.bam
    # else ref_alt is empty and we use hg19 reference 
    else
      java -Xmx3000m -jar /usr/gitc/picard.jar \
        SamToFastq \
        INPUT=${input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        CLIPPING_ATTRIBUTE=XT \
        CLIPPING_ACTION=2 \
        NON_PF=true |
      /usr/gitc/${bwa_commandline} /dev/stdin - |
      /usr/gitc/samtools view -1 - > ${output_bam_basename}.bam
    fi
  >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "14 GB"
    cpu: "16"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }
}

task MergeBamAlignment {
  File unmapped_bam
  String bwa_commandline
  String bwa_version
  File aligned_bam
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict  
  Int disk_size

  command {

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    java -Xmx3000m -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ALIGNED_BAM=${aligned_bam} \
      UNMAPPED_BAM=${unmapped_bam} \
      OUTPUT=${output_bam_basename}.bam \
      REFERENCE_SEQUENCE=${ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="${bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "3500 MB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
  }  
}

task SortAndFixTags {
  File input_bam
  String output_bam_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size
  
  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
    SortSam \
    INPUT=${input_bam} \
    OUTPUT=/dev/stdout \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=false \
    CREATE_MD5_FILE=false | java -Xmx500m -jar /usr/gitc/picard.jar \
    SetNmAndUqTags \
    INPUT=/dev/stdin \
    OUTPUT=${output_bam_basename}.bam \
    CREATE_INDEX=true \
    CREATE_MD5_FILE=true \
    REFERENCE_SEQUENCE=${ref_fasta}

  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "5000 MB"
    preemptible: 3
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

task MarkDuplicates {
  Array[File] input_bams
  String output_bam_basename
  String metrics_filename
  Int disk_size

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This workes because the output of BWA is query-grouped, and thus so is the output of MergeBamAlignment. 
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=${sep=' INPUT=' input_bams} \
      OUTPUT=${output_bam_basename}.bam \
      METRICS_FILE=${metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname"
      CREATE_MD5_FILE=true
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File duplicate_metrics = "${metrics_filename}"
  }
}

task CreateSequenceGroupingTSV {
  File ref_dict

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.  It outputs to stdout
  # where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]

    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]

    print tsv_string
    CODE
  >>>
  runtime {
    docker: "python:2.7"
    memory: "2 GB"
    preemptible: 3
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv(stdout())
  }
}

task BaseRecalibrator {
  File input_bam
  File input_bam_index
  String recalibration_report_filename
  Array[String] sequence_group_interval
  File dbSNP_vcf
  File dbSNP_vcf_index
  File known_snps_sites_vcf
  File known_snps_sites_vcf_index
  File known_indels_sites_vcf
  File known_indels_sites_vcf_index
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Dsamjdk.use_async_io=false -Xmx4000m \
      -jar /usr/gitc/GATK4.jar \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${recalibration_report_filename} \
      -knownSites ${dbSNP_vcf} \
      -knownSites ${known_snps_sites_vcf} \
      -knownSites ${known_indels_sites_vcf} \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "6 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File recalibration_report = "${recalibration_report_filename}"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}

task ApplyBQSR {
  File input_bam
  File input_bam_index
  String output_bam_basename
  File recalibration_report
  Array[String] sequence_group_interval
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  command {
    java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx3000m \
      -jar /usr/gitc/GATK4.jar \
      ApplyBQSR \
      --createOutputBamMD5 \
      --addOutputSAMProgramRecord \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      -SQQ 10 -SQQ 20 -SQQ 30 -SQQ 40 \
      --emit_original_quals \
      -L ${sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "3500 MB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File recalibrated_bam = "${output_bam_basename}.bam"
    File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
    #this output is only for GOTC STAGING to give some GC statistics to the GATK4 team
    #File gc_logs = "gc_log.log"
  }
}

task GatherBqsrReports {
  Array[File] input_bqsr_reports
  String output_report_filename
  Int disk_size

  command {
    java -Xmx3000m -jar /usr/gitc/GATK4.jar \
      GatherBQSRReports \
      -I ${sep=' -I ' input_bqsr_reports} \
      -O ${output_report_filename}
    }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "3500 MB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File output_bqsr_report = "${output_report_filename}"
  }
}

task GatherBamFiles {
  Array[File] input_bams
  File input_unmapped_reads_bam
  String output_bam_basename
  Int disk_size

  command {
    java -Xmx2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=${sep=' INPUT=' input_bams} \
      INPUT=${input_unmapped_reads_bam} \
      OUTPUT=${output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true

    }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    File output_bam_md5 = "${output_bam_basename}.bam.md5"
  }
}

task ScatterIntervalList {
  File interval_list
  Int scatter_count
  Int break_bands_at_multiples_of

  command <<<
    mkdir out
    java -Xmx1g -jar /usr/gitc/picard.jar \
    IntervalListTools \
    SCATTER_COUNT=${scatter_count} \
    SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    UNIQUE=true \
    SORT=true \
    BREAK_BANDS_AT_MULTIPLES_OF=${break_bands_at_multiples_of} \
    INPUT=${interval_list} \
    OUTPUT=out
    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    for i, interval in enumerate(glob.glob("out/*/*.interval_list")):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i) + filename)
      os.rename(interval, newName)
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
  }

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "2 GB"
    preemptible: 3
  }
}

task HaplotypeCaller {
  File input_bam
  File input_bam_index
  Float contamination
  File interval_list
  String gvcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Int disk_size

  # tried to find lowest memory variable where it would still work, might change once tested on JES
  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx8000m \
      -jar /usr/gitc/GenomeAnalysisTK.jar \
      -T HaplotypeCaller \
      -R ${ref_fasta} \
      -o ${gvcf_basename}.vcf.gz \
      -I ${input_bam} \
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ${contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
  output {
    File output_gvcf = "${gvcf_basename}.vcf.gz"
    File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  }
}

task GatherVCFs {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String output_vcf_name
  Int disk_size

  # using MergeVcfs instead of GatherVcfs so we can create indicies
  # WARNING	2015-10-28 15:01:48	GatherVcfs	Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xmx2g -jar /usr/gitc/picard.jar \
    MergeVcfs \
    INPUT=${sep=' INPUT=' input_vcfs} \
    OUTPUT=${output_vcf_name}
  }
  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }
}

task ConvertToCram {
  File input_bam
  File ref_fasta
  File ref_fasta_index
  String output_basename
  Int disk_size

  # Note that we are not activating pre-emptible instances for this step yet,
  #  but we should if it ends up being fairly quick
  command <<<
      /usr/gitc/samtools view -C -T ${ref_fasta} ${input_bam} > ${output_basename}.cram
      /usr/gitc/samtools index ${output_basename}.cram
      mv ${output_basename}.cram.crai ${output_basename}.crai
  >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.0.0"
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram = "${output_basename}.cram"
    File output_cram_index = "${output_basename}.crai"
  }
}

workflow PairedEndSingleSampleWorkflow {
  String sample_name
  String final_gvcf_name
  Array[File] flowcell_unmapped_bams
  String unmapped_bam_suffix
  File calling_interval_list
  File dbSNP_vcf
  File dbSNP_vcf_index
  File known_snps_sites_vcf
  File known_snps_sites_vcf_index
  File known_indels_sites_vcf
  File known_indels_sites_vcf_index
  File haplotype_database_file
  Int  haplotype_scatter_count
  Int  break_bands_at_multiples_of

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt

  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  Int flowcell_small_disk
  Int flowcell_medium_disk
  Int agg_small_disk
  Int agg_medium_disk
  Int agg_large_disk
  Int agg_extra_large_disk

  String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 $bash_ref_fasta"

  String recalibrated_bam_basename = sample_name + ".aligned.duplicates_marked.recalibrated"

  call GetBwaVersion
  
  # PRE-PROCESSING PHASE 1
  # Align flowcell-level unmapped input bams in parallel
  
  scatter (unmapped_bam in flowcell_unmapped_bams) {
  
    # Because of a wdl/cromwell bug the line below is not currently valid so we have to sub(sub()) in each task
    # String base_name = sub(sub(unmapped_bam, "gs://.*/", ""), ".unmapped.bam$", "")

    String sub_strip_path = "gs://.*/"
    String sub_strip_unmapped = unmapped_bam_suffix + "$"
    
    call SamToFastqAndBwaMem {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".unmerged",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        disk_size = flowcell_medium_disk
    }

    call MergeBamAlignment {
      input:      
        unmapped_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        bwa_version = GetBwaVersion.version,
        aligned_bam = SamToFastqAndBwaMem.output_bam,
        output_bam_basename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        disk_size = flowcell_medium_disk
    }   

    call SortAndFixTags as SortAndFixReadGroupBam {
      input:
      input_bam = MergeBamAlignment.output_bam,
      output_bam_basename = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "") + ".sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,  
      disk_size = flowcell_medium_disk
    }
  }
 
  # PRE-PROCESSING PHASE 2
  # Aggregate aligned+merged flowcell bams, mark duplicates, and recalibrate qualities
  
  call MarkDuplicates {
    input:
      input_bams = MergeBamAlignment.output_bam,
      output_bam_basename = sample_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_name + ".duplicate_metrics",
      disk_size = agg_extra_large_disk
  }

  call SortAndFixTags as SortAndFixSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = sample_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,  
      disk_size = agg_extra_large_disk
  }

  call CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
  
    call BaseRecalibrator {
      input:
        input_bam = SortAndFixSampleBam.output_bam,
        input_bam_index = SortAndFixSampleBam.output_bam_index,
        recalibration_report_filename = sample_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_snps_sites_vcf = known_snps_sites_vcf,
        known_snps_sites_vcf_index = known_snps_sites_vcf_index,
        known_indels_sites_vcf = known_indels_sites_vcf,
        known_indels_sites_vcf_index = known_indels_sites_vcf_index,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = agg_medium_disk
    }
    
    call ApplyBQSR {
      input:
        input_bam = SortAndFixSampleBam.output_bam,
        input_bam_index = SortAndFixSampleBam.output_bam_index,
        output_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        disk_size = agg_medium_disk
    }
  }
  
  call GatherBqsrReports{
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = sample_name + ".recal_data.csv",
      disk_size = agg_small_disk
  }
  
  # This ApplyBQSR uses "-L unmapped" as its interval to grab all paired unmapped reads so they can be added to the final bam
  Array[String] unmapped_group_interval = ["unmapped"]
  
  call ApplyBQSR as ApplyBQSRToUnmappedReads {
    input:
      input_bam = SortAndFixSampleBam.output_bam,
      input_bam_index = SortAndFixSampleBam.output_bam_index,
      output_bam_basename = recalibrated_bam_basename,
      recalibration_report = GatherBqsrReports.output_bqsr_report,
      sequence_group_interval = unmapped_group_interval,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = agg_medium_disk
  }
  
  ## TODO when capability of adding elements to arrays, can just have one array as an input and add the output of the above task to the scattered printreads bams
  call GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      input_unmapped_reads_bam = ApplyBQSRToUnmappedReads.recalibrated_bam,
      output_bam_basename = sample_name,
      disk_size = agg_extra_large_disk
  }
  
  # Make a CRAM file for archival
  call ConvertToCram {
    input:
      input_bam = GatherBamFiles.output_bam,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = sample_name,
      disk_size = agg_large_disk
  }
  
  # VARIANT CALLING PHASE 1 
  # Phase 2 happens in a separate WDL
  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  
  call ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }
  
  scatter (subInterval in ScatterIntervalList.out) {
     call HaplotypeCaller {
       input:
         input_bam = GatherBamFiles.output_bam,
         input_bam_index = GatherBamFiles.output_bam_index,
         contamination = PreBqsrCheckContamination.contamination,
         interval_list = subInterval,
         gvcf_basename = sample_name,
         ref_dict = ref_dict,
         ref_fasta = ref_fasta,
         ref_fasta_index = ref_fasta_index,
         disk_size = agg_medium_disk
     }
  }
  
  call GatherVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_vcf_name = final_gvcf_name,
      disk_size = agg_medium_disk
  }

  output {
    MarkDuplicates.duplicate_metrics
    GatherBqsrReports.*
    ConvertToCram.*
    GatherVCFs.*
  }
}