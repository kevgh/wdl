task TileDBHistogram {
  Int size_per_column_partition
  Array[File] gvcf_list
  Array[File] gvcf_index_list
  File vid_mapping_file
  File ref_fasta
  File ref_fasta_index

  command <<<
    python <<CODE
    from collections import OrderedDict
    import json

    w = open("/usr/gitc/callset.json", "w")
    array_of_vcfs = ["${sep='","' gvcf_list}"]

    file_division = [array_of_vcfs]

    vcf_counter = 0
    callsets = OrderedDict()
    for vcf in array_of_vcfs:
      vcf_data = OrderedDict()
      vcf_data["row_idx"] = vcf_counter
      vcf_data["idx_in_file"] = 0
      vcf_data["filename"] = vcf
      sample_name = vcf.split("/")[-1].split(".")[0]
      callsets[sample_name] = vcf_data
      vcf_counter += 1

    final_json = OrderedDict()
    final_json["callsets"] = callsets
    final_json["file_division"] = file_division

    w.write(json.JSONEncoder().encode(final_json))
    w.close()

    w = open("/usr/gitc/loader.json", "w")

    column_partitions = [ {'begin': 0} ]

    final_json = OrderedDict()
    final_json["column_partitions"] = column_partitions
    final_json["callset_mapping_file"] = "/usr/gitc/callset.json"
    final_json["size_per_column_partition"] = ${size_per_column_partition}
    final_json["vid_mapping_file"] = "${vid_mapping_file}"
    final_json["max_histogram_range"] = 3999999999
    final_json["num_bins"] = 100000
    final_json["num_converter_processes"] = 1
    final_json["num_parallel_vcf_files"] = 16

    w.write(json.JSONEncoder().encode(final_json))
    w.close()
    CODE

    (/usr/bin/time -v /usr/gitc/vcf_histogram /usr/gitc/loader.json )

    >>>

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:tiledb"
    memory: "4 GB"
    cpu: "16"
    disks: "local-disk 1000 SSD"
  }
  output {
  }

}

task CreateTileDBPartitions {
  Int num_of_tiledb_partitions

  command <<<
    python <<CODE
    import math

    GENOME_SIZE = 3150000000

    partitions = ${num_of_partitions}
    base_split = int(math.floor(GENOME_SIZE / partitions))
    if partitions == 1 :
        partitions_json_string = "{'begin': 0}"
    else:
        partitions_json_string = "{'begin': 0, 'end': " + str(base_split) + "}\n"
        for partition in range(1, partitions - 1):
            current_partition_json = "{'begin': "
            current_partition_json += str(base_split * partition + 1) + ", 'end': "
            current_partition_json += str(base_split * (partition + 1)) + "}\n"
            partitions_json_string += current_partition_json

        partitions_json_string += "{'begin': " + str(base_split * (partitions - 1 ) + 1) + "}"
    print partitions_json_string
    CODE
    >>>

  runtime {
      docker: "python:2.7"
      memory: "2 GB"
    }

  output {
    Array[String] partition_json = read_lines(stdout())
    }
}

task TileDBCombineGVCF {
  String output_gvcf_filename
  String partition

  Array[File] gvcf_list
  Array[File] gvcf_index_list

  File ref_fasta
  File ref_fasta_index
  File vid_mapping_file
  File vcf_header

  Int size_per_column_partition

    command <<<
    python <<CODE
    from collections import OrderedDict
    import json

    w = open("/usr/gitc/callset.json", "w")
    array_of_vcfs = ["${sep='","' gvcf_list}"]

    vcf_counter = 0
    callsets = OrderedDict()
    for vcf in array_of_vcfs:
      vcf_data = OrderedDict()
      vcf_data["row_idx"] = vcf_counter
      vcf_data["idx_in_file"] = 0
      vcf_data["filename"] = vcf
      sample_name = vcf.split("/")[-1].split(".")[0]
      callsets[sample_name] = vcf_data
      vcf_counter += 1

    final_json = OrderedDict()
    final_json["callsets"] = callsets

    w.write(json.JSONEncoder().encode(final_json))
    w.close()

    w = open("/usr/gitc/loader.json", "w")

    column_partitions = [ ${partition} ]

    final_json = OrderedDict()
    final_json["column_partitions"] = column_partitions
    final_json["callset_mapping_file"] = "/usr/gitc/callset.json"
    final_json["size_per_column_partition"] = ${size_per_column_partition}
    final_json["treat_deletions_as_intervals"] = True
    final_json["vcf_header_filename"] = "${vcf_header}"
    final_json["vid_mapping_file"] = "${vid_mapping_file}"
    final_json["reference_genome"] = "${ref_fasta}"
    final_json["num_parallel_vcf_files"] = 1
    final_json["do_ping_pong_buffering"] = True
    final_json["offload_vcf_output_processing"] = True
    final_json["discard_vcf_index"] = True
    final_json["produce_combined_vcf"] = True
    # zipped output
    final_json["vcf_output_format"] = "z"
    final_json["vcf_output_filename"] = "${output_gvcf_filename}"

    w.write(json.JSONEncoder().encode(final_json))
    w.close()
    CODE

    (/usr/bin/time -v /usr/gitc/vcf2tiledb /usr/gitc/loader.json )

    >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:tiledb"
    memory: "4 GB"
    cpu: "2"
    disks: "local-disk 1000 SSD"
  }
  output {
    File output_gvcf = "${output_gvcf_filename}"
  }
}

task TabixAndGenotypeGVCF {
  String unfiltered_vcf_filename
  String unfiltered_vcf_index_filename

  File gvcf

  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbsnp_vcf
  File dbsnp_vcf_index

  command <<<
    /usr/gitc/tabix ${gvcf}

    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
      -jar /usr/gitc/GenomeAnalysisTK.jar \
      -T GenotypeGVCFs \
      -R ${ref_fasta} \
      -o ${unfiltered_vcf_filename} \
      -V ${gvcf} \
      -L ${interval_list} \
      -D ${dbsnp_vcf}
    >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.1010_with_gatk3.5"
    memory: "4 GB"
    cpu: "1"
    disks: "local-disk 300 SSD"
  }
  output {
    File unfiltered_vcf = "${unfiltered_vcf_filename}"
    File unfiltered_vcf_index = "${unfiltered_vcf_index_filename}"
  }
}

## TODO - add useInbreedingCoefficientFilter as workflow input and use true and false wdl function to input correct command

# filterExpression and filtername should only be added to the commandline if useInbreedingCoefficientFilter is true

task VariantFiltration {
  String variant_filtered_vcf_filename
  String variant_filtered_vcf_index_filename

  File unfiltered_vcf
  File unfiltered_vcf_index
  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
      -jar /usr/gitc/GenomeAnalysisTK.jar \
      -T VariantFiltration \
      --filterExpression "InbreedingCoeff < -0.3" \
      --filterName InbreedingCoeff \
      -R ${ref_fasta} \
      -o ${variant_filtered_vcf_filename} \
      -V ${unfiltered_vcf} \
      -L ${interval_list}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.1010_with_gatk3.5"
    memory: "4 GB"
    cpu: "1"
    disks: "local-disk 300 SSD"
  }
  output {
    File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "${variant_filtered_vcf_index_filename}"
  }
}

task MakeSitesOnlyVcf {
  String sites_only_vcf_filename
  String sites_only_vcf_index_filename

  File input_vcf
  File input_vcf_index

  command {
    java -Xmx4000m -jar /usr/gitc/picard.jar \
      MakeSitesOnlyVcf \
      INPUT=${input_vcf} \
      OUTPUT=${sites_only_vcf_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.1010_with_gatk3.5"
    memory: "4 GB"
    disks: "local-disk 300 SSD"
  }
  output {
    File sites_only_vcf = "${sites_only_vcf_filename}"
    File sites_only_vcf_index = "${sites_only_vcf_index_filename}"
  }
}
## TODO - add forceDeterminism as workflow output and use true and false wdl function to input correct command
# num_threads is 2 if forceDeterminism is false which is not captured here, set it to 2 for now

## TODO - Combine following two Variant Recalibrators into one once map iteration similar to array 'sep' or looping is implemented
# the -resource file inputs need to be localised and pulled into the workflow and currently this is the only viable way of doing it
# until one of the above mentioned syntax is implemented

task SNPSVariantRecalibrator {
  String recalibration_filename
  String recalibration_index_filename
  String tranches_filename
  String recalibration_plots_rscript_filename
  String mode

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  Array[File] sites_only_variant_filtered_vcf_list
  Array[File] sites_only_variant_filtered_vcf_index_list

  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index


  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16000m \
      -jar /usr/gitc/GenomeAnalysisTK.jar \
      -T VariantRecalibrator \
      -R ${ref_fasta} \
      -input ${sep=' -input 'sites_only_variant_filtered_vcf_list} \
      -L ${interval_list} \
      --num_threads 2 \
      -recalFile ${recalibration_filename} \
      -tranchesFile ${tranches_filename} \
      -allPoly \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode ${mode} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf} \
      -rscriptFile ${recalibration_plots_rscript_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.1010_with_gatk3.5"
    memory: "16 GB"
    cpu: "2"
    disks: "local-disk 300 SSD"
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_index_filename}"
    File tranches = "${tranches_filename}"
    File recalibration_plots_rscript = "${recalibration_plots_rscript_filename}"
    File recalibration_plots_pdf = "${recalibration_plots_rscript_filename}.pdf"
  }
}

task IndelsVariantRecalibrator {
  String recalibration_filename
  String recalibration_index_filename
  String tranches_filename
  String recalibration_plots_rscript_filename
  String mode

  Array[String] recalibration_tranche_values
  Array[String] recalibration_annotation_values

  Array[File] sites_only_variant_filtered_vcf_list
  Array[File] sites_only_variant_filtered_vcf_index_list

  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File dbsnp_resource_vcf
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index
  File dbsnp_resource_vcf_index

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16000m \
      -jar /usr/gitc/GenomeAnalysisTK.jar \
      -T VariantRecalibrator \
      -R ${ref_fasta} \
      -input ${sep=' -input 'sites_only_variant_filtered_vcf_list} \
      -L ${interval_list} \
      --num_threads 2 \
      -recalFile ${recalibration_filename} \
      -tranchesFile ${tranches_filename} \
      -allPoly \
      -tranche ${sep=' -tranche ' recalibration_tranche_values} \
      -an ${sep=' -an ' recalibration_annotation_values} \
      -mode ${mode} \
      -resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf} \
      -rscriptFile ${recalibration_plots_rscript_filename}
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.1010_with_gatk3.5"
    memory: "16 GB"
    cpu: "2"
    disks: "local-disk 300 SSD"
  }
  output {
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_index_filename}"
    File tranches = "${tranches_filename}"
    File recalibration_plots_rscript = "${recalibration_plots_rscript_filename}"
    File recalibration_plots_pdf = "${recalibration_plots_rscript_filename}.pdf"
  }
}

task ApplyRecalibration {
  String filtered_vcf_filename
  String filtered_vcf_index_filename
  Array[File] variant_filtered_vcf_list
  Array[File] variant_filtered_vcf_index_list
  File interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File recalibration
  File recalibration_index
  File tranches

  String mode

  Float filter_level

  command {
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m \
      -jar /usr/gitc/GenomeAnalysisTK.jar \
      -T ApplyRecalibration \
      -R ${ref_fasta} \
      -o ${filtered_vcf_filename} \
      -input ${sep=' -input 'variant_filtered_vcf_list} \
      -L ${interval_list} \
      -recalFile ${recalibration} \
      -tranchesFile ${tranches} \
      -ts_filter_level ${filter_level} \
      -mode ${mode}

  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:1.1010_with_gatk3.5"
    memory: "4 GB"
    cpu: "1"
    disks: "local-disk 300 SSD"
  }
  output {
    File filtered_vcf = "${filtered_vcf_filename}"
    File filtered_vcf_index = "${filtered_vcf_index_filename}"
  }
}

workflow JointGenotyping {
  String callset_name

  Array[String] snp_recalibration_tranche_values
  Array[String] snp_recalibration_annotation_values
  Array[String] indel_recalibration_tranche_values
  Array[String] indel_recalibration_annotation_values

  Array[File] gvcf_list
  # need inferenced index files of the input gvcfs for GenotypeGVCF and because of how input files are handled,
  # this puts the index files in the same workflow-inputs directory as the gvcf files
  Array[File] gvcf_index_list

  File call_interval_list
  File eval_interval_list
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File dbsnp_vcf
  File dbsnp_vcf_index
  File hapmap_resource_vcf
  File omni_resource_vcf
  File one_thousand_genomes_resource_vcf
  File dbsnp_resource_vcf = dbsnp_vcf
  File mills_resource_vcf
  File axiomPoly_resource_vcf
  File hapmap_resource_vcf_index
  File omni_resource_vcf_index
  File one_thousand_genomes_resource_vcf_index
  File dbsnp_resource_vcf_index = dbsnp_vcf_index
  File mills_resource_vcf_index
  File axiomPoly_resource_vcf_index

  Float snp_filter_level
  Float indel_filter_level

  # this maps the genomic position based on contigs to the tiledb interpretation of the same position
  # chr 3 pos 300 for tiledb would be ((size of chr1 + size of chr2 + 300) - 1)  (0-based)
  File vid_mapping_file
  File vcf_header
  Int num_of_tiledb_partitions
  Int size_per_column_partition

  # this only has to be called once per callset to get appropriately sized partition intervals.  appropriate here
  # means that each interval will have equivalent numbers of SNPS
#  call TileDBHistogram {
#    input:
#      size_per_column_partition = size_per_column_partition,
#      vid_mapping_file = vid_mapping_file,
#      gvcf_list = gvcf_list,
#      gvcf_index_list = gvcf_index_list,
#      ref_fasta = ref_fasta,
#      ref_fasta_index = ref_fasta_index
#  }

  call CreateTileDBPartitions {
    input:
      num_of_tiledb_partitions = num_of_tiledb_partitions
  }

  scatter (partition in CreateTileDBPartitions.partition_json) {
    call TileDBCombineGVCF {
      input:
        output_gvcf_filename = "output.vcf.gz",
        gvcf_list = gvcf_list,
        gvcf_index_list = gvcf_index_list,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        partition = partition,
        size_per_column_partition = size_per_column_partition,
        vid_mapping_file = vid_mapping_file,
        vcf_header = vcf_header
    }

    call TabixAndGenotypeGVCF {
      input:
        unfiltered_vcf_filename = callset_name + ".unfiltered.vcf.gz",
        unfiltered_vcf_index_filename = callset_name + ".unfiltered.vcf.gz.tbi",
        gvcf = TileDBCombineGVCF.output_gvcf,
        interval_list = call_interval_list,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index
    }

    call VariantFiltration {
      input:
        variant_filtered_vcf_filename = callset_name + ".variant_filtered.vcf.gz",
        variant_filtered_vcf_index_filename = callset_name + ".variant_filtered.vcf.gz.tbi",
        unfiltered_vcf = TabixAndGenotypeGVCF.unfiltered_vcf,
        unfiltered_vcf_index = TabixAndGenotypeGVCF.unfiltered_vcf_index,
        interval_list = call_interval_list,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict
    }

    call MakeSitesOnlyVcf {
      input:
        input_vcf = VariantFiltration.variant_filtered_vcf,
        input_vcf_index = VariantFiltration.variant_filtered_vcf_index,
        sites_only_vcf_filename = callset_name + ".sites_only.variant_filtered.vcf.gz",
        sites_only_vcf_index_filename = callset_name + ".sites_only.variant_filtered.vcf.gz.tbi"
    }
  }

  call SNPSVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf_list = MakeSitesOnlyVcf.sites_only_vcf,
      sites_only_variant_filtered_vcf_index_list = MakeSitesOnlyVcf.sites_only_vcf_index,
      interval_list = call_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      recalibration_filename = callset_name + ".snps.recal",
      recalibration_index_filename = callset_name + ".snps.recal.idx",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_plots_rscript_filename = callset_name + ".snps.recalibration_plots.rscript",
      recalibration_tranche_values = snp_recalibration_tranche_values,
      recalibration_annotation_values = snp_recalibration_annotation_values,
      mode = "SNP",
      hapmap_resource_vcf = hapmap_resource_vcf,
      hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      omni_resource_vcf = omni_resource_vcf,
      omni_resource_vcf_index = omni_resource_vcf_index,
      one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
  }

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf_list = MakeSitesOnlyVcf.sites_only_vcf,
      sites_only_variant_filtered_vcf_index_list = MakeSitesOnlyVcf.sites_only_vcf_index,
      interval_list = call_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      recalibration_filename = callset_name + ".indels.recal",
      recalibration_index_filename = callset_name + ".indels.recal.idx",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_plots_rscript_filename = callset_name + ".indels.recalibration_plots.rscript",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mode = "INDEL",
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
  }

  call ApplyRecalibration as IndelsApplyRecalibration {
    input:
      filtered_vcf_filename = callset_name + ".tmp_indel_filtered.vcf.gz",
      filtered_vcf_index_filename = callset_name +".tmp_indel_filtered.vcf.gz.tbi",
      variant_filtered_vcf_list = VariantFiltration.variant_filtered_vcf,
      variant_filtered_vcf_index_list = VariantFiltration.variant_filtered_vcf_index,
      interval_list = call_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      recalibration = IndelsVariantRecalibrator.recalibration,
      recalibration_index = IndelsVariantRecalibrator.recalibration_index,
      tranches = IndelsVariantRecalibrator.tranches,
      mode = "INDEL",
      filter_level = indel_filter_level
  }

  call ApplyRecalibration as SNPSApplyRecalibration {
    input:
      filtered_vcf_filename = callset_name + ".filtered.vcf.gz",
      filtered_vcf_index_filename = callset_name + ".filtered.vcf.gz.tbi",
      variant_filtered_vcf_list = IndelsApplyRecalibration.filtered_vcf,
      variant_filtered_vcf_index_list = IndelsApplyRecalibration.filtered_vcf_index,
      interval_list = call_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      recalibration = SNPSVariantRecalibrator.recalibration,
      recalibration_index = IndelsVariantRecalibrator.recalibration_index,
      tranches = SNPSVariantRecalibrator.tranches,
      mode = "SNP",
      filter_level = snp_filter_level
  }

  output {
    SNPSApplyRecalibration.*
  }
}