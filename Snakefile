configfile: "config.yaml"

import os
from os.path import basename, join

WD = config['wd']

FASTQ_DIR = join(WD, config['fastq_dir'])
REFGEN_DIR = join(WD, config['refgen_dir'])
GENCODE_DIR = join(WD, config['gencode_dir'])
VCF_DIR = join(WD, config['vcf_dir'])

SNP_DIR = join(WD, config['snp_dir'])
SAM_DIR = join(WD, config['sam_dir'])
READ_GROUPS_DIR = join(WD, config['add_rg_dir'])
MARK_DUPS_DIR = join(WD, config['mark_dups_dir'])
FIND_SNPS_DIR = join(WD, config['find_snps_dir'])
REMAP_DIR = join(WD, config['remap_dir'])
REMAP_RG_DIR = join(WD, config['remap_add_rg_dir'])
REMAP_MD_DIR = join(WD, config['remap_mark_dups_dir'])
FILT_REMAPPED_READS_DIR = join(WD, config['filter_remapped_reads'])
MERGED_DIR = join(WD, config['merged_dir'])
SORTED_DIR = join(WD, config['sorted_dir'])
READ_COUNTS_DIR = join(WD, config['read_counts_dir'])


WC_fastqs = glob_wildcards(join(FASTQ_DIR, '{sample}.{pair}.fq.gz'))
WC_refgens = glob_wildcards(join(REFGEN_DIR, '{genome}.fa'))
WC_gencodes = glob_wildcards(join(GENCODE_DIR, '{gencode}.gtf'))

SAMPLES = set(WC_fastqs.sample)
PAIRS = ['R1', 'R2']

GENOMES = set(WC_refgens.genome)
GENCODES = set(WC_gencodes.gencode)


rule all:
    input:
        expand('{read_counts_dir}/{{sample}}.{{genome}}.tsv'.format(read_counts_dir=READ_COUNTS_DIR),
        sample=SAMPLES, genome=GENOMES)
    run:
        print("GENOMIC_ASECOUNT FINISHED WITH NO EXCEPTIONS!")


rule bowtie2_index_genome:
    input:
        "{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
    output:
        "{refgen_dir}/{{genome}}.fa.1.bt2".format(refgen_dir=REFGEN_DIR)
    threads: config['bowtie2_build_threads']
    shell:
        "bowtie2-build --threads {threads} -f {input} {input}"


rule bowtie2_map:
    input:
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
        idx="{refgen_dir}/{{genome}}.fa.1.bt2".format(refgen_dir=REFGEN_DIR),
        f1='{fastq_dir}/{{sample}}.R1.fq.gz'.format(fastq_dir=FASTQ_DIR),
        f2='{fastq_dir}/{{sample}}.R2.fq.gz'.format(fastq_dir=FASTQ_DIR)
    output:
        '{sam_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.sam'.format(sam_dir=SAM_DIR)
    threads: config['bowtie2_map_threads']
    shell:
        "bowtie2 -x {input.refgen} -1 {input.f1} -2 {input.f2} -S {output} -p {threads}"


rule add_read_groups:
    input:
        '{sam_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.sam'.format(sam_dir=SAM_DIR)
    output:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=READ_GROUPS_DIR)
    run:
        command = "{picard} AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=id " \
                  "RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample; ".format(picard=config['picard_path'],
                  input=input, output=output)

        print(command)
        shell(command)


rule mark_duplicates:
    input:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=READ_GROUPS_DIR)
    output:
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=MARK_DUPS_DIR),
        metrics='{mark_dups_dir}/{{sample}}.{{genome}}.metrics'.format(mark_dups_dir=MARK_DUPS_DIR)
    run:
        command = "{picard} MarkDuplicates I={input} O={output_bam} CREATE_INDEX=true " \
                  "VALIDATION_STRINGENCY=SILENT M={output_metrics}".format(picard=config['picard_path'],
                  input=input, output_bam=output.bam, output_metrics=output.metrics)

        print(command)
        shell(command)


rule make_snp_dir:
    input:
        '{vcf_dir}/{{sample}}.{{genome}}.vcf'.format(vcf_dir=VCF_DIR)
    output:
        '{snp_dir}/{{sample}}.{{genome}}'.format(snp_dir=SNP_DIR)
    run:
        import vcf, gzip
        os.makedirs(str(output))
        reader = vcf.Reader(filename=str(input))


        snp_files = {}
        for rec in reader:
            if rec.CHROM not in snp_files.keys():
                snp_files[rec.CHROM] = gzip.open(join(str(output), rec.CHROM+".snps.txt.gz"), 'wt')

            snp_files[rec.CHROM].write(' '.join(map(str, [rec.POS, rec.REF, rec.ALT[0]]))+'\n')

        for chrom, out in snp_files.items():
            out.close()


rule find_intersecting_snps:
    input:
        snps='{snp_dir}/{{sample}}.{{genome}}'.format(snp_dir=SNP_DIR),
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=MARK_DUPS_DIR)
    output:
        dir='{find_snps_dir}/{{sample}}.{{genome}}'.format(find_snps_dir=FIND_SNPS_DIR),
        keep_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.keep.bam'.format(find_snps_dir=FIND_SNPS_DIR),
        remap_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.to.remap.bam'.format(find_snps_dir=FIND_SNPS_DIR),
        fq1='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq1.gz'.format(find_snps_dir=FIND_SNPS_DIR),
        fq2='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq2.gz'.format(find_snps_dir=FIND_SNPS_DIR),
        single='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.single.fq.gz'.format(find_snps_dir=FIND_SNPS_DIR)
    conda:
        "envs/wasp.yaml"
    shell:
        "python scripts/find_intersecting_snps.py --is_paired_end --is_sorted "
        "--output_dir {output.dir} --snp_dir {input.snps} {input.bam}"


rule bowtie2_remap:
    input:
        refgen="{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR),
        idx="{refgen_dir}/{{genome}}.fa.1.bt2".format(refgen_dir=REFGEN_DIR),
        f1='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq1.gz'.format(find_snps_dir=FIND_SNPS_DIR),
        f2='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.remap.fq2.gz'.format(find_snps_dir=FIND_SNPS_DIR)
    output:
        '{remap_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.sam'.format(remap_dir=REMAP_DIR)
    threads: config['bowtie2_map_threads']
    shell:
        "bowtie2 -x {input.refgen} -1 {input.f1} -2 {input.f2} -S {output} -p {threads}"


rule remap_add_read_groups:
    input:
        '{remap_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.sam'.format(remap_dir=REMAP_DIR)
    output:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=REMAP_RG_DIR)
    run:
        command = "{picard} AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=id " \
                  "RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample; ".format(picard=config['picard_path'],
                  input=input, output=output)

        print(command)
        shell(command)


rule remap_mark_duplicates:
    input:
        '{add_rg_dir}/{{sample}}.{{genome}}.bam'.format(add_rg_dir=REMAP_RG_DIR)
    output:
        bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=REMAP_MD_DIR),
        metrics='{mark_dups_dir}/{{sample}}.{{genome}}.metrics'.format(mark_dups_dir=REMAP_MD_DIR)
    run:
        command = "{picard} MarkDuplicates I={input} O={output_bam} CREATE_INDEX=true " \
                  "VALIDATION_STRINGENCY=SILENT M={output_metrics}".format(picard=config['picard_path'],
                  input=input, output_bam=output.bam, output_metrics=output.metrics)

        print(command)
        shell(command)


rule filter_remapped_reads:
    input:
        to_remap_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.to.remap.bam'.format(find_snps_dir=FIND_SNPS_DIR),
        remap_bam='{mark_dups_dir}/{{sample}}.{{genome}}.bam'.format(mark_dups_dir=REMAP_MD_DIR),
    output:
        '{filt_remapped_reads_dir}/{{sample}}.{{genome}}.bam'.format(filt_remapped_reads_dir=FILT_REMAPPED_READS_DIR)
    conda:
        "envs/wasp.yaml"
    shell:
        'python scripts/filter_remapped_reads.py {input.to_remap_bam} {input.remap_bam} {output}'


rule samtools_merge:
    input:
        filtered_bam='{filt_remapped_reads_dir}/{{sample}}.{{genome}}.bam'.format(filt_remapped_reads_dir=FILT_REMAPPED_READS_DIR),
        keep_bam='{find_snps_dir}/{{sample}}.{{genome}}/{{sample}}.{{genome}}.keep.bam'.format(find_snps_dir=FIND_SNPS_DIR)
    output:
        '{merged_dir}/{{sample}}.{{genome}}.bam'.format(merged_dir=MERGED_DIR)
    shell:
        'samtools merge {output} {input.filtered_bam} {input.keep_bam}'


rule samtools_sort:
    input:
        '{merged_dir}/{{sample}}.{{genome}}.bam'.format(merged_dir=MERGED_DIR)
    output:
        '{sorted_dir}/{{sample}}.{{genome}}.bam'.format(sorted_dir=SORTED_DIR)
    shell:
        'samtools sort {input} > {output}'


rule samtools_index:
    input:
        '{sorted_dir}/{{sample}}.{{genome}}.bam'.format(sorted_dir=SORTED_DIR)
    output:
        '{sorted_dir}/{{sample}}.{{genome}}.bam.bai'.format(sorted_dir=SORTED_DIR)
    shell:
        'samtools index {input}'


rule index_genome:
    input:
        "{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        "{refgen_dir}/{{genome}}.fa.fai".format(refgen_dir=REFGEN_DIR)
    shell:
        "samtools faidx {input}"


rule create_genome_dictionary:
    input:
        "{refgen_dir}/{{genome}}.fa".format(refgen_dir=REFGEN_DIR)
    output:
        "{refgen_dir}/{{genome}}.dict".format(refgen_dir=REFGEN_DIR)
    run:
        command = "{picard} CreateSequenceDictionary R={input} O={output}".format(picard=config['picard_path'],
                  input=input, output=output)

        print(command)
        shell(command)


rule ase_read_counter:
    input:
        bam='{sorted_dir}/{{sample}}.{{genome}}.bam'.format(sorted_dir=SORTED_DIR),
        idx='{sorted_dir}/{{sample}}.{{genome}}.bam.bai'.format(sorted_dir=SORTED_DIR),
        ref='{refgen_dir}/{{genome}}.fa'.format(refgen_dir=REFGEN_DIR),
        refidx="{refgen_dir}/{{genome}}.fa.fai".format(refgen_dir=REFGEN_DIR),
        refdict="{refgen_dir}/{{genome}}.dict".format(refgen_dir=REFGEN_DIR),
        vcf='{vcf_dir}/{{sample}}.{{genome}}.vcf'.format(vcf_dir=VCF_DIR)
    output:
        '{read_counts_dir}/{{sample}}.{{genome}}.tsv'.format(read_counts_dir=READ_COUNTS_DIR)
    run:
        command = "{gatk} -T ASEReadCounter -I {input} -o {output} -sites {vcf} -R {refgen} -U ALLOW_N_CIGAR_READS -minDepth -1 " \
                  "--minBaseQuality -1 --minMappingQuality -1 -drf DuplicateRead".format(
                  gatk=config['gatk_path'], input=str(input.bam), vcf=input.vcf, output=str(output),
                  refgen=str(input.ref))

        print(command)
        shell(command)
