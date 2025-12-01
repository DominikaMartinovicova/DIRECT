
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = 'data/'
log_dir = 'logs/'
TMA = ['slide_3', 'slide_4', 'slide_5', 'slide_6']
Njob = 4
phenotyping_level = ['Salcher_celltypes', 'Neutro_Epi_extImm']
#'extNeutro_extEpi_extImm',
#'extNeutro_Epi_extImm',
#'Neutro_Epi_extImm',
#'extEpi_extImm',
#'Epi_extImm']

#-------------------------------------------------------------------------------------------------------------------
# 0.1 specify target rules
rule all:
    input:
        #checked_adatas = data_dir + "combined/{phenotyping_level}_checked_adatas.h5ad"
        expand(data_dir + "combined/{phenotyping_level}_checked_adatas.h5ad", phenotyping_level=phenotyping_level)
        #expand(data_dir + "phenotyped/{TMA}.zarr/zmetadata", TMA=TMA)
        #data_dir + "combined/checked_adatas.h5ad"
        #expand(data_dir + "preprocessed/{TMA}.zarr/zmetadata", TMA=TMA),
        # expand(data_dir + "phenotyped/{TMA}.zarr/zmetadata", TMA=TMA),
        # data_dir + "combined/combined_adatas.h5ad",
        # data_dir + "combined/checked_adatas.h5ad"

#++++++++++++++++++++++++++++++++++++++++++ 1 Preprocess Xenium data +++++++++++++++++++++++++++++++++++++++++	
# 1.1 Compress Xenium data to zarr
rule compress_Xenium:
    input:
        raw_Xenium = data_dir + "raw/{TMA}/cell_feature_matrix.h5"
    output:
        compressed_Xenium = data_dir + "compressed/{TMA}.zarr/zmetadata"
    params:
        in_dir = data_dir + "raw/{TMA}/",
        out_dir = data_dir + "compressed/{TMA}.zarr/"
    threads: Njob
    #conda:
    #    "envs/squidpy.yaml"
    shell:
        """
        python3 scripts/python/compress_Xenium.py \
        -i {input.raw_Xenium} \
        --input_dir {params.in_dir} \
        --output_dir {params.out_dir} \
        --threads {threads} \
        -o {output.compressed_Xenium} 
        """

# 1.2 Do QA and add experimental info (i.e. these cells belong to core1=T.....=biopsy)
rule QC_Xenium:
    input:
        compressed_Xenium = data_dir + "compressed/{TMA}.zarr/zmetadata",   # input compressed Xenium data
        coordinates = data_dir + "coordinates/coordinates_{TMA}.csv",       # coordinates to assign cells to their coresponding samples
        metadata = data_dir + "raw/{TMA}/metadata_{TMA}.csv"                # metadata regarding clinical/experimental info about samples
    output:
        preprocessed_Xenium = data_dir + "preprocessed/{TMA}.zarr/zmetadata",   # output preprocessed Xenium data
        quality_plots = directory('plots/preprocessing/{TMA}/QA/'),             # output plots for QA
        spatial_plots = directory('plots/preprocessing/{TMA}/spatial/')         # output spatial plots to check metadata assignment
    params:
        in_dir = data_dir + "compressed/{TMA}.zarr/",
        out_dir = data_dir + "preprocessed/{TMA}.zarr/"
    #conda:
    #    "envs/squidpy.yaml"
    shell:
        """
        python3 scripts/python/QC_Xenium.py \
        -i {input.compressed_Xenium} \
        --input_coords {input.coordinates} \
        --input_meta {input.metadata} \
        --input_dir {params.in_dir} \
        --output_dir {params.out_dir} \
        -o {output.preprocessed_Xenium} \
        --output_plot_QA {output.quality_plots} \
        --output_plot_spatial {output.spatial_plots}
        """


# 1.3 Run TACCO for cell types
rule Tacco:
    input:
        preprocessed_Xenium = data_dir + "preprocessed/{TMA}.zarr/zmetadata",
        scRNAseq_atlas = '/net/beegfs/groups/tgac/dmartinovicova_new/NSCLC/scRNAseq/data/final_scRNAseq_atlas_Salcher.h5ad'
    output:
        phenotyped_Xenium = data_dir + "phenotyped/{phenotyping_level}/{TMA}.zarr/zmetadata",
        plots_dir = directory('plots/tacco/{phenotyping_level}/{TMA}/')
    params:
        in_dir = data_dir + "preprocessed/{TMA}.zarr/",
        out_dir = data_dir + "phenotyped/{phenotyping_level}/{TMA}.zarr/",
        phen_level = "{phenotyping_level}"
    #conda:
    #    "envs/tacco.yaml"
    shell:
        """
        python3 scripts/python/run_tacco.py \
        -i {input.preprocessed_Xenium} \
        --input_atlas {input.scRNAseq_atlas} \
        --input_dir {params.in_dir} \
        --output_dir {params.out_dir} \
        --phen_level {params.phen_level} \
        -o {output.phenotyped_Xenium} \
        --output_dir_plot {output.plots_dir}
        """

# 1.4 Combine adatas from all TMAs into one adata
rule combine_adatas:
    input:
        phenotyped_Xenium = expand(data_dir + "phenotyped/{phenotyping_level}/{TMA}.zarr/zmetadata", TMA = TMA, phenotyping_level = phenotyping_level)
    output:
        combined_adatas = data_dir + "combined/{phenotyping_level}_combined_adatas.h5ad",
        output_plots = directory('plots/combined/{phenotyping_level}/')
    #conda:
    #    "envs=""
    params:
        in_dir = data_dir + "phenotyped/{phenotyping_level}/",
    shell:
        """
        python3 scripts/python/combine_adatas.py \
	    --input_dir {params.in_dir} \
        -o {output.combined_adatas} \
        --output_plot {output.output_plots}
        """

# 1.5 Check phenotyping based on marker genes
rule check_phenotyping:
    input:
        combined_adatas = data_dir + "combined/{phenotyping_level}_combined_adatas.h5ad"
    output:
        checked_adatas = data_dir + "combined/{phenotyping_level}_checked_adatas.h5ad"
    #conda:
    #    "envs=""
    params:
        out_plot_dir = 'plots/combined/{phenotyping_level}/',
        phen_level = "{phenotyping_level}"
    shell:
        """
        python3 scripts/python/check_phenotyping.py \
	    -i {input.combined_adatas} \
        --phen_level {params.phen_level} \
        -o {output.checked_adatas} \
        --output_plot {params.out_plot_dir}
        """
#++++++++++++++++++++++++++++++++++++++++++ 2 Downstream analysis +++++++++++++++++++++++++++++++++++++++++	



