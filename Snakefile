# Configuration of the pipeline
# configfile: "config.yaml"
data_dir = config["data_dir"]
out_dir = config["out_dir"]
distances_file = config["distances_file"]
# Use glob_wildcards to find all VCF files in the directory
# This will match files like '/path/to/data/file1.vcf.gz', '/path/to/data/file2.vcf.gz', etc.
(names,) = glob_wildcards(f"{data_dir}{{name}}.vcf.gz")
# We require one map file in plink format for each VCF file
# We expect map files being in the same directory as the VCF files
suffix_map = config.get("map_file_suffix", ".map")
(map_files,) = glob_wildcards(f"{data_dir}{{name}}{suffix_map}")

contig_lengths = config.get("contig_lengths")
if contig_lengths is not None:
    assert all(
        isinstance(x, float) for x in contig_lengths
    ), "Contig lengths must be a list of floats"
    outfiles = [
        f"{out_dir}/ibd_dispersal_data.csv",
        f"{out_dir}/constant_density_mle.csv",
    ]
else:
    outfiles = [f"{out_dir}/ibd_dispersal_data.csv"]


rule all:
    input:
        outfiles,


# IBD detection
rule download_hapibd:
    output:
        "resources/hap-ibd.jar",
    shell:
        """
        curl -L -o {output} https://faculty.washington.edu/browning/hap-ibd.jar
        """


rule download_merge_ibd:
    output:
        "resources/merge-ibd.jar",
    shell:
        """
        curl -L -o {output} https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar
        """


# Run HapIBD detection software
rule run_hapibd:
    input:
        vcf=f"{data_dir}{{name}}.vcf.gz",
        map=f"{data_dir}{{name}}{suffix_map}",
        jar="resources/hap-ibd.jar",
    output:
        "steps/{name}.ibd.gz",
        "steps/{name}.hbd.gz",
    log:
        "logs/hapibd/{name}.log",
    shadow:
        "minimal"
    threads: 2
    conda:
        "conda_environment.yaml"
    params:
        extra=config.get("params", ""),
    shell:
        """
        # java -jar hap-ibd.jar gt=$file map=plink.chr$CHR.GRCh38.map  out=IBD/$PREFIX
        java -jar {input.jar} \
            gt={input.vcf} map={input.map} out=steps/{wildcards.name} \
            nthreads={threads} {params.extra}
        mv steps/{wildcards.name}.log {log}
        """


# Run post-processing of IBD
rule post_processing_ibd:
    input:
        vcf=f"{data_dir}{{name}}.vcf.gz",
        map=f"{data_dir}{{name}}{suffix_map}",
        ibd="steps/{name}.ibd.gz",
        merge_jar="resources/merge-ibd.jar",
    output:
        "steps/{name}.postprocessed.ibd.gz",
    log:
        "logs/postprocessing_ibd/{name}.log",
    shadow:
        "minimal"
    params:
        gap=config.get("gap", 0.6),  # in cM
        discord=config.get("discord", 1),  # at most one discordant homozygote
    conda:
        "conda_environment.yaml"
    shell:
        """
        echo "Processing {wildcards.name}" > {log}
        echo "Processing IBD file" >> {log}
        gunzip -c {input.ibd} | tee >(wc -l > {log}.ibd_lines) | \
            java -jar {input.merge_jar} {input.vcf} {input.map} \
            {params.gap} {params.discord} > steps/{wildcards.name}.postprocessed.ibd 2>> {log}
        # Count final output lines
        wc -l steps/{wildcards.name}.postprocessed.ibd > {log}.final_lines
        # Log everything
        echo "IBD lines: $(cat {log}.ibd_lines)" >> {log}
        echo "Final lines: $(cat {log}.final_lines)" >> {log}
        # Compression
        gzip steps/{wildcards.name}.postprocessed.ibd
        echo "Done" >> {log}
        # Clean up temp files
        rm {log}.ibd_lines {log}.final_lines
        """


# Postprocessing of the IBD data and the individual distances
rule postprocess:
    input:
        ibds=expand("steps/{name}.postprocessed.ibd.gz", name=names),
        distances=distances_file,
    output:
        f"{out_dir}/ibd_dispersal_data.csv",
    log:
        "logs/ibd_dispersal_data.log",
    shell:
        """
        julia -e '
        import Pkg
        Pkg.activate(; temp = true)
        Pkg.add(Pkg.PackageSpec(name="CSV", version="0.10.15"))
        Pkg.add(Pkg.PackageSpec(name="DataFrames", version="1.7.0"))
        Pkg.add(Pkg.PackageSpec(name="IdentityByDescentDispersal", version="0.1.2"))
        using CSV, DataFrames, IdentityByDescentDispersal
        colnames = ["ID1", "HAP1", "ID2", "HAP2", "CHR", "START", "END", "SCORE", "LENGTH"]
        # Read IBD files and concatenate them:
        files = split("{input.ibds}")
        println("Number of IBD files: ", length(files))
        dfs = [CSV.read(f, DataFrame; header=colnames) for f in files]
        df_ibds = vcat(dfs...)
        df_ibds.span = df_ibds.LENGTH ./ 100
        # Read distance files
        df_distances = CSV.read("{input.distances}", DataFrame)
        # Merge both dataframes with default bins
        bins, min_threshold = default_ibd_bins()
        df_preprocessed = preprocess_dataset(df_ibds, df_distances, bins, min_threshold);
        CSV.write("{output}", df_preprocessed)
        ' &> {log}
        """


# MLE estimate example assuming a constant density


rule mle_constant_density:
    input:
        f"{out_dir}/ibd_dispersal_data.csv",
    output:
        f"{out_dir}/constant_density_mle.csv",
    log:
        "logs/constant_density_mle",
    params:
        contig_lengths=contig_lengths,
    shell:
        """
        julia -e '
        # Parse contig lengths
        contig_lengths = "{params.contig_lengths}"
        contig_lengths = parse.(Float64, split(contig_lengths, " "))
        println("Parsed contig lengths:")
        println(contig_lengths)
        import Pkg
        Pkg.activate(; temp = true)
        Pkg.add(Pkg.PackageSpec(name="CSV", version="0.10.15"))
        Pkg.add(Pkg.PackageSpec(name="DataFrames", version="1.7.0"))
        Pkg.add(Pkg.PackageSpec(name="Turing", version="0.39.1"))
        Pkg.add(Pkg.PackageSpec(name="StatsBase", version="0.34.5"))
        Pkg.add(Pkg.PackageSpec(name="IdentityByDescentDispersal", version="0.1.2"))
        using CSV, DataFrames, IdentityByDescentDispersal, Turing, StatsBase
        df = CSV.read("{input}", DataFrame)
        @model function constant_density(df, contig_lengths)
            D ~ Uniform(0, 1e8)
            σ ~ Uniform(0, 1e8)
            Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
        end
        m = constant_density(df, contig_lengths);
        mle_estimate = maximum_likelihood(m)
        println("MLE Estimate:")
        println(mle_estimate)
        coef_table = mle_estimate |> coeftable |> DataFrame
        CSV.write("{output}", coef_table)
        ' &> {log}
        """
