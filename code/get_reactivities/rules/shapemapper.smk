"""
Description:
Specific arrangement of rules to create the workflow

Usage:
"""

import pandas as pd
import glob

rule calculate_length_fasta:
    input:
        handle_file = os.path.join(config['project_path'], config['data_folder'], config['trna_reference_original']),
    output:
        handle_file_len = os.path.join(config['project_path'], config['data_folder'], config['trna_reference_original_len']),
    run:
        header = None
        length = 0
        with open(output.handle_file_len, 'w') as outfile:
            with open(input.handle_file, 'r') as fasta:
                for line in fasta:
                    line = line.rstrip()
                    if line.startswith('>'):
                        # If we captured one before, print it now
                        if header is not None:
                            outfile.write(str(header)+"\t"+str(length)+"\n")
                            length = 0
                        header = line[1:]
                    else:
                        length += len(line)
            # Don't forget the last one
            if length:
                outfile.write(str(header)+"\t"+str(length)+"\n")

rule convert_sam_old:
    input:
        sam_file = os.path.join(config['project_path'], config['data_folder'], "{sample}.sam"),
    output:
        sam_file = os.path.join(config['project_path'], config['data_folder'], config['old_sam'], "{sample}.sam"),
    shell:
        """
        reformat.sh in={input.sam_file} out={output.sam_file} sam=1.3
        """

rule snakemake_create_list_references:
    """
    Take reference and get header names as a list
    """
    input:
        target_fasta = os.path.join(config['project_path'], config['data_folder'], config['trna_reference']),
    output:
        list_references = os.path.join(config['project_path'], config['data_folder'], "list_references.txt"),
    shell:
        """
        grep -o '^>.*' {input.target_fasta} | cut -c 2- | sort -u > {output.list_references}
        """

rule snakemake_split_reference_sequence:
    """
    Split tRNA reference by individual fasta files
    """
    input:
        target_fasta = os.path.join(config['project_path'], config['data_folder'], config['trna_reference']),
    output:
        success_file = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "success_split_files.txt"),
    params:
        reference_folder = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder']),
    shell:
        """
        faSplit byname {input.target_fasta} {params.reference_folder}/ && echo "0" > {output.success_file} || echo "Command failed"
        """

rule split_sam_by_reference:
    """
    Split SAM files by reference (column 3) using samtools + awk.
    Each reference will produce a separate SAM file named <reference>.sam
    inside the specified split folder.
    """
    input:
        sam_file = os.path.join(config['project_path'], config['data_folder'], config['old_sam'], "{sample}.sam"),
    output:
        success_file = os.path.join(config['project_path'], config['data_folder'],
                                    config['split_alignment_file'], "success_split_sam_{sample}.txt"),
    params:
        split_folder = os.path.join(config['project_path'], config['data_folder'],
                                    config['split_alignment_file']),
    log:
        os.path.join(config['project_path'], config['data_folder'], "split_sam_{sample}.log"),
    threads: 2
    shell:
        r"""
        mkdir -p {params.split_folder}
        cd {params.split_folder}

        # Split SAM records by reference name (field 3)
        samtools view {input.sam_file} | \
            awk -F'\t' '{{print $0 > $3"_{wildcards.sample}.sam"}}' &> {log}

        # Check for success
        if [ $? -eq 0 ]; then
            echo "0" > {output.success_file}
        else
            echo "1" > {output.success_file}
            exit 1
        fi
        """

#rule snakemake_split_alignment_files:
    #"""
    #Split SAM files (non-extended format) files based on query calculated references
    #"""
    #input:
        #list_references = os.path.join(config['project_path'], config['data_folder'], "list_references.txt"),
        #sam_file = os.path.join(config['project_path'], config['data_folder'], "{sample}.sam"),
    #output:
        #success_file = os.path.join(config['project_path'], config['data_folder'], config['split_alignment_file'], "success_split_sam_{sample}.txt"),
    #params:
        #reference_folder = os.path.join(config['project_path'], config['data_folder'], config['split_alignment_file']),
    #run:
        #import os
        #import subprocess
        #import shlex

        ## Create reference folder if it does not exist
        #os.makedirs(params.reference_folder, exist_ok=True)

        ## Read the list of references
        #with open(input.list_references) as ref_file:
            #references = [line.strip() for line in ref_file]

        #success = True
        #for ref in references:
            #temp_output = os.path.join(params.reference_folder, f"{ref}_{wildcards.sample}.sam")
            #awk_script = "$3 == ref {print $0}"
            #command = (
                #f"awk -v ref={ref_q} "
                #f"'{awk_script}' "
                #f"{infile_q} > {out_q}"
            #)
            ##command = f"awk -v ref={shlex.quote(ref)} '$3 == ref {{print $0}}' {shlex.quote(input.sam_file)} > {shlex.quote(temp_output)}"
            ##command = f"awk -v ref={ref} '$3 == ref {{print $0}}' {input.sam_file} > {temp_output}"
            #result = subprocess.run(command, shell=True)
            #if result.returncode != 0:
                #success = False
                #break

        ## Create the success file if all references were processed successfully
        #if success:
            #with open(output.success_file, 'w') as f:
                #f.write("0")

# rule snakemake_split_alignment_files:
#     """
#     Split SAM files (non-extended format) files based on query calculated references
#     """
#     input:
#         list_references = os.path.join(config['project_path'], config['data_folder'], "list_references.txt"),
#         sam_file = os.path.join(config['project_path'], config['output_folder'], config['mapping_folder_5'], "{sample}.sam"),
#     output:
#         success_file = os.path.join(config['project_path'], config['data_folder'], config['split_alignment_file_v2'], "success_split_sam_{sample}.txt"),
#     params:
#         reference_folder = os.path.join(config['project_path'], config['data_folder'], config['split_alignment_file_v2']),
#     shell:
#         """
#         mkdir -p {params.reference_folder}
#         while read -r line; do
#            awk -v ref=$line '$3 == ref {{print $0}}' {input.sam_file} > {params.reference_folder}/${{line}}_{wildcards.sample}.sam && echo "0" > {output.success_file} || echo "Failed"
#         done < {input.list_references}
#         """

rule reduce_reference_list:
    """
    Get complete list of references and print those with SAM file from DMSO >0
    """
    input:
        list_references = os.path.join(config['project_path'], config['data_folder'], "list_references.txt"),
        list_references_short = os.path.join(config['project_path'], config['data_folder'], "list_references_short.txt"),
    output:
        list_references_short_temp = os.path.join(config['project_path'], config['data_folder'], "list_references_short_temp.txt"),
    params:
        reference_folder = os.path.join(config['project_path'], config['data_folder'], config['split_alignment_file']),
    run:
        import os 

        seen_references = set()
        with open(input.list_references, 'r') as file, open(output.list_references_short_temp, 'w') as out_file:
            for ref in file:
                ref = ref.strip()
                if ref in seen_references:
                    continue

                for tto in DMSO:
                    sam_file = f"{ref}_{tto}.sam"
                    complete_sam_file = os.path.join(params.reference_folder, sam_file)

                    if os.path.exists(complete_sam_file) and os.path.getsize(complete_sam_file) > 0:
                        out_file.write(ref + "\n")
                        seen_references.add(ref)
                        break

        with open(output.list_references_short_temp, 'r') as temp_file, open(input.list_references_short, 'a') as out_file_concat:
            for ref in temp_file:
                ref = ref.strip()
                out_file_concat.write(ref+"\n")
        

rule shapemapper_parse:
    """
    NOTE: Some references didn't get a mapping reads. Such cases must be check and create a
    empty file to not broke snakemake.
        $program_parser --input_is_unpaired --in $align_file1 -o $outfolder/$out_namem.mut -m 0 -d $outfolder/$out_namem.log
          --input_is_unpaired                   specify that reads are unpaired (as 
                                                opposed to paired and/or unmerged 
                                                paired reads)
          -m [ --min_mapq ] arg (=30)           minimum reported mapping quality to 
                                                allow
          Here the sam quality should not be a discriminant becasuse segemenhl didn't report that as bowtie2-> -m = 0
          NOTE: important to look only references with mapped reads and add --latency-wait 20 when run snakemake
    """
    input:
        alignment_file = os.path.join(config['project_path'], config['data_folder'], config['split_alignment_file'],"{reference_seq}_{sample}.sam"),
    output:
        mutations_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}","{reference_seq}_{sample}.mut"),
        log_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}","{reference_seq}_{sample}.log"),
    run:
        import subprocess
        import os
        
        if os.path.exists(input.alignment_file) and os.path.getsize(input.alignment_file) > 0:
            call_args = ['shapemapper_mutation_parser',
                         '--input_is_unpaired',
                         '--in', input.alignment_file,
                         '-o', output.mutations_file,
                         '-m', '0',
                         '-w',
                         '-d', output.log_file]
            subprocess.run(call_args)
        else:
            m = open(output.mutations_file, "x")
            l = open(output.log_file, "x")
            m.close()
            l.close()

# Index Sequence length file
def build_len_dict(files):
    seq_len = {}
    with open(files, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                key, value = line.split("\t")
                seq_len[key] = int(value)
    return seq_len

rule clean_countings_length:
    """
    Given the extension context, remove those reads that reported mappings greater than sequence length,
    the reference mapping didn't included N_s.
    """
    input:
        mutations_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}","{reference_seq}_{sample}.mut"),
        length_references = os.path.join(config['project_path'], config['data_folder'], config['trna_reference_original_len'])
    output:
        result_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}","{reference_seq}_{sample}_modified.mut"),
    run:
        import subprocess
        import os

        seq_len_dict = build_len_dict(input.length_references)
        if os.path.exists(input.mutations_file) and os.path.getsize(input.mutations_file) > 0 and wildcards.reference_seq in seq_len_dict:
            reference_seq_length = str(seq_len_dict[wildcards.reference_seq])
            call_args = ['awk', '-v', f'threshold={reference_seq_length}', '$4 < threshold {print $0}', input.mutations_file]
            with open(output.result_file, 'w') as output_file:
                subprocess.run(call_args, stdout=output_file, check=True)
        else:
            r = open(output.result_file, "x")
            r.close()

rule shapemapper_count:
    """
    $program_counter --in $outfolder/$out_namem.mut -c $outfolder/$out_namem.out --hist $outfolder/$out_namem.hist --separate_ambig_counts
    """
    input:
        mutations_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}","{reference_seq}_{sample}_modified.mut"),
        length_references = os.path.join(config['project_path'], config['data_folder'], config['trna_reference_original_len'])
    output:
        result_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}","{reference_seq}_{sample}.out"),
        #histogram_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}","{reference_seq}_{sample}.hist"),
    run:
        import subprocess
        import os
        
        seq_len_dict = build_len_dict(input.length_references)
        if os.path.exists(input.mutations_file) and os.path.getsize(input.mutations_file) > 0 and wildcards.reference_seq in seq_len_dict:
            reference_seq_length = str(seq_len_dict[wildcards.reference_seq])
            call_args = ['shapemapper_mutation_counter',
                         '-i', input.mutations_file,
                         '-n', reference_seq_length,
                         '-c', output.result_file,
                         '-w',
                         '--separate_ambig_counts']
            subprocess.run(call_args)
        else:
            r = open(output.result_file, "x")
            r.close()


def get_input_files_control_mutation(wildcards):
    import re

    file_list = []
    folder = wildcards.treatment
    sample = wildcards.sample
    ref_seq = wildcards.reference_seq
    cell = wildcards.cell
    path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], folder) 
    files_in_folder = glob_wildcards(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], folder, "{ref}.out")).ref
    for f in files_in_folder:
        # Asp-GTC-1(Asp-GTC)_C_i1.out
        name = f.split('_')
        #C
        f0 = name[0]
        part1 = name[1]
        #_i2.out
        part2 = name[2]
        part2 = re.sub(r'.out','',part2)
        #C_i1
        namem = str(part1)+"_"+str(part2)
        selected_sample = []
        if folder == "NAI":
            selected_sample = NAI
        elif folder == "DMS":
            selected_sample = DMS
        elif folder == "DMSO":
            selected_sample = DMSO

        selected_cell = SAMPLES_DICT[namem]['Cell']
        if namem in selected_sample and f0 == ref_seq and cell == selected_cell:
            concatenate = os.path.join(str(path_all),str(f)+".out")
            file_list.append(concatenate)
    return file_list

rule normalize_mutation_count_from_replicates:
    input:
        mut_file = get_input_files_control_mutation
    output:
        result_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder_summary'], "{treatment}", "{reference_seq}_{sample}_{cell}_summarised.out"),
    run:
        import pandas as pd
        dfs = []
        for replicate_file in input.mut_file:
            if os.path.getsize(replicate_file) == 0:
                continue
            with open(replicate_file, 'r') as f:
                    tto_control = str(f).split('_')[1].replace('.out','')
                    df = pd.read_csv(f, sep="\s+")
                    df['TTO'] = str(tto_control)
            if not df.empty:
                dfs.append(df)
        if dfs:
            complete = pd.concat(dfs, axis=0).reset_index(drop=False)
            complete_sum = complete.groupby('index').mean(numeric_only=True).reset_index()
            complete_sum = complete_sum.drop(columns='index')
            if not complete_sum.empty:
                complete_sum.to_csv(output.result_file, index=False, na_rep='nan', sep="\t")
            else:
                m = open(output.result_file, 'x')
                m.close()
        else:
            m = open(output.result_file, 'x')
            m.close()

rule shapemapper_reactivity_profiles:
    """
    ./make_reactivity_profiles.py --fa ${sequences_fasta}/$references.fa --counts $outfolder/$out_namem2.out $outfolder/$out_namem.out --out $outfolder/$references.reactivities
    Compare control with TTOs and get reactivities
    Control is harcoded to get the comparison reference: folder DMSO and file T3_i1
    """
    input:
        target_fasta = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "{reference_seq}.fa" ),
        control = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "DMSO", "{reference_seq}_{sample1}.out"),
        shape = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}", "{reference_seq}_{sample2}.out"),
    output:
        reactivities = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_reactivities.txt"),
    params: 
        shape_treatment = "{treatment}"
    run:
        import subprocess
        import os
        
        if os.path.exists(input.control) and os.path.getsize(input.control) > 0 and os.path.exists(input.shape) and os.path.getsize(input.shape) > 0:
            if params.shape_treatment == "DMS":
                call_args = ['make_reactivity_profiles.py',
                             '--fa', input.target_fasta,
                             '--counts', input.shape, input.control,
                             '--dms',
                             '--out', output.reactivities]
                subprocess.run(call_args)
            elif params.shape_treatment == "NAI":
                call_args = ['make_reactivity_profiles.py',
                             '--fa', input.target_fasta,
                             '--counts', input.shape, input.control,
                             '--out', output.reactivities]
                subprocess.run(call_args)
            else:
                sys.exit(1)
        else:
            r = open(output.reactivities, "x")
            r.close()

rule shapemapper_reactivity_profiles_summary:
    # Comparisons summarised vs. summarised
    input:
        target_fasta = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "{reference_seq}.fa" ),
        control = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder_summary'], "DMSO", "{reference_seq}_{sample1}_{cell}_summarised.out"),
        shape = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder_summary'], "{treatment}", "{reference_seq}_{sample2}_{cell}_summarised.out"),
    output:
        reactivities = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}_{cell}-{sample2}_{cell}_reactivities.txt"),
    params: 
        shape_treatment = "{treatment}"
    run:
        import subprocess
        import os
        
        if os.path.exists(input.control) and os.path.getsize(input.control) > 0 and os.path.exists(input.shape) and os.path.getsize(input.shape) > 0:
            if params.shape_treatment == "DMS":
                call_args = ['make_reactivity_profiles_rnahaydn.py',
                             '--fa', input.target_fasta,
                             '--counts', input.shape, input.control,
                             '--dms',
                             '--out', output.reactivities]
                subprocess.run(call_args)
            elif params.shape_treatment == "NAI":
                call_args = ['make_reactivity_profiles_rnahaydn.py',
                             '--fa', input.target_fasta,
                             '--counts', input.shape, input.control,
                             '--out', output.reactivities]
                subprocess.run(call_args)
            else:
                sys.exit(1)
        else:
            r = open(output.reactivities, "x")
            r.close()

use rule shapemapper_reactivity_profiles_summary as shapemapper_reactivity_profiles_summary_ind with:
    # Comparisons summarised vs. individual
    input:
        target_fasta = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "{reference_seq}.fa" ),
        control = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder_summary'], "DMSO", "{reference_seq}_{sample1}_{cell}_summarised.out"),
        shape = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['mutation_folder'], "{treatment}", "{reference_seq}_{sample2}.out"),
    output:
        reactivities = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}_{cell}-{sample2}_{cell}_reactivities.txt"),
        #reactivities = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}_reactivities.txt"),

import os
from snakemake.io import glob_wildcards

def get_input_files_reactivities(wildcards):
    import re

    react_folder = config['reactivity_folder']
    folder = wildcards.treatment
    sample2 = wildcards.sample2

    path_all = os.path.join(
        config['project_path'],
        config['output_folder'],
        "shapemapper2",
        react_folder,
        folder
    )

    # Find all files ending with _reactivities.txt
    files_in_folder = glob_wildcards(
        os.path.join(path_all, "{ref}_reactivities.txt")
    ).ref

    # Define allowed samples per treatment
    if folder == "NAI":
        selected_sample = NAI
    elif folder == "DMS":
        selected_sample = DMS
    else:
        selected_sample = []

    file_list = []

    for f in files_in_folder:
        # Skip unwanted variants
        if "_norm_reactivities" in f or "OLDRE" in f:
            continue

        # Extract the last part after the last dash
        sample_in_file = f.split("-")[-1]

        print(f"DEBUG: file={f}, extracted sample={sample_in_file}, sample2={sample2}")

        if sample_in_file == sample2 and sample2 in selected_sample:
            fullpath = os.path.join(path_all, f + "_reactivities.txt")
            file_list.append(fullpath)

    print("Selected files:", file_list)
    return file_list

from snakemake.io import glob_wildcards

def get_input_files_reactivities(wildcards):
    import re

    react_folder = config['reactivity_folder']
    folder = wildcards.treatment
    sample2 = wildcards.sample2

    path_all = os.path.join(
        config['project_path'],
        config['output_folder'],
        "shapemapper2",
        react_folder,
        folder
    )

    # Find all files ending with _reactivities.txt
    files_in_folder = glob_wildcards(
        os.path.join(path_all, "{ref}_reactivities.txt")
    ).ref

    # Define allowed samples per treatment
    if folder == "NAI":
        selected_sample = NAI
    elif folder == "DMS":
        selected_sample = DMS
    else:
        selected_sample = []

    file_list = []

    for f in files_in_folder:
        # Skip unwanted variants
        if "_norm_reactivities" in f or "OLDRE" in f:
            continue

        # Extract the last part after the last dash
        sample_in_file = f.split("-")[-1]

        #print(f"DEBUG: file={f}, extracted sample={sample_in_file}, sample2={sample2}")

        if sample_in_file == sample2 and sample2 in selected_sample:
            fullpath = os.path.join(path_all, f + "_reactivities.txt")
            file_list.append(fullpath)

    #print("Selected files:", file_list)
    return file_list

# Create concatenation of reactivity tables based on shape-TTO.
rule concatenate_reactivities:
    input:
        files = get_input_files_reactivities
    output:
        final = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}", "{sample2}_concatenated_react.txt"),
    run:
        with open(output.final, 'w') as outfile:
            lines = []
            head = False
            for file in input.files:
                with open(file) as infile:
                    for line in infile:
                        # Add header only once when reading all files
                        if line.startswith('Nucleotide') and not head:
                            head = True 
                            lines.append(line.strip()) 
                        elif not line.startswith('Nucleotide'):
                            lines.append(line.strip())
            outfile.write('\n'.join(lines))

def get_input_files_reactivities_sum(wildcards):
    import re
    file_list = []
    react_folder = config['reactivity_folder_summarised']
    folder = wildcards.treatment
    sample2 = wildcards.sample2
    cell = wildcards.cell
    path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", react_folder, folder) 
    files_in_folder = glob_wildcards(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", react_folder, folder, "{ref}_reactivities.txt")).ref
    for f in files_in_folder:
        #Leu-AAG-1-1(Leu-AAG)_DMSO_HAP1-DMS_HAP1 _reactivities.txt
        #Asp-GTC-1(Asp-GTC)[H]_T3_i1-T3_i2_reactivities.txt
        name = f.split('_')
        part1 = name[2] #HAP1-DMS
        part1 = re.sub(r'.*-','',part1) #DMS
        part2 = name[3] #HAP1
        namem = str(part1)+"_"+str(part2)
        selected_sample = []
        if folder == "NAI":
            selected_sample = NAI
        elif folder == "DMS":
            selected_sample = DMS
        #NAI KO NAI_KO NAI
        if part1 == sample2 and cell == part2:
            if not ("_norm_reactivities.txt" or "OLDRE") in f:
                concatenate = os.path.join(str(path_all),str(f)+"_reactivities.txt")
                file_list.append(concatenate)
    return file_list

def get_input_files_reactivities_sum_ind(wildcards):
    import re
    file_list = []
    react_folder = config['reactivity_folder_summarised_ind']
    folder = wildcards.treatment
    sample2 = wildcards.sample2
    cell = wildcards.cell
    path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", react_folder, folder) 
    files_in_folder = glob_wildcards(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", react_folder, folder, "{ref}_reactivities.txt")).ref
    for f in files_in_folder:
        #Leu-AAG-1-1(Leu-AAG)_DMSO_HAP1-DMS_HAP1 _reactivities.txt
        #Asp-GTC-1(Asp-GTC)[H]_T3_i1-T3_i2_reactivities.txt
        name = f.split('_')
        part1 = name[2] #HAP1-DMS
        part1 = re.sub(r'.*-','',part1) #DMS
        part2 = name[3] #HAP1
        namem = str(part1)+"_"+str(part2)
        selected_sample = []
        if folder == "NAI":
            selected_sample = NAI
        elif folder == "DMS":
            selected_sample = DMS
        #NAI KO NAI_KO NAI
        #print(cell, part2, namem, sample2)
        if namem == sample2:
            if not ("_norm_reactivities.txt" or "OLDRE") in f:
                concatenate = os.path.join(str(path_all),str(f)+"_reactivities.txt")
                file_list.append(concatenate)
    return file_list

use rule concatenate_reactivities as concatenate_reactivities_sum with:
    input:
        files = get_input_files_reactivities_sum
    output:
        final = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}", "{sample2}_{cell}_concatenated_react.txt"),

use rule concatenate_reactivities as concatenate_reactivities_sum_ind with:
    input:
        files = get_input_files_reactivities_sum_ind
    output:
        final = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}", "{sample2}_{cell}_concatenated_react.txt"),

rule get_normalization_values:
    input:
        concatenated_file = ancient(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}", "{sample2}_concatenated_react.txt")),
    output:
        values_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}", "{sample2}_react_normalization_values.txt"),
    params:
        shape_treatment = "{treatment}"
    run:
        import subprocess
        import os
        if params.shape_treatment == "DMS":
            dms_tto = 1
        else:
            dms_tto = 0
        call_args = ['infer_normalization_scores.py',
                     input.concatenated_file,
                     str(dms_tto),
                     output.values_file]
        subprocess.run(call_args)

use rule get_normalization_values as get_normalization_values_summ with:
    input:
        concatenated_file = ancient(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}", "{sample2}_concatenated_react.txt")),
    output:
        values_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}", "{sample2}_react_normalization_values.txt"),
    params:
        shape_treatment = "{treatment}"

use rule get_normalization_values as get_normalization_values_summ_ind with:
    input:
        concatenated_file = ancient(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}", "{sample2}_concatenated_react.txt")),
    output:
        values_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}", "{sample2}_react_normalization_values.txt"),
    params:
        shape_treatment = "{treatment}"

rule normalize_reactivity:
    input:
        reactivities = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_reactivities.txt"),
        values_file = ancient(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}", "{sample2}_react_normalization_values.txt")),
    output:
        reactivities_normalized = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_react_norm.txt"),
    run:
        import pickle
        import sys
        import pandas as pd
        import numpy as np

        reactivities = input.reactivities
        normalization = input.values_file

        if os.path.exists(reactivities) and os.path.getsize(reactivities) > 0 and os.path.exists(input.values_file) and os.path.getsize(input.values_file) > 0:
            with open(normalization, 'rb') as fp:
                values_normalization = pickle.load(fp)
                df = pd.read_csv(reactivities, engine="python",delim_whitespace=True)
                # Reeplace whitespace to nan.
                df2 = df.replace(r'^\s*$', np.nan, regex=True)
                df2['Norm_profile'] = df2.apply(lambda row: row['Reactivity_profile'] / values_normalization[row['Sequence']] if pd.notna(row['Reactivity_profile']) else "nan", axis=1)
                df2['Norm_stderr'] = df2.apply(lambda row: row['Std_err'] / values_normalization[row['Sequence']] if pd.notna(row['Std_err']) else "nan", axis=1)
                #df2['Norm_profile'] = df2.apply(lambda row: row['HQ_profile'] / values_normalization[row['Sequence']] if pd.notna(row['HQ_profile']) else "nan", axis=1)
                #df2['Norm_stderr'] = df2.apply(lambda row: row['HQ_stderr'] / values_normalization[row['Sequence']] if pd.notna(row['HQ_stderr']) else "nan", axis=1)
                df2.to_csv(output.reactivities_normalized, index=False, sep="\t", na_rep='nan')
        else:
            r = open(output.reactivities_normalized, "x")
            r.close()

use rule normalize_reactivity as normalize_reactivity_summ with:
    input:
        reactivities = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}-{sample2}_reactivities.txt"),
        values_file = ancient(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}", "{sample2}_react_normalization_values.txt")),
    output:
        reactivities_normalized = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}-{sample2}_react_norm.txt"),

use rule normalize_reactivity as normalize_reactivity_ind with:
    input:
        reactivities = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}_reactivities.txt"),
        values_file = ancient(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}", "{sample2}_react_normalization_values.txt")),
    output:
        reactivities_normalized = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}_react_norm.txt"),


rule shapemapper_generate_tables:
    """
    ./tab_to_shape.py --infile $outfolder/${references}_norm.reactivities --shape $outfolder/$references.shape --map $outfolder/$references.map --varna $outfolder/$references.varna
    """
    input:
        reactivities_normalized = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_react_norm.txt"),
    output:
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        map_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.map"),
        varna_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.varna"),
    run:
        import subprocess
        import os
        
        if os.path.exists(input.reactivities_normalized) and os.path.getsize(input.reactivities_normalized) > 0:
            call_args = ['tab_to_shape.py',
                         '--infile', input.reactivities_normalized,
                         '--shape', output.shape_file,
                         '--map', output.map_file,
                         '--varna', output.varna_file]
            subprocess.run(call_args)
        else:
            s = open(output.shape_file, "x")
            m = open(output.map_file, "x")
            v = open(output.varna_file, "x")
            s.close()
            m.close()
            v.close()

use rule shapemapper_generate_tables as shapemapper_generate_tables_summ with:
    input:
        reactivities_normalized = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}-{sample2}_react_norm.txt"),
    output:
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        map_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}-{sample2}.map"),
        varna_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}-{sample2}.varna"),

use rule shapemapper_generate_tables as shapemapper_generate_tables_summ_ind with:
    input:
        reactivities_normalized = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}_react_norm.txt"),
    output:
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        map_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}.map"),
        varna_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}.varna"),

rule shapemapper_generate_plots:
    """
    ./render_figures.py --infile $outfolder/${references}_norm.reactivities --plot $outfolder/${references}.pdf --hist $outfolder/${references}.hist --title $references
    """
    input:
        reactivities_normalized = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_react_norm.txt"),
    output:
        plot_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['plots_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_shapemapper.pdf"),
        stats_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['plots_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_stats.pdf"),
    params:
        shape_treatment = "{treatment}"
    run:
        import subprocess
        import os
        
        threshold_coverage = 1000
        if os.path.exists(input.reactivities_normalized) and os.path.getsize(input.reactivities_normalized) > 0:
            title = {wildcards.reference_seq}
            if params.shape_treatment == "NAI":
                call_args = ['render_figures.py',
                             '--infile', input.reactivities_normalized,
                             '--plot', output.plot_file,
                             '--hist', output.stats_file,
                             '--mindepth', threshold_coverage,
                             '--title', str(title)]
                subprocess.run(call_args)
            elif params.shape_treatment == "DMS":
                call_args = ['render_figures.py',
                             '--infile', input.reactivities_normalized,
                             '--plot', output.plot_file,
                             '--hist', output.stats_file,
                             '--mindepth', threshold_coverage,
                             '--dms',
                             '--title', str(title)]
                subprocess.run(call_args)
        else:
            s = open(output.plot_file, "x")
            l = open(output.stats_file, "x")
            s.close()
            l.close()

rule shapemapper_plot_ss:
    input:
        fasta = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "{reference_seq}.fa" ),
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        varna_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.varna"),
    output:
        plot_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['plots_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.svg"),
    params:
        output_path = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['plots_folder'], "{treatment}") 
    run:
        import subprocess
        import os
        
        if os.path.exists(input.varna_file) and os.path.getsize(input.varna_file) > 0:
            call_args = ['overlap_ss_shape_V2.py',
                         input.fasta,
                         input.varna_file,
                         input.shape_file,
                         wildcards.sample1,
                         wildcards.sample2,
                         params.output_path]
            subprocess.run(call_args)
        else:
            s = open(output.plot_file, "x")
            s.close()

rule svg_to_pdf:
    input:
        svg = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['plots_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.svg"),
    output:
        pdf = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['plots_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.pdf"),
    shell:
        """
        convert {input.svg} {output.pdf}
        """



rule trnas_to_table:
    """
    From the fasta that contains all tRNAs from tRNADB, generate a 
    reference table.
    TODO: The fasta file hsa_all_trnas.fa is downloaded manually
    """
    input:
        all_trna = os.path.join(config['project_path'], config['data_folder'], "hsa_all_trnas.fa"),
    output:
        table = os.path.join(config['project_path'], config['data_folder'], "hsa_all_trnas.fa.table"),
    run:
        import subprocess
        import os

        call_args = ['create_table_from_fasta_str.py',
                     input.all_trna,
                     output.table]
        subprocess.run(call_args)

rule create_natural_structures:
    input:
        fasta_trna = os.path.join(config['project_path'], config['data_folder'], config['trna_reference_original']),
        # Files from human tRNAs
        # TODO: The manual source of that is tRNADB <03-11-23, cavelandiah> -
        # (http://trnadb.bioinf.uni-leipzig.de)
        nuclear_hsa = os.path.join(config['project_path'], config['data_folder'], "hsa_all_nt.fa"),
        mt_hsa = os.path.join(config['project_path'], config['data_folder'], "hsa_mt.fa"),
        table = os.path.join(config['project_path'], config['data_folder'], "hsa_all_trnas.fa.table"),
    output:
        natural_mod_table = os.path.join(config['data_folderc'], config['natural_structures_table'])
    run:
        import subprocess
        import os

        call_args = ['get_relation_most_expressed_natural_trnas.py',
                     input.fasta_trna,
                     input.nuclear_hsa,
                     input.mt_hsa,
                     input.table,
                     output.natural_mod_table]
        subprocess.run(call_args)

rule shapemapper_experiments_folding:
    input:
        fasta = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "{reference_seq}.fa" ),
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        natural_mod_table = os.path.join(config['data_folderc'], config['natural_structures_table'])
    output:
        text_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_temp.txt"),
    params:
        output_path = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "{treatment}") 
    run:
        import subprocess
        import os
        
        if os.path.exists(input.shape_file) and os.path.getsize(input.shape_file) > 0:
            call_args = ['overlap_table_shape.py',
                         input.fasta,
                         input.shape_file,
                         wildcards.sample1,
                         wildcards.sample2,
                         params.output_path,
                         wildcards.treatment,
                         input.natural_mod_table]
            subprocess.run(call_args)
        else:
            m = open(output.text_file, "x")
            m.close()

use rule shapemapper_experiments_folding as shapemapper_experiments_distances_summ with:
    input:
        fasta = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "{reference_seq}.fa" ),
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        natural_mod_table = os.path.join(config['data_folderc'], config['natural_structures_table'])
    output:
        text_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder_sum'], "{treatment}","{reference_seq}_{sample1}-{sample2}_temp.txt"),
    params:
        output_path = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder_sum'], "{treatment}") 

use rule shapemapper_experiments_folding as shapemapper_experiments_distances_summ_id with:
    input:
        fasta = os.path.join(config['project_path'], config['data_folder'], config['split_reference_folder'], "{reference_seq}.fa" ),
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        natural_mod_table = os.path.join(config['data_folderc'], config['natural_structures_table'])
    output:
        text_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder_sum_ind'], "{treatment}","{reference_seq}_{sample1}-{sample2}_temp.txt"),
    params:
        output_path = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder_sum_ind'], "{treatment}") 


rule shapemapper_experiments_distances:
    input:
        shape_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}.shape"),
        text_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_temp.txt"),
    output:
        text = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "{treatment}","{reference_seq}_{sample1}-{sample2}_temp_distances.txt")
    run:
        import subprocess
        import os
        
        if os.path.exists(input.shape_file) and os.path.getsize(input.shape_file) > 0:
            if os.path.exists(input.text_file) and os.path.getsize(input.text_file) > 0:
                call_args = ['calculate_pair_probsV2.py',
                             input.text_file,
                             input.shape_file,
                             output.text]
                #with open(output.text, "w+") as f:
                subprocess.run(call_args)
            else:
                m = open(output.text, "x")
                m.close()
        else:
            m = open(output.text, "x")
            m.close()

def get_input_files(wildcards):
    """
    Concatenate all foldings into one file. This will include DMSO, Natural, and
    SHAPE treatements
    """
    file_list = []
    folder = wildcards.treatment
    reference = wildcards.reference_seq
    path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], folder)
    files_in_folder = glob_wildcards(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], folder, "{ref}_temp.txt")).ref
    for f in files_in_folder:
        name = f.split('_')
        # Here I can subset the concatenation based on name rules
        #Ser-GCT-6-1(Ser-GCT)_T1_i1-T2_i3_temp.txt
        #namem = str(name[0])+"_"+str(name[1])+"_"+str(name[2])
        #modified = reference+"_T3_i1-T3" 
        namem = str(name[0])
        if namem == reference:
            concatenate = os.path.join(str(path_all),str(f)+"_temp.txt")
            file_list.append(concatenate)
    #print(file_list)
    return file_list
 
rule concat_sort_unique:
    input:
        files = get_input_files
    output:
        final = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "{treatment}", "{reference_seq}_tab.txt"),
    run:
        with open(output.final, 'w') as outfile:
            lines = set()
            for file in input.files:
                with open(file) as infile:
                    for line in infile:
                        lines.add(line.strip())
            sorted_lines = sorted(lines)
            outfile.write('\n'.join(sorted_lines))

def get_input_files_dist(wildcards):
    file_list = []
    folder = wildcards.treatment
    reference = wildcards.reference_seq
    path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], folder)
    files_in_folder = glob_wildcards(os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], folder, "{ref}_temp_distances.txt")).ref
    for f in files_in_folder:
        name = f.split('_')
        #namem = str(name[0])+"_"+str(name[1])+"_"+str(name[2])
        #modified = reference+"_T3_i1-T3" 
        namem = str(name[0])
        if namem == reference:
            concatenate = os.path.join(str(path_all),str(f)+"_temp_distances.txt")
            file_list.append(concatenate)
    return file_list

#INPUT_FILE2 = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "{treatment}","{reference_seq}_T3_i1-{sample2}_temp_distances.txt")
#SAMPLE2=list(set(glob_wildcards(INPUT_FILE2).sample2))

rule concat_sort_unique_distances:
    input:
        files = get_input_files_dist
    output:
        final = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "{treatment}", "{reference_seq}_dist.txt"),
    run:
        with open(output.final, 'w') as outfile:
            lines = set()
            for file in input.files:
                with open(file) as infile:
                    for line in infile:
                        lines.add(line.strip())
            sorted_lines = sorted(lines)
            outfile.write('\n'.join(sorted_lines))

## Creation of plots 

def get_input_files_plots(wildcards):
    """
    Concatenate all reactivities from a specific reference into one file. This will 
    include DMSO and SHAPE treatements. Selection should 
    be cell specific. Based on the control, subset the corresponding tretatments.
    """
    import sys
    import re

    file_list = []
    target_treatments = ""
    reference = wildcards.reference_seq 
    reference_controls = wildcards.reference_control
    # Based on the reference control, get the cell and then build the Valid
    # Get SHAPE ttos
    control_cell = SAMPLES_DICT[reference_controls]['Cell'] 
    folders = SHAPE_NO_CONTROL
    target_treatments = SAMPLE_BY_CELL[SAMPLE_BY_CELL['Cell'] == control_cell]['Sample_id'].tolist()

    for folder in folders:
        path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder'], folder)
        files_in_folder = glob_wildcards(os.path.join(path_all, "{ref}_react_norm.txt")).ref
        for f in files_in_folder:
            #Ala-AGC-1-1(Ala-AGC)_T1_i1-T2_i1_react_norm.txt
            basename = os.path.basename(f)
            name = basename.split('_')
            tto_name = str(name[2])+"_"+str(name[3])
            tto_name = re.sub(r'.*-','', tto_name)
            name[2] = re.sub(r'-.*','',name[2])
            namem = str(name[0])+"_"+str(name[1])+"_"+str(name[2])
            ref_control = str(reference)+"_"+str(reference_controls)
            # Concatenate only those with same comparison
            if namem == ref_control and tto_name in target_treatments[0]:
                concatenate = os.path.join(str(path_all),str(f)+"_react_norm.txt")
                file_list.append(concatenate)
    return file_list

rule get_input_merge_plot:
    input:
        files = get_input_files_plots
    output:
        final_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}-{reference_control}_data_reactivities.txt"),
    params:
        reference = lambda wc: wc.get("reference_seq")
    run:
        import sys
        import os
        import re
        import pandas as pd

        def get_treatment_name(experiment):
            """
            Get name experiment based on label
            """
            if experiment in DMS:
                return "DMS"
            elif experiment in NAI:
                return "NAI"
            else:
                sys.error("Not defined experiment")

        list_dfs = []
        for f in input.files:
            name = os.path.basename(f)
            name = name.split('_')
            namem = str(name[0])
            if namem == params.reference:
                part1 = name[2]
                part1 = re.sub(r'.*-','',part1)
                part2 = name[3]
                exp = str(part1)+"_"+str(part2)
                if (os.path.exists(f) and os.path.getsize(f) > 0):
                    df = pd.read_csv(f, sep='\t', header=0)
                    df['Exp'] = exp
                    shape_chem = get_treatment_name(exp)
                    df['TTOs'] = shape_chem
                    df['trna'] = params.reference
                    list_dfs.append(df)
        if list_dfs:
            final_df = pd.concat(list_dfs, axis=0, ignore_index=False)
            final_df.to_csv(output.final_file, sep='\t', header=False, index=False)
        else:
            m = open(output.final_file, 'x')
            m.close()

def get_input_files_plots_summ(wildcards):
    """
    Concatenate all reactivities from a specific reference into one file. This will 
    include DMSO and SHAPE treatements. Selection should 
    be cell specific. Based on the control, subset the corresponding tretatments.
    Leu-AAG-1-1(Leu-AAG)_DMSO_KO-NAI_KO_react_norm.txt
    """
    import sys
    import re

    file_list = []
    target_treatments = ""
    reference_controls = wildcards.reference_control
    cell = wildcards.cell
    folders = SHAPE_NO_CONTROL
    for folder in folders: 
        path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised'], folder)
        files_in_folder = glob_wildcards(os.path.join(path_all, "{ref}_react_norm.txt")).ref
        for f in files_in_folder:
            #Leu-AAG-1-1(Leu-AAG)_DMSO_KO-NAI_KO_react_norm.txt
            basename = os.path.basename(f)
            name = basename.split('_')
            control_name = name[1]+"_"+name[2]
            control_name = re.sub(r'-.*','', control_name) #DMSO_KO
            tto_name = re.sub(r'.*-','', name[2])
            cell_name = name[3]
            # Concatenate only those with same comparison
            if tto_name == folder and cell_name == cell:
                concatenate = os.path.join(str(path_all),str(f)+"_react_norm.txt")
                file_list.append(concatenate)
    return file_list

rule get_input_merge_plot_summ:
    input:
        files = get_input_files_plots_summ
    output:
        final_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_summ_folder'], "{reference_seq}-{reference_control}_{cell}_data_reactivities.txt"),
    params:
        reference = lambda wc: wc.get("reference_seq")
    run:
        import sys
        import os
        import re
        import pandas as pd

        list_dfs = []
        for f in input.files:
            # Leu-AAG-1-1(Leu-AAG)_DMSO_KO-DMS_KO_react_norm.txt
            #Leu-AAG-1-1(Leu-AAG)_DMSO_HAP1-D_i2_HAP1_react_norm.txt
            name_file = os.path.basename(f)
            name = name_file.split('_')
            namem = str(name[0])
            # Match only anticodon family
            mref = re.sub(r'\(.*','',str(params.reference))
            mnamem = re.sub(r'\(.*\)','',namem)
            if mnamem == mref:
                exp = str(name[2]) # KO-DMS
                shape_chem = re.sub(r'.*-','',name[2])
                if (os.path.exists(f) and os.path.getsize(f) > 0):
                    df = pd.read_csv(f, sep='\t', header=0)
                    df['Exp'] = exp
                    df['TTOs'] = shape_chem
                    df['trna'] = mref
                    list_dfs.append(df)
        if list_dfs:
            final_df = pd.concat(list_dfs, axis=0, ignore_index=False)
            final_df.to_csv(output.final_file, sep='\t', header=False, index=False)
        else:
            m = open(output.final_file, 'x')
            m.close()

def get_input_files_plots_ind(wildcards):
    """
    Concatenate all reactivities from a specific reference into one file. This will 
    include DMSO and SHAPE treatements. Selection should 
    be cell specific. Based on the control, subset the corresponding tretatments.
    Leu-AAG-1-1(Leu-AAG)_DMSO_KO-NAI_KO_react_norm.txt
    """
    import sys
    import re

    file_list = []
    target_treatments = ""
    reference_controls = wildcards.reference_control
    cell = wildcards.cell
    folders = SHAPE_NO_CONTROL
    #target_treatments = SAMPLE_BY_CELL[SAMPLE_BY_CELL['Cell'] == control_cell]['Sample_id'].tolist()
    for folder in folders: 
        path_all = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['reactivity_folder_summarised_ind'], folder)
        files_in_folder = glob_wildcards(os.path.join(path_all, "{ref}_react_norm.txt")).ref
        for f in files_in_folder:
            #Leu-AAG-1-1(Leu-AAG)_DMSO_KO-NAI_KO_react_norm.txt
            #Leu-AAG-1-1(Leu-AAG)_DMSO_HAP1-D_i2_HAP1_react_norm.txt
            basename = os.path.basename(f)
            name = basename.split('_')
            control_name = name[1]+"_"+name[2]
            control_name = re.sub(r'-.*','', control_name) #DMSO_KO
            joinTTO = name[2]+"_"+name[3]
            tto_name = re.sub(r'.*-','', joinTTO)
            cell_name = name[4]
            # Concatenate only those with same comparison
            if cell_name == cell:
                concatenate = os.path.join(str(path_all),str(f)+"_react_norm.txt")
                file_list.append(concatenate)
    return file_list

rule get_input_merge_plot_ind:
    input:
        files = get_input_files_plots_ind
    output:
        final_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_ind_folder'], "{reference_seq}-{reference_control}_{cell}_data_reactivities.txt"),
    params:
        reference = lambda wc: wc.get("reference_seq")
    run:
        import sys
        import os
        import re
        import pandas as pd

        list_dfs = []
        for f in input.files:
            #Leu-AAG-1-1(Leu-AAG)_DMSO_KO-DMS_KO_react_norm.txt
            #Leu-AAG-1-1(Leu-AAG)_DMSO_HAP1-D_i2_HAP1_react_norm.txt
            name_file = os.path.basename(f)
            name = name_file.split('_')
            namem = str(name[0])
            # Match only anticodon family
            mref = re.sub(r'\(.*','',str(params.reference))
            mnamem = re.sub(r'\(.*\)','',namem)
            if mnamem == mref:
                exp = str(name[2])+"_"+str(name[3]) # HAP1-D_i2
                tto = re.sub(r'.*-','',exp) # Search for D_i2
                shape_chem = SAMPLES_DICT[tto]['Treatment'] 
                if (os.path.exists(f) and os.path.getsize(f) > 0):
                    df = pd.read_csv(f, sep='\t', header=0)
                    df['Exp'] = tto
                    df['TTOs'] = shape_chem
                    df['trna'] = mref
                    list_dfs.append(df)
        if list_dfs:
            final_df = pd.concat(list_dfs, axis=0, ignore_index=False)
            final_df.to_csv(output.final_file, sep='\t', header=False, index=False)
        else:
            m = open(output.final_file, 'x')
            m.close()
rule get_natural_ss:
    input:
        ss_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['tables_folder'], "NAI", "{reference_seq}_tab.txt"),
    output:
        final_ss_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}_ss.txt"),
    run:
        import pandas as pd
        import os

        if (os.path.exists(input.ss_file) and os.path.getsize(input.ss_file) > 0):
            df = pd.read_csv(input.ss_file, sep='\s+')
            # DMSO T1_i1 .......(((((...((((...(.(.((((((...((((((.....))))))....)))))).).)))))...))))) GUCUCUGUGGCGCAAUUGGUUAGCGCGUUCGGUUGUUAACCGUAAAGGUUGGUGGUUCGAGCCCACCCAGGAACGCCA -22.299999237060547 
            df.columns = ['TTO', 'EXP', "SS", "SEQ", "ENERGY_FOLD"]
            natural = df[df['TTO'] == "Natural"]
            # Extract the 'SS' column value
            ss_value = natural['SS'].iloc[0]
            filtered_ss = [char for char in ss_value]
            output_filt = '\n'.join(filtered_ss)
            with open(output.final_ss_file, 'w') as out:
                out.write(output_filt+"\n")
        else:
            m = open(output.final_ss_file, "x")
            m.close()

rule merge_plot:
    input:
        final_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}-{reference_control}_data_reactivities.txt"),
        final_ss_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}_ss.txt"),
        alignment_reference = os.path.join(config['data_folderc'], config['alignment_path'], config['trna_alignment_CM']),
        metadata = os.path.join(config["data_folder"], config["samples_file"]),
        natural_modifications = os.path.join(config['data_folderc'], "natural_modifications_"+str(config['species_name'])+"_relation_modomics.json"),
    output:
        final_plot = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}-{reference_control}_reactivities.pdf"),
    params:
        reference = lambda wc: wc.get("reference_seq"),
        control_tto = lambda wc: wc.get("reference_control"),
        mode = "all"
    run:
        import subprocess
        
        call_args = ['plot_reactivities_comparison_V4.py',
                     input.final_file,
                     input.final_ss_file,
                     params.reference,
                     input.alignment_reference,
                     params.control_tto,
                     input.metadata,
                     input.natural_modifications,
                     params.mode,
                     output.final_plot]
        subprocess.run(call_args)

use rule merge_plot as merge_plot_summary_ind with:
    input:
        final_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_ind_folder'], "{reference_seq}_{reference_control}_data_reactivities.txt"),
        final_ss_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}_ss.txt"),
        alignment_reference = os.path.join(config['data_folderc'], config['alignment_path'], config['trna_alignment_CM']),
        metadata = os.path.join(config["data_folder"], config["samples_file"]),
        natural_modifications = os.path.join(config['data_folderc'], "natural_modifications_"+str(config['species_name'])+"_relation_modomics.json"),
    output:
        final_plot = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_ind_folder'], "{reference_seq}_{reference_control}_reactivities.pdf"),
    params:
        reference = lambda wc: wc.get("reference_seq"),
        control_tto = lambda wc: wc.get("reference_control"),
        mode = "summary_ind"

use rule merge_plot as merge_plot_summary_sum with:
    input:
        final_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_summ_folder'], "{reference_seq}_{reference_control}_data_reactivities.txt"),
        final_ss_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}_ss.txt"),
        alignment_reference = os.path.join(config['data_folderc'], config['alignment_path'], config['trna_alignment_CM']),
        metadata = os.path.join(config["data_folder"], config["samples_file"]),
        natural_modifications = os.path.join(config['data_folderc'], "natural_modifications_"+str(config['species_name'])+"_relation_modomics.json"),
    output:
        final_plot = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_summ_folder'], "{reference_seq}_{reference_control}_reactivities.pdf"),
    params:
        reference = lambda wc: wc.get("reference_seq"),
        control_tto = lambda wc: wc.get("reference_control"),
        mode = "summary_all"

rule merge_plot_summary:
    input:
        final_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}-{reference_control}_data_reactivities.txt"),
        final_ss_file = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}_ss.txt"),
        alignment_reference = os.path.join(config['data_folderc'], config['alignment_path'], config['trna_alignment_CM']),
        metadata = os.path.join(config["data_folder"], config["samples_file"]),
        natural_modifications = os.path.join(config['data_folderc'], "natural_modifications_"+str(config['species_name'])+"_relation_modomics.json"),
    output:
        final_plot = os.path.join(config['project_path'], config['output_folder'], "shapemapper2", config['merged_plots_folder'], "{reference_seq}-{reference_control}_reactivities_summary.pdf"),
    params:
        reference = lambda wc: wc.get("reference_seq"),
        control_tto = lambda wc: wc.get("reference_control"),
        mode = "all"
    run:
        import subprocess
        
        call_args = ['plot_reactivities_comparison_V4_merged.py',
                     input.final_file,
                     input.final_ss_file,
                     params.reference,
                     input.alignment_reference,
                     params.control_tto,
                     input.metadata,
                     input.natural_modifications,
                     output.final_plot]
        subprocess.run(call_args)
