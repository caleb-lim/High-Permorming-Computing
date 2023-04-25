#! /usr/bin/env nextflow

version='1.0' // your submitted version should be 1.0
date='24/04/2023'  // update to the date that you last changed this file
author="Caleb Lim" // Change to your name

log.info """\
         PHYS4004 workflow assignment
         ============================
         version      : ${version} - ${date}
         author       : ${author}
         param.table  : ${params.table}
         --
         run as       : ${workflow.commandLine}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()

seeds =  Channel.from(5..95).filter{it%5==0}
ncores = Channel.of(1,2,4,7)
input_ch = seeds.combine(ncores).view()

println "params.image: ${params.image}"

process find {

        input:
        tuple(val(seed), val(cores)) from input_ch
        // the following images are constant across all versions of this process
        // so just use a 'static' or 'ad hoc' channel
        path(image) from file("${params.image}")
        path(bkg) from file("${params.bkg}")
        path(rms) from file("${params.rms}")

        output:
        file('*.csv') into files_ch

        // indicate that this process should be allocated a specific number of cores
        cpus "${cores}"
        
        script:
        """
        aegean ${image} --background=${bkg} --noise=${rms} --table=out.csv --seedclip=${seed} --cores=${cores}
        mv out_comp.csv table_${seed}_${cores}.csv
        """
}


process count {

        input:
        // The input should be all the files provided by the 'find' process
        // they are provided through the files_ch channel
        path(files) from files_ch.collect()

        output:
        file('results.csv') into counted_ch

	// don't use singularity for a bunch of bash commands, it's a waste
 	// (and also not all commands work in my container for some reason!)
	container = ''

        // Since we are using bash variables a lot and no nextflow variables
        // we use "shell" instead of "script" so that bash variables don't have
        // to be escaped
        shell:
        '''
        echo "seed,ncores,nsrc" > results.csv
	files=($(ls table*.csv))
        echo ${files[@]}
        for f in ${files[@]}; do
          seed_cores=($(echo ${f} | tr '_.' ' ' | awk '{print $2 " " $3}'))
          seed=${seed_cores[0]}
          cores=${seed_cores[1]}
          nsrc=$(echo "$(cat ${f}  | wc -l)-1" | bc -l)
          echo "${seed},${cores},${nsrc}" >> results.csv
        done
        '''
}

process plot {
        input:
        path(table) from counted_ch
        path(plot) from file("${params.plot}")

        output:
        file('*.png') into final_ch

	cpus 4

        script:
        // forloop
        // """
        // cores=\$(tail -n +2 results.csv | cut -d',' -f2 | sort | uniq)
        // echo "Cores unique: \${cores}\n"

        // for i in \${cores[@]}
        // do
        //     python ${plot} --infile ${table} --outfile core\${i}_${params.outfile} --cores \${i}
        // done
        // """

        //xargs
        """
        cores=\$(tail -n +2 results.csv | cut -d',' -f2 | sort | uniq)
        echo "Cores unique: \${cores}\n"

        echo "\${cores}" | xargs -n 1 -I {} python ${plot} --infile ${table} --outfile core{}_${params.outfile} --cores {}
        """
}
