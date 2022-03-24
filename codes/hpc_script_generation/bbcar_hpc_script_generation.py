import os
import sys
import re
import subprocess
import glob

class BBCarScript:
    '''
    This class helps write scripts based on template scripts in the same directory (double check directory requirements).

    Example usage:
        import os
        import sys
        sys.path.append('bbcar/codes/hpc_script_generation')
        from bbcar_hpc_script_generation import BBCarScript
        for category in ['tissue', 'germline']:
            rawdn=f'/projects/b1131/saya/bbcar/raw/{category}'
            for patid in os.listdir(rawdn):
                if patid.isnumeric():
                    generator = BBCarScript(patid=int(patid), category=category, data_path=f'{rawdn}/{patid}')
                    generator.write_bwa(tmp_file='bbcar/codes/hpc_script_generation/bwa_tmp.sh', script_dir=f'bbcar/codes/somatic_mutations/02_processing/01_alignment/{category}')
                    # generator.write_gatk(temp_file='../gatk_tmp.sh', script_dir='../gatk_scripts')
    '''
    def __init__(self, patid, data_path, category, A='b1042', p='genomics', mail='sayarenedennis@northwestern.edu', mailtype='END,FAIL'):
        # data_path must be absolute path. e.g. '/projects/b1122/Zexian/Alignment/BBCAR/RAW_data/1419/'
        # here we set universal slurm settings
        self.patid = patid
        self.category = category
        self.allocation = A
        self.partition = p
        self.mail = mail
        self.mailtype = mailtype
        self.data_path = data_path
    
    
    def write_bwa(self, tmp_file, script_dir, mem='50G', time='12:00:00'): # tmp_file and script_dir can be relative or absolute path
        # slurm settings
        self.nodes = 1
        # self.cores = cores
        self.mem = mem
        self.time = time
        self.dnout = '/home/srd6051/bbcar/out'
        self.fnout = f'bwa_{self.patid}_{self.category}.out'
        self.jobname = f'bwa_{self.patid}_{self.category}'

        flist = glob.glob(self.data_path + '/*.fastq.gz')
        plist = [] # list all unique patterns aside from the replicates part from flist (e.g. '1419_S36_L005_')
        for fn in flist:
            pattern_wo_rep = fn[:re.search(r'_R[0-9]', fn).span()[0]]
            if pattern_wo_rep in plist:
                continue
            else:
                plist.append(pattern_wo_rep)
        pairs_dict = {} # will be dict with keys = matching pattern (S&L) // values = list of R1 and R2 filenames
        for i in range(len(plist)):
            fpair = []
            for fn in flist:
                if plist[i] in fn:
                    fpair.append(fn)
                else:
                    continue
            pairs_dict[plist[i]] = fpair # elements of flist that matches plist pattern (the two replicates)
        
        self.cores=len(pairs_dict) # we can parallelize up to the number of BAMs that are generated 

        fn = 'bwa_' + str(self.patid) + '.sh' # script filename

        with open(tmp_file, 'r') as fsh_tmp: # set directory? 
            with open(os.path.join(script_dir, fn), 'w') as fsh: # set directory? 
                for ln in fsh_tmp.readlines():
                    ln = ln.rstrip()
                    if ln == '[settings]':
                        fsh.write(f'#SBATCH -A {self.allocation}\n')
                        fsh.write(f'#SBATCH -p {self.partition}\n')
                        fsh.write(f'#SBATCH -t {self.time}\n')
                        fsh.write(f'#SBATCH -N {self.nodes}\n')
                        fsh.write(f'#SBATCH -n {self.cores}\n')
                        fsh.write(f'#SBATCH --mem={self.mem}\n')
                        fsh.write(f'#SBATCH --mail-user={self.mail}\n')
                        fsh.write(f'#SBATCH --mail-type={self.mailtype}\n')
                        fsh.write(f'#SBATCH --output={self.dnout}/{self.fnout}\n')
                        fsh.write(f'#SBATCH --job-name={self.jobname}\n')
                    elif ln == '[set output directories]':
                        fsh.write(f'metrics=\'/projects/b1131/saya/bbcar/01_alignment/{self.category}/metrics\'\n')
                        fsh.write(f'rec_tables=\'/projects/b1131/saya/bbcar/01_alignment/{self.category}/recal_tables\'\n')
                        fsh.write(f'interim=\'/projects/b1131/saya/bbcar/01_alignment/{self.category}/interim\'\n')
                        fsh.write(f'aligned=\'/projects/b1131/saya/bbcar/01_alignment/{self.category}/aligned\'\n')
                    elif ln == '[bwa and picard]':
                        for i in range(len(pairs_dict)):
                            key = list(pairs_dict.keys())[i] # make sure that the order of files to be processed is ok with this
                            fsh.write(f'R1={pairs_dict[key][0]}\n') # '$tissue_fq/{}/' removed glob() outputs entire path
                            fsh.write(f'R2={pairs_dict[key][1]}\n') # '$tissue_fq/{}/' removed same as above
                            fsh.write('\n')
                            fsh.write('header=$(zcat $R1 | head -n 1)\n')
                            fsh.write('id=$(echo $header | head -n 1 | cut -f 3-4 -d\':\' | sed \'s/@//\' | sed \'s/:/_/g\')\n')
                            fsh.write(f'echo \'@RG\\tID:$id\\tSM:{self.patid}_{self.category}\\tLB:library1\\tPL:ILLUMINA\\tPU:$id\'\n')
                            fsh.write('\n')
                            fsh.write(f'bwa mem -M -t 24 -R $(echo \'@RG\\tID:$id\\tSM:{self.patid}_{self.category}\\tLB:library1\\tPL:ILLUMINA\\tPU:$id\') $FA $R1 $R2 > $interim/{self.patid}_{i+1}.sam\n')
                            fsh.write('\n')
                            fsh.write(f'$pic SortSam I=$interim/{self.patid}_{i+1}.sam O=$interim/{self.patid}_{i+1}_sorted.bam SORT_ORDER=coordinate\n\n')
                    elif ln == '[merge samples]':
                        fsh.write('$pic MergeSamFiles CREATE_INDEX=true ')
                        for i in range(len(pairs_dict)):
                            fsh.write(f'I=$interim/{self.patid}_{i+1}_sorted.bam ')
                        fsh.write(f'O=$interim/{self.patid}_final_sorted.bam USE_THREADING=true\n')
                    elif ln == '[remove sam reps]':
                        for i in range(len(pairs_dict)):
                            fsh.write(f'rm -f $interim/{self.patid}_{i+1}.sam\n')
                    elif ln == '[remove bam reps]':
                        for i in range(len(pairs_dict)):
                            fsh.write(f'rm -f $interim/{self.patid}_{i+1}_sorted.bam\n')
                    elif ln == '[mark duplicates]':
                        fsh.write(f'$pic MarkDuplicates I=$interim/{self.patid}_final_sorted.bam O=$interim/{self.patid}_dup.bam M=$metrics/{self.patid}_{self.category}_reads.mdup.metrics.txt\n')
                    elif ln == '[remove final sorted]':
                        fsh.write(f'rm -f $interim/{self.patid}_final_sorted.bam\n')
                        fsh.write(f'rm -f $interim/{self.patid}_final_sorted.bai\n')
                    elif ln == '[base recalibrate]':
                        fsh.write(f'gatk BaseRecalibrator -I $interim/{self.patid}_dup.bam -R $FA --known-sites $dbsnp --known-sites $gold1000Indel -O $rec_tables/{self.patid}_{self.category}_recal_data.table\n')
                    elif ln == '[index dup]':
                        fsh.write(f'samtools index $interim/{self.patid}_dup.bam\n')
                    elif ln == '[ApplyBQSR]':
                        fsh.write(f'gatk ApplyBQSR -R $FA -I $interim/{self.patid}_dup.bam --bqsr-recal-file $rec_tables/{self.patid}_{self.category}_recal_data.table -L $interval -O $aligned/{self.patid}_bqsr.bam\n')
                    elif ln == '[remove dups]':
                        fsh.write(f'rm -f $interim/{self.patid}_dup.bam\n')
                        fsh.write(f'rm -f $interim/{self.patid}_dup.bam.bai\n')
                    else:
                        fsh.write(ln + '\n')

    def write_gatk(self, tmp_file, script_dir, nodes=1, cores=24, mem='110G', time='00:30:00'): # tmp_file and script_dir can be relative or absolute path
        # slurm settings
        self.nodes = nodes
        self.cores = cores
        self.mem = mem
        self.time = time
        self.dnout = '/output/directory'
        self.fnout = 'gatk_collectreads_{}.out'.format(self.patid)
        self.jobname = 'gatk_' + str(self.patid)
        
        fn = 'gatk_' + str(self.patid) + '.sh' # script filename

        with open(tmp_file, 'r') as fsh_tmp:
            with open(os.path.join(script_dir, fn), 'w') as fsh:
                for ln in fsh_tmp.readlines():
                    ln = ln.rstrip()
                    if ln == '[settings]':
                        fsh.write('#SBATCH -A {}\n'.format(self.allocation))
                        fsh.write('#SBATCH -p {}\n'.format(self.partition))
                        fsh.write('#SBATCH -t {}\n'.format(self.time))
                        fsh.write('#SBATCH -N {}\n'.format(self.nodes))
                        fsh.write('#SBATCH -n {}\n'.format(self.cores))
                        fsh.write('#SBATCH --mem={}\n'.format(self.mem))
                        fsh.write('#SBATCH --mail-user={}\n'.format(self.mail))
                        fsh.write('#SBATCH --mail-type={}\n'.format(self.mailtype))
                        fsh.write('#SBATCH --output={}\n'.format(os.path.join(self.dnout, self.fnout)))
                        fsh.write('#SBATCH -J {}\n'.format(self.jobname))
                    elif '[patid]' in ln:
                        lsplit = ln.split('[patid]')
                        newline = ''
                        for i in range(len(lsplit)-1):
                            newline = newline + lsplit[i] + str(self.patid)
                        newline = newline + lsplit[-1] + '\n'
                        fsh.write(newline)
                    else:
                        fsh.write(ln + '\n')
