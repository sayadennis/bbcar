import os
import sys
import re
import subprocess
import glob

class BBCarScript:
    """
    This class helps write scripts based on template scripts in the same directory (double check directory requirements).

    Example usage:
        from bbcar_script_generation import BBCarScript
        generator = BBCarScript(patid=1419, data_path="/path/to/fastq/data")
        generator.write_bwa(tmp_file="../bwa_tmp.sh", script_dir="../bwa_scripts")
        generator.write_gatk(temp_file="../gatk_tmp.sh", script_dir="../gatk_scripts")
    """
    def __init__(self, patid, data_path, A="b1042", p="genomics", mail="gannon.cottone@northwestern.edu", mailtype="END,FAIL"): # data_path must be absolute path. e.g. "/projects/p30007/Zexian/Alignment/BBCAR/RAW_data/1419/"
        # universal slurm settings
        self.patid = patid
        self.allocation = A
        self.partition = p
        self.mail = mail
        self.mailtype = mailtype
        self.data_path = data_path
    
    
    def write_bwa(self, tmp_file, script_dir, nodes=1, cores=12, mem="50G", time="12:00:00"): # tmp_file and script_dir can be relative or absolute path
        # slurm settings
        self.nodes = nodes
        self.cores = cores
        self.mem = mem
        self.time = time
        self.dnout = "/projects/p30007/gannon/bbcar/out"
        self.fnout = "bwa_{}.out".format(self.patid)
        self.jobname = "bwa_" + str(self.patid)

        flist = glob.glob(self.data_path + "/*")
        plist = [] # list all unique patterns aside from the replicates part from flist (e.g. "1419_S36_L005_")
        for fn in flist:
            if fn[-25:-11] in plist:
                continue
            else:
                plist.append(fn[-25:-11])
        pairs_dict = {} # will be dict with keys = matching pattern (S&L) // values = list of R1 and R2 filenames
        for i in range(len(plist)):
            fpair = []
            for fn in flist:
                if plist[i] in fn:
                    fpair.append(fn)
                else:
                    continue
            pairs_dict[plist[i]] = fpair # elements of flist that matches plist pattern (the two replicates)

        fn = "bwa_" + str(self.patid) + ".sh" # script filename

        with open(tmp_file, "r") as fsh_tmp: # set directory? 
            with open(os.path.join(script_dir, fn), "w") as fsh: # set directory? 
                for ln in fsh_tmp.readlines():
                    ln = ln.rstrip()
                    if ln == "[settings]":
                        fsh.write("#SBATCH -A {}\n".format(self.allocation))
                        fsh.write("#SBATCH -p {}\n".format(self.partition))
                        fsh.write("#SBATCH -t {}\n".format(self.time))
                        fsh.write("#SBATCH -N {}\n".format(self.nodes))
                        fsh.write("#SBATCH -n {}\n".format(self.cores))
                        fsh.write("#SBATCH --mem={}\n".format(self.mem))
                        fsh.write("#SBATCH --mail-user={}\n".format(self.mail))
                        fsh.write("#SBATCH --mail-type={}\n".format(self.mailtype))
                        fsh.write("#SBATCH --output={}\n".format(os.path.join(self.dnout, self.fnout)))
                        fsh.write("#SBATCH -J {}\n".format(self.jobname))
                    elif ln == "[bwa and picard]":
                        for i in range(len(pairs_dict)):
                            key = list(pairs_dict.keys())[i] # make sure that the order of files to be processed is ok with this
                            fsh.write("R1={}\n".format(pairs_dict[key][0])) # "$tissue_fq/{}/" removed –– glob() outputs entire path
                            fsh.write("R2={}\n".format(pairs_dict[key][1])) # "$tissue_fq/{}/" removed –– same as above
                            fsh.write("\n")
                            fsh.write("header=$(zcat $R1 | head -n 1)\n")
                            fsh.write("id=$(echo $header | head -n 1 | cut -f 3-4 -d\":\" | sed 's/@//' | sed 's/:/_/g')\n")
                            fsh.write("echo \"@RG\\tID:$id\\tSM:{}_germ\\tLB:library1\\tPL:ILLUMINA\\tPU:$id\"\n".format(self.patid))
                            fsh.write("\n")
                            fsh.write("bwa mem -M -t 24 -R $(echo \"@RG\\tID:$id\\tSM:{}_germ\\tLB:library1\\tPL:ILLUMINA\\tPU:$id\") $FA $R1 $R2 > {}_{}.sam\n".format(self.patid, self.patid, i+1))
                            fsh.write("\n")
                            fsh.write("$pic SortSam I={}_{}.sam O={}_{}_sorted.bam SORT_ORDER=coordinate\n\n".format(self.patid, i+1, self.patid, i+1))
                    elif ln == "[merge samples]":
                        fsh.write("$pic MergeSamFiles CREATE_INDEX=true ")
                        for i in range(len(pairs_dict)):
                            fsh.write("I={}_{}_sorted.bam ".format(self.patid, i+1))
                        fsh.write("O={}_final_sorted.bam USE_THREADING=true\n".format(self.patid))
                    elif ln == "[remove sam reps]":
                        for i in range(len(pairs_dict)):
                            fsh.write("rm -f {}_{}.sam\n".format(self.patid, i+1))
                    elif ln == "[remove bam reps]":
                        for i in range(len(pairs_dict)):
                            fsh.write("rm -f {}_{}_sorted.bam\n".format(self.patid, i+1))
                    elif ln == "[mark duplicates]":
                        fsh.write("$pic MarkDuplicates I={}_final_sorted.bam O={}_dup.bam M=$metrics/{}_tissue_reads.mdup.metrics.txt\n".format(self.patid, self.patid, self.patid))
                    elif ln == "[remove final sorted]":
                        fsh.write("rm -f {}_final_sorted.bam\n".format(self.patid))
                        fsh.write("rm -f {}_final_sorted.bai\n".format(self.patid))
                    elif ln == "[base recalibrate]":
                        fsh.write("gatk BaseRecalibrator -I {}_dup.bam -R $FA --known-sites $dbsnp --known-sites $gold1000Indel -O $rec_tables/{}_tissue_recal_data.table\n".format(self.patid, self.patid))
                    elif ln == "[index dup]":
                        fsh.write("samtools index {}_dup.bam\n".format(self.patid))
                    elif ln == "[ApplyBQSR]":
                        fsh.write("gatk ApplyBQSR -R $FA -I {}_dup.bam --bqsr-recal-file $rec_tables/{}_tissue_recal_data.table -L $interval -O {}t_bqsr.bam\n".format(self.patid, self.patid, self.patid))
                    elif ln == "[remove dups]":
                        fsh.write("rm -f {}_dup.bam\n".format(self.patid))
                        fsh.write("rm -f {}_dup.bam.bai\n".format(self.patid))
                    else:
                        fsh.write(ln + "\n")

    def write_gatk(self, tmp_file, script_dir, nodes=1, cores=24, mem="110G", time="00:30:00"): # tmp_file and script_dir can be relative or absolute path
        # slurm settings
        self.nodes = nodes
        self.cores = cores
        self.mem = mem
        self.time = time
        self.dnout = "/output/directory"
        self.fnout = "gatk_collectreads_{}.out".format(self.patid)
        self.jobname = "gatk_" + str(self.patid)
        
        fn = "gatk_" + str(self.patid) + ".sh" # script filename

        with open(tmp_file, "r") as fsh_tmp:
            with open(os.path.join(script_dir, fn), "w") as fsh:
                for ln in fsh_tmp.readlines():
                    ln = ln.rstrip()
                    if ln == "[settings]":
                        fsh.write("#SBATCH -A {}\n".format(self.allocation))
                        fsh.write("#SBATCH -p {}\n".format(self.partition))
                        fsh.write("#SBATCH -t {}\n".format(self.time))
                        fsh.write("#SBATCH -N {}\n".format(self.nodes))
                        fsh.write("#SBATCH -n {}\n".format(self.cores))
                        fsh.write("#SBATCH --mem={}\n".format(self.mem))
                        fsh.write("#SBATCH --mail-user={}\n".format(self.mail))
                        fsh.write("#SBATCH --mail-type={}\n".format(self.mailtype))
                        fsh.write("#SBATCH --output={}\n".format(os.path.join(self.dnout, self.fnout)))
                        fsh.write("#SBATCH -J {}\n".format(self.jobname))
                    elif "[patid]" in ln:
                        lsplit = ln.split("[patid]")
                        newline = ""
                        for i in range(len(lsplit)-1):
                            newline = newline + lsplit[i] + str(self.patid)
                        newline = newline + lsplit[-1] + "\n"
                        fsh.write(newline)
                    else:
                        fsh.write(ln + "\n")
