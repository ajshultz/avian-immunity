#!/usr/bin/env python

#For use with Python 3.

import re
import sys
import os
import argparse
from subprocess import Popen,PIPE
from time import sleep
import datetime

#Extract Sample, SRA and genome info from config file. Sample data will be stored as a dictionary with sample ID as keys and a list of SRA accessions as values. Returns a dictionary of this info.
def extract_config(config_filename):
    print("Opening %s"%config_filename)
    config_file = open(config_filename,"r")
    
    #Create dictionaries
    sample_ncbi_dict = {}
    sample_bioproj_dict = {}
    sample_ref_dict = {}
    bioproj_list = []
    
    for line in config_file:
        if line[0] == "#":
            pass
        elif line == "\n":
            pass
        elif line[0:7] == 'Species':
            pass
        else:
            line=line.strip()
            line = line.split(",")
            
            if line[5] not in sample_ncbi_dict:
                sample_ncbi_dict[line[5]] = [line[7]]
            elif line[5] in sample_ncbi_dict:
                sample_ncbi_dict[line[5]].append(line[7])
            
            if line[5] not in sample_bioproj_dict:
                sample_bioproj_dict[line[5]] = line[4]
            
            if line[5] not in sample_ref_dict:
                sample_ref_dict[line[5]] = line[1]
                
            if line[4] not in bioproj_list:
                bioproj_list.append(line[4])
            
    config_file.close()
        
    #Return objects    
    return(sample_ncbi_dict,sample_bioproj_dict,sample_ref_dict,bioproj_list)


#Check if a directory exists. If not, create it.
def directory_create(test_dir):
    dir = os.path.dirname("%s/"%test_dir)
    if not os.path.exists(dir):
        os.makedirs(dir)

#Create generic slurm script
def script_create():
    slurm_script = '''#!/bin/bash\n#SBATCH -p {partition}\n#SBATCH -t {time}\n#SBATCH --mem {mem}\n#SBATCH -n {cores}\n#SBATCH -N {nodes}\n#SBATCH -J {jobid}\n#SBATCH -o {sp_dir}/logs/{jobid}_%j.out\n#SBATCH -e {sp_dir}/logs/{jobid}_%j.err\n\n{cmd}
    '''
    return(slurm_script)

    
#Submit filename to slurm with sbatch, returns job id number
def sbatch_submit(filename):
    proc = Popen('sbatch %s'%filename,shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr = proc.communicate()
    stdout = stdout.decode("utf-8","ignore")
    stdout = stdout.strip()
    stdout = stdout.strip('Submitted batch job ')
    return(stdout)


#Check job status of specific jobid: returns job status
def jobid_status(jobid,date):
    proc = Popen('sacct --format state --noheader -j %d -S %s'%(int(jobid),date),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr=proc.communicate()
    if proc.returncode != 0:
        raise Exception('Error running sacct: %s'%stderr)
    if stdout.strip() == '':
        return("No Status")
    lines = stdout.split()
    return(lines[0].decode("utf-8","ignore"))
    
    
#Check the status of all jobs, returns dictionary of jobid:status. Ignores jobs with ".ba+" and ".ex+" file extensions.
def all_jobs_status(date):
    proc = Popen('sacct --format jobid,state --noheader -S %s'%(date),shell=True,stdout=PIPE,stderr=PIPE)
    stdout,stderr=proc.communicate()
    stdout=stdout.decode("utf-8","ignore")
    stderr=stderr.decode("utf-8","ignore")
    if proc.returncode != 0:
        raise Exception('Error running sacct: %s'%stderr)
    if stdout.strip() == '':
        return("No Status")
    lines = stdout.split("\n")
    status_dict = {}
    for line in lines:
        line = line.split()
        if len(line) > 1:
            if not re.search('.ba+',line[0]) and not re.search('.ex+',line[0]):
                status_dict[line[0]]=line[1]
    return(status_dict)
    

def num_pend_run(job_id_list,date):
    count = 0
    status_dict = all_jobs_status(date)
    for job in job_id_list:
        if status_dict[job] == "PENDING" or status_dict[job] == "RUNNING":
            count += 1
    return(count)


#Create an sbatch file for a given set of SRAs and split into fastq files. Returns a list of new sbatch filenames
def ncbi_sra_download_sbatch(working_dir,sample_ncbi_dict,sample_bioproj_dict):
    slurm_script = script_create()
    sra_dl_sbatch_filenames = []
    path_to_sratools = "/n/home13/ashultz/sw/progs/sratoolkit.2.8.2-1-centos_linux64/bin/"

    for sample in sample_ncbi_dict.keys():
        for sra in sample_ncbi_dict[sample]:
            
            bioproj_dir = "%s/%s"%(working_dir,sample_bioproj_dict[sample])

            #First check if fastq file is already present (already downloaded). If it is, print statment and continue with next sample. 
            fastq_1_filename = '%s/fastq/%s_1.fastq.gz'%(bioproj_dir,sra)
            if os.path.isfile(fastq_1_filename):
                print('%s_1.fastq.gz already present, skipping'%(sra))
            else:
                print('Will download %s for sample %s and split'%(sra,sample))
    
                #Load modules and get versions for all programs used
                cmd_1 = 'wget -P %s/sra/ ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra'%(bioproj_dir,sra[0:3],sra[0:6],sra,sra)
                cmd_2 = r'%sfastq-dump --outdir %s/fastq --gzip --split-files %s/sra/%s.sra'%(path_to_sratools,bioproj_dir,bioproj_dir,sra)
            
                final_cmd = "%s\n\n%s"%(cmd_1,cmd_2)
    
        #Format sbatch script
                sra_script = slurm_script.format(partition="shared",time="0-12:00",mem="4000",cores="1",nodes="1",jobid="SRA",sp_dir=bioproj_dir,cmd=final_cmd)
                out_filename = "%s/scripts/01_sra_download_parse_%s.sbatch"%(bioproj_dir,sra)
                out_file = open(out_filename,"w")
                out_file.write(sra_script)
                out_file.close
                sra_dl_sbatch_filenames.append(out_filename)
    
    return(sra_dl_sbatch_filenames)


#Create sbatch files for running Kallisto given paired-end sequencing.
def  kallisto_paired_sbatch(working_dir,sample, sra_1_list, sra_2_list,bioproj,index_name):
	print('Will run Kallisto for sample %s in paired-end mode'%(sample))
	
	output_dir = '%s/%s/%s'%(working_dir,bioproj,sample)
	
	sra_combo_list = []
	for i in range(0,len(sra_1_list)):
	    sra_combo_list.append(sra_1_list[i])
	    sra_combo_list.append(sra_2_list[i])
	    
	sra_combo = " ".join(sra_combo_list)

	#Load necessary modules
	cmd_1 = 'module load gcc/4.8.2-fasrc01 kallisto/0.43.1-fasrc01'

    #Run Kallisto
	cmd_2 = 'kallisto quant -i ./indexes/%s -o %s %s'%(index_name,output_dir,sra_combo)

	cmd_list = [cmd_1,cmd_2]

	final_cmd = "\n\n".join(cmd_list)

	#Format sbatch script and write file
	slurm_script = script_create()
	genome_script = slurm_script.format(partition="shared",time="0-6:00",mem="4000",cores="1",nodes="1",jobid="Kallisto",sp_dir="%s/%s"%(working_dir,bioproj),cmd=final_cmd)

	out_filename = "%s/%s/scripts/02_kallisto_%s.sbatch"%(working_dir,bioproj,sample)
	out_file = open(out_filename,"w")
	out_file.write(genome_script)
	out_file.close

	return(out_filename)


#Create sbatch files for running Kallisto given single-end sequencing.
def  kallisto_single_sbatch(working_dir,sample, sra_1_list,bioproj,index_name):
	print('Will run Kallisto for sample %s in single-end mode'%(sample))
	
	output_dir = '%s/%s/%s'%(working_dir,bioproj,sample)
	    
	sra_combo = " ".join(sra_1_list)

	#Load necessary modules
	cmd_1 = 'module load gcc/4.8.2-fasrc01 kallisto/0.43.1-fasrc01'

    #Run Kallisto
	cmd_2 = 'kallisto quant -i ./indexes/%s -o %s --single -l 250 -s 50 %s'%(index_name,output_dir,sra_combo)

	cmd_list = [cmd_1,cmd_2]

	final_cmd = "\n\n".join(cmd_list)

	#Format sbatch script and write file
	slurm_script = script_create()
	genome_script = slurm_script.format(partition="shared",time="0-6:00",mem="4000",cores="1",nodes="1",jobid="Kallisto",sp_dir="%s/%s"%(working_dir,bioproj),cmd=final_cmd)

	out_filename = "%s/%s/scripts/02_kallisto_%s.sbatch"%(working_dir,bioproj,sample)
	out_file = open(out_filename,"w")
	out_file.write(genome_script)
	out_file.close

	return(out_filename)





def main():
    #Get config file from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file specifying samples")
    parser.add_argument("--clade", help="birds or mammals")
    args = parser.parse_args()
    config_filename = args.config
    clade = args.clade
    
    now = datetime.datetime.now()
    print('Staring work: %s'%now)
    start_date = now.strftime("%Y-%m-%d")
    
    #Open config file and get Sample, SRA and Genome attributes    
    sample_ncbi_dict,sample_bioproj_dict,sample_ref_dict,bioproj_list = extract_config(config_filename)
    
    #####Check if bioproject, sample, logs, and scripts directories exit, if not create them
    
    working_dir = "./%s"%clade
    
    print("\nOutput will be written to %s\n"%working_dir)
    
    for proj in bioproj_list:
        #Create bioproject directory
        proj_dir = "%s/%s"%(working_dir,proj)
        directory_create(proj_dir)
        
        #Create general directories
        logs_dir = "%s/logs"%(proj_dir)
        scripts_dir = "%s/scripts"%(proj_dir)
        sra_dir = "%s/sra"%(proj_dir)
        fastq_dir = "%s/fastq"%(proj_dir)
        
        directory_create(logs_dir)
        directory_create(scripts_dir)
        directory_create(sra_dir)
        directory_create(fastq_dir)
        
    

    #####Create sbatch files to download SRA files and use fastq-dump to split
    
    #Create sbatch files
    ncbi_sra_dl_sbatch_filenames = []
    
    #Create sbatch files for ncbi
    if len(sample_ncbi_dict) > 0:
        ncbi_sra_dl_sbatch_filenames = ncbi_sra_download_sbatch(working_dir,sample_ncbi_dict,sample_bioproj_dict)
    
    #Combine sbatch filenames into single object:
    sra_dl_sbatch_filenames = ncbi_sra_dl_sbatch_filenames
    
    '''      
    #Submit SRA read sbatch files, only allow 20 SRA jobs to run (or pend) at a time (set max_jobs)  
    max_jobs = 20
    sra_dl_jobids = []
    completed_jobids = {}
    job_count = 0
    sra_dl_jobid_filename_dict = {}
    #First submit up to the maximum number of jobs quickly
    if len(sra_dl_sbatch_filenames) > max_jobs:
        for i in range(0,max_jobs):
            job_id = sbatch_submit(sra_dl_sbatch_filenames[job_count])
            sra_dl_jobids.append(job_id)
            sra_dl_jobid_filename_dict[job_id] = sra_dl_sbatch_filenames[job_count]
            job_count += 1
            sleep(1)
    else:
        for i in range(0,len(sra_dl_sbatch_filenames)):
            job_id = sbatch_submit(sra_dl_sbatch_filenames[job_count])
            sra_dl_jobids.append(job_id)
            sra_dl_jobid_filename_dict[job_id] = sra_dl_sbatch_filenames[job_count]
            job_count += 1
            sleep(1)
    #Add an extra sleep to give sacct a chance to catch up
    sleep(20)
    #Then, enter while loop that will continue until the number of completed jobs matches the. number of sbatch files
    while len(completed_jobids) < len(sra_dl_sbatch_filenames):
        num_running = num_pend_run(sra_dl_jobids,start_date)
        while num_running < max_jobs and job_count < (len(sra_dl_sbatch_filenames)):
            job_id = sbatch_submit(sra_dl_sbatch_filenames[job_count])
            sra_dl_jobids.append(job_id)
            sra_dl_jobid_filename_dict[job_id] = sra_dl_sbatch_filenames[job_count]
            job_count += 1
            sleep(20)
            num_running = num_pend_run(sra_dl_jobids,start_date)
        job_statuses = all_jobs_status(start_date)
        for job in sra_dl_jobids:
            if job not in completed_jobids:
                try:
                    if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                        completed_jobids[job] = job_statuses[job]
                        print("Job %s completed"%job)
                except:
                    pass
        sleep(30)
    
    #After all jobs have finished, report which jobs failed
    for job in completed_jobids:
        if completed_jobids[job] != "COMPLETED":
            print("SRA download and parse job %s failed with code: %s for file %s"%(job,completed_jobids[job],sra_dl_jobid_filename_dict[job]))
    
    '''
    
    #####Final processing and mapping
    #Create output dir for that biosample within the bioproject directory   
    #Check if paired end or single end
    #Run Kallisto with parameters appropriate to single end or paired end data
    #Before submitting jobs, check to see that fastq files are there. If not, print info statement.
    mapping_jobids = []
    mapping_completed_jobids = {}
    mapping_jobid_filename_dict = {}
    
    for sample in sample_ncbi_dict:
        out_dir = "%s/%s/%s"%(working_dir,sample_bioproj_dict[sample],sample)
        
        sra_1_list = []
        sra_2_list = []
        for sra in sample_ncbi_dict[sample]:
            sra_1 = "%s/%s/fastq/%s_1.fastq.gz"%(working_dir,sample_bioproj_dict[sample],sra)
            sra_2 = "%s/%s/fastq/%s_2.fastq.gz"%(working_dir,sample_bioproj_dict[sample],sra)
            sra_1_list.append(sra_1)
            sra_2_list.append(sra_2)

        #Check if any missing first end files, if so, add to vec
        sra_1_missing = []
        for file in sra_1_list:
            if os.path.isfile(file):
                sra_1_missing.append(0)
            else:
                sra_1_missing.append(1)
                print("Oops, file %s is missing for sample %s"%(file,sample))
                
        #If no missing files, continue, check if missing second end files
        if sum(sra_1_missing) == 0:
            sra_2_missing = []
            for file in sra_2_list:
                if os.path.isfile(file):
                    sra_2_missing.append(0)
                else:
                    sra_2_missing.append(1)
            
            #If no missing second end files, means there is paired end data, run Kallisto in paired-end mode
            if sum(sra_2_missing) == 0:
                if os.path.isfile("%s/%s/%s/abundances.h5"%(working_dir,sample_bioproj_dict[sample],sample)):
                    print("Output already present for sample %s, skipping"%sample)
                else:
                    mapping_filename = kallisto_paired_sbatch(working_dir,sample, sra_1_list, sra_2_list,bioproj=sample_bioproj_dict[sample],index_name=sample_ref_dict[sample])
                    #mapping_jobid = sbatch_submit(mapping_filename)
                    #mapping_jobids.append(mapping_jobid)
                    #mapping_jobid_filename_dict[mapping_jobid] = mapping_filename
                    sleep(1)
            
            #If number of missing files is equal to the number of files, assume single end data, run Kallisto in single-end mode
            elif sum(sra_2_missing) == len(sra_2_list):
                if os.path.isfile("%s/%s/%s/abundances.h5"%(working_dir,sample_bioproj_dict[sample],sample)):
                    print("Output already present for sample %s, skipping"%sample)
                else:
                    mapping_filename = kallisto_single_sbatch(working_dir,sample, sra_1_list,bioproj=sample_bioproj_dict[sample],index_name=sample_ref_dict[sample])
                    #mapping_jobid = sbatch_submit(mapping_filename)
                    #mapping_jobids.append(mapping_jobid)
                    #mapping_jobid_filename_dict[mapping_jobid] = mapping_filename
                    sleep(1)  
    '''             
    sleep(60)
       
    while len(mapping_completed_jobids) < len(mapping_jobids):
        job_statuses = all_jobs_status(start_date)
        for job in mapping_jobids:
            if job not in mapping_completed_jobids:
                if job_statuses[job] != "PENDING" and job_statuses[job] != "RUNNING":
                    mapping_completed_jobids[job] = job_statuses[job]
                    print("Job %s completed"%job)
        sleep(60)

    #After all jobs have finished, report which jobs failed
    for job in mapping_completed_jobids:
        if mapping_completed_jobids[job] != "COMPLETED":
            print("Kallisto mapping job %s failed with code: %s for file %s"%(job,mapping_completed_jobids[job],mapping_jobid_filename_dict[mapping_jobid]))
            
    #Check that the final Kallisto output is available, if so, remove intermediate files (SRA, fastq)
    for sample in sample_ncbi_dict]:
        for sra in sample_ncbi_dict[sample]:
            if os.path.isfile("%s/%s/%s/abundances.h5"%workign_dir,sample_bioproj_dict[sample],sample):
                if os.path.isfile('%s/%s/fastq/%s_1.fastq.gz'%(working_dir,sample_bioproj_dict[sample],sra)):
                    proc = Popen('rm %s/%s/fastq/%s*'%(working_dir,sample_bioproj_dict[sample],sra),shell=True)
                if os.path.isfile('%s/%s/sra/%s.sra'%(working_dir,sample_bioproj_dict[sample],sra)):
                    proc = Popen('rm %s/%s/sra/%s.sra'%(working_dir,sample_bioproj_dict[sample],sra),shell=True)
            else:
                print("Something happened with SRA: %s for sample: %s"%(sra,sample))   
    '''
    
    now = datetime.datetime.now()
    print('Scripted finished: %s'%now)     

if __name__ == "__main__":
    main()