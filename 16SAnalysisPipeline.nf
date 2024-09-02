params.read=/disk6/11.SH_MouseGut/1.Raw_data


workflow{

}
process Pear{


}

import sys,re,os,string,glob
import threading,time

#################################
# main
#################################

#Step 1: Trimmomatic to remove adapter sequences (in order to generate clean read)
Pear = "/disk6/11.SH_MouseGut/Program/pear-0.9.11-linux-x86_64/bin/pear"
my_path = "/disk6/11.SH_MouseGut/"
data_path = "1.Raw_data"

#os.mkdir(my_path + '2.Pear')
count = 1;

outfile = open("Pear_cmd.log","w")

os.system("mkdir 2.Pear")

for id_left in glob.glob("./"+data_path+"/*_1.fastq.gz"):
	print(id_left)
	id_right = id_left.replace("_1.fastq.gz", "_2.fastq.gz")
	ID  = os.path.split(id_left)[1].replace("_1.fastq.gz", "")

	cmd = "%s -f %s -r %s -j 12 -o ./2.Pear/%s > ./2.Pear/%s.log"%(Pear,id_left,id_right,ID,ID)
	print(cmd)
	outfile.write(cmd+"\n")
	os.system(cmd)


	count = count +1
