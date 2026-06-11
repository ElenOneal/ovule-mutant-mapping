#!/bin/bash
#rename sequencing files
ls -1 JW*gz > sequencing_files.txt
#take second column from library layout
cut -f2 ovule_plate_layout.txt > library_names.txt
#produce new file copying each line from library_names.txt 4x to make new file
perl -ne 'print "$_" x4' library_names.txt > library_names_4x
#strip first two "_" to get file type from first column of sequencing_files.txt
cut -d'_' -f4 sequencing_files.txt > file_type
sed -i 's/I1/.i7/g' file_type 
sed -i 's/R1/.1/g' file_type 
sed -i 's/R2/.i5/g' file_type 
sed -i 's/R3/.2/g' file_type 
paste library_names_4x file_type > test
#replace tab with ""
sed -i 's/\t//g' test
#add .fq.gz to each in test
sed -e 's/$/.fq.gz/' -i test
#paste columns sequencing_files.txt library_names_4x test > test2
paste sequencing_files.txt test > mvcommand
#prefix mvfile with mv command
sed -i 's/JW/mv JW/g' mvcommand
#replace tab with space (because it's ugly not because it matters)
sed -i 's/\t/ /g' mvcommand
echo '#!/bin/bash' >> mv.sh
echo '#' >> mv.sh
echo '#$ -S /bin/bash' >> mv.sh 
cat mv.sh mvcommand > mvscript.sh
chmod 755 mvscript.sh
mv library*gz 01_rmdup 
./mvscript.sh