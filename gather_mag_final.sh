#!/bin/bash

rm magmom*

#!/bin/bash

# Check if the correct number of arguments were provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 Natoms outfile"
    exit 1
fi

Natoms="$1"
outfile="$2"


echo "NATOMS: $NATOMS"
echo "FILE: $FILE"



b=$(( Natoms+3 ))


	
bzgrep -A $b 'magnetization (x)' $outfile | tail -$Natoms | awk '{print $NF}' >mag_x
bzgrep -A $b 'magnetization (y)' $outfile | tail -$Natoms | awk '{print $NF}' >mag_y
bzgrep -A $b 'magnetization (z)' $outfile | tail -$Natoms | awk '{print $NF}' >mag_z
paste mag_x mag_y mag_z >magmom

awk '{print sqrt($1*$1+$2*$2+$3*$3)}' magmom >> magmom_size


rm mag_*


