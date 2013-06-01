#!/bin/bash

modules=`pwd`'/source/modules'
main=`pwd`'/source/main.f90'


# first create the mod objects

for compLevel in 'dmod1' 'dmod2' 'dmod3' 'dmod4' 'dmod5' ; do
   echo 'Level ' $compLevel
   for file in $modules/*.$compLevel ; do
   
      f90filename=`echo "$file" | cut -d'.' -f1`.f90
      echo ' file ' $file

      cp $file $f90filename
      ifort -g -c $f90filename
      rm $f90filename
   done
done

# now compile main.f90
echo 'Compiling main...'
ifort $main *.o -I`pwd` -o Grid.x
   

# clean a bit
rm *.o
rm *.mod

