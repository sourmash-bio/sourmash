#! /bin/bash
for i in $(cat urls.txt)
do
filename=$(basename $i)
if [ \! -f $filename ]; then
  echo DNE, downloading $filename
  curl -L $i | gunzip -c | head -4000000 | gzip -9c > $filename
else
  echo $filename exists, not re-getting
fi
   
done
