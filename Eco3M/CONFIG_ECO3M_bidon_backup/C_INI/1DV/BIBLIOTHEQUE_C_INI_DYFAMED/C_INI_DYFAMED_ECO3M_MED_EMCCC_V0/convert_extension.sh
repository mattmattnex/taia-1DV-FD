#!/bin/bash
old_ext='DAT_OK'
for file in *.$old_ext
do
	mv -v "$file" "${file%.$old_ext}.DAT"
done;

