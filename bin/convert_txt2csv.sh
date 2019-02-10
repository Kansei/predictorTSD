#!/bin/zsh

awk -F '[\t]' -v 'OFS=,' '{split($21,barcode,"_Grn");print $1,$3,$6,$8,$11,barcode[1]}' $1 | sort -u -r > $2
