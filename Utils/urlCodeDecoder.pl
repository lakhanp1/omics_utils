#!/usr/bin/perl -w

use strict;
use warnings;

while(<>){
	$_ =~ s/\%([A-Fa-f0-9]{2})/pack('C', hex($1))/seg;
    #$_ =~ s/\+/ /g;
	print $_;
}