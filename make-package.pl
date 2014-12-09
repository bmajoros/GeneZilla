#!/usr/bin/perl
use strict;

my $date=`date +"%m-%d-%y"`;
chomp $date;
my $outfile="genezilla.tar";
my $T="genezilla";

my $tarCmd="cd ..; tar hcf $T/$outfile $T/*.model $T/makefile $T/*.[HC] BOOM/*.[HCch] perl/*.pl perl/*.pm $T/train.sh $T/*.pl $T/*.iso $T/*.cfg $T/*.top $T/web/*.html $T/web/*.png $T/web/*.gif $T/xgraph.tar.gz";
print "Tarring into $outfile...\n$tarCmd\n";
system("$tarCmd");
print "Zipping...\n";
system("gzip $outfile");
print "done.  archive is in $outfile.gz\n";

print "\n###########################################################\n";
print "Now follow these steps:\n";
print "  1) ftp this file to thrawn, put in directory:\n";
print "     /local2/Netscape/prod/ftp_deck/ftp/pub/software/GeneZilla\n";
print "  2) telnet to thrawn and execute these commands:\n";
print "     cd  /local2/Netscape/prod/ftp_deck/ftp/pub/software/GeneZilla\n";
print "     conman -Ff pub/software/GeneZilla/genezilla.tar.gz\n";
print "###########################################################\n\n";











