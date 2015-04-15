#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <mu> <delta> <epsilon>\n" unless @ARGV==3;
my ($mu,$delta,$eps)=@ARGV;

print "\tU\tN\tI\tE0\tE1\tE2\n";
print "U\t1\t1\t1\t1\t1\t1\n";
print "N\t1\t1\t$eps\t$eps\t$eps\t$eps\n";
print "I\t1\t$eps\t1\t$delta\t$delta\t$delta\n";
print "E0\t1\t$eps\t$mu\t1\t$eps\t$eps\n";
print "E1\t1\t$eps\t$mu\t$eps\t1\t$eps\n";
print "E2\t1\t$eps\t$mu\t$eps\t$eps\t1\n";


