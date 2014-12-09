#!/usr/bin/perl
use strict;

my $usage="$0 <mdd.model>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my @nodes;
open(IN,$filename) || die "can't open $filename\n";
#<IN>;
#my $line=<IN>;
#my @fields=split/\s+/$line;

while(<IN>)
  {
    if(/internal/)
      {
	my $line=<IN>;
	my @fields=split/\s+/,$line;
	my ($pos,$l,$symbols)=@fields;
	push @nodes,{type=>"internal",pos=>$pos,sym=>$symbols};
      }
    elsif(/leaf/)
      {
	push @nodes,{type=>"leaf"};
      }
  }
close(IN);

my $n=@nodes;
my $depth=0;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $node=$nodes[$i];
    if($node->{type} eq "leaf"){print "L"}else{print "+"}
  }
print "\n";

dumpTree(0);
sub dumpTree
  {
    my ($depth)=@_;
    my $node=shift @nodes;
    return if($node->{type} eq "leaf");

    my $margin="\t"x$depth;
    my $pos=$node->{pos};
    my $sym=$node->{sym};
    print "$margin$pos=$sym\n";

    dumpTree($depth+1);
    dumpTree($depth+1);
  }

