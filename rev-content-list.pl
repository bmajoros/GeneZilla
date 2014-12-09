#!/usr/bin/perl
while(<STDIN>)
{
    @a=split/\s+/;
    next unless @a==3;
    next unless($_=~/^\d+\s+\S+/);
    $pos=4924-$a[0];
    $cons=$a[1];
    if($cons=~/NEG-/){$cons=~s/NEG-//}
    else{$cons="NEG-$cons"}
    push @out,"$pos $cons $a[2]\n";
}
while(@out)
{
    $x=pop @out;
    print $x;
}
