#!/usr/bin/perl
while(<STDIN>)
{
    @a=split/\s+/;
    next unless @a==3;
    $p=4924-$a[0];
    $s=$a[1];
    $s=~s/-//g;
    $l=length($s);
    if($l>3) {$l=6}
    $p-=$l;
    $cons=$a[1];
    if($cons=~/-/){$cons=~s/-//}
    else{$cons="-$cons"}
    push @out,"$p $cons $a[2]\n";
}
while(@out)
{
    $x=pop @out;
    print $x;
}
