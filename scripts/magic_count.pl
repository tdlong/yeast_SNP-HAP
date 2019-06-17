use strict;
use warnings;
use List::Util qw/first/;

while (my $line = <STDIN>){
	chomp $line;
	my @Line = split("\t",$line);
	my $POS = $Line[0]."_".$Line[1]."_".$Line[3];
	my @Info = split(";",$Line[7]);
#	print join("xx",@Info),"\n";
	my $MM = first { $_ =~ /DP4/ } @Info;
	my $DP4 = (split("=",$MM))[1];
	my @vv = split(",",$DP4);
	my $tot=0;
	foreach (@vv) {
		$tot += $_;
		}
	my $nRef = $vv[0] + $vv[1];
	print "$POS\t$DP4\t$nRef\t$tot\n";
	}

