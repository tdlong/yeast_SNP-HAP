use warnings;
use strict;

# perl scripts/make.muttable.general.pl "May1/mut/" ".newmut.vcf" <$input >temp.vcf

my $locationOfMUTs = $ARGV[0];
my $fileExtension = $ARGV[1];

print("chr\tpos\tsample\tref\talt\tNref\tNalt\n");
while (my $name = <STDIN>){
	chomp $name;
	my $filename = $locationOfMUTs.$name.$fileExtension;
	open (FH, "$filename") or die "$filename: no go";
	while (my $line = <FH>){
		chomp $line;
		if ($line !~ m/^#/){
			my @ff = split("\t",$line);
			if((length($ff[3]) == 1) and (length($ff[4]) == 1)){
				my @info = split('\:',$ff[9]);
				if ($info[4] !~ m/\,/ ){
					print("$ff[0]\t$ff[1]\t$name\t$ff[3]\t$ff[4]\t");
					print("$info[2]\t$info[4]\n");
					}
				}
			}
		}
	}

