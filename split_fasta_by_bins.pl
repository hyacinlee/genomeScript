#!/usr/bin/perl  -w 
use strict;
die "Usage: < In fasta > < sample name > < Bin size >\n" unless(@ARGV==3);
my (%fasta,$id,@id);
`ln -sf $ARGV[0] $ARGV[1].Contig.fasta`;
open FA,"<$ARGV[0]" or die;
while(<FA>){
	chomp;
	if(/>/){
		$id=$_;
		$id=~s/>//g;
		push @id,$id;
		next;
	}
	$fasta{$id}.=$_;
}
close FA;

my $num=0;
open ID,">$ARGV[1].bin.id";
open LEN,">$ARGV[1].Contig.len";
open BIN,">$ARGV[1].Bins.fasta";
foreach my $fa(@id){
	my $len=length $fasta{$fa};
	print LEN "$fa\t$len\n";
	my $id=1;
	my $start=0;
	while($start<$len){
		my $step=$ARGV[2];
		$step=$len-$start if($start+$ARGV[2]>$len);	
		my $seq=substr($fasta{$fa},$start,$step);
#		next if ($seq=~m/^(N*)$/);
#		next if ($seq=~m/^(n*)$/);
##next if ($seq=~m/^n(.*)$/)
                my $end = $start + $step;
		my $ss  = $start+1;
                print ID "$fa\.$id\t1\t$step\t$fa\t$ss\t$end\t+\t$num\n";
		print BIN ">$fa\_$id\n$seq\n";
		$start+=$step;
		$id++;
		$num++;
	}
}
close ID;
close LEN;
close BIN;
