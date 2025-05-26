use strict;
use Bio::SeqIO;

die "perl $0 <in.fa> <window> <out>" if ($#ARGV<2);

my $window=$ARGV[1];
open IN,"<$ARGV[0]" or die "1";
open OUT,">$ARGV[2]" or die "2";

my $seqin=Bio::SeqIO->new(-format=>"fasta",-file=>"$ARGV[0]");
while (my $obj=$seqin->next_seq){
	my $id=$obj->id;
        $id=~s/\|/\_/g;
        my $len=$obj->length;
        my $seq=$obj->seq;
	my$i;
	for($i=0;$i<$len;$i+=$window){
		my $subseq=substr($seq,$i,$window);
		my $start=$i+1;
		my $end=$i+$window;
		$end=$len if($end>$len);
		my $gc=$subseq=~tr/GCgc/GCgc/;
		my $gc_con=$gc/($end-$start);
		print OUT "$id\t$start\t$end\t$gc_con\n";
		}
	}

close OUT;
