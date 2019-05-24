#!/usr/bin/env perl -w
#
# a quick&dirty script to convert "yass -d 1 " output files into
# - axt format
# - fasta format
# - blast format
# usage: yass2blast.pl {-o <outputFile>} [-blast|-fasta|-axt] <yassOutputFiles>\n";
#
use strict;


#
# Print an alignment according to $format = {"fasta",blast","blast header","axt"}
#

my $timeMod    = 0;
my %hash       = ();
my $hash_count = 0;

my $fh           = *STDOUT;

sub PrintAlign() {
    # 1) read parameters
    my ($nbAlignments,$nbSSequences,$pos_q_str,$pos_q_end,$pos_s_str,$pos_s_end,$name_q,$name_s,$size_q,$size_s,$reverse,$evalue,$score,$bitscore,$bias,$ts,$tv,$proba,$entropy,$fasta1,$middle,$fasta2,$format,$newchunk_q,$newchunk_s) = @_ ;

    # 2) reverse and complement sequences to keep the "query order"
    if ($reverse eq "-") {
        # change first file positions (increasing order)
        my $tmp = $pos_q_str;
        $pos_q_str = $pos_q_end;
        $pos_q_end = $tmp;
        # sequence complement
        $fasta1 =~ tr /ATUGCRYKMSWVBDHatugcrykmswvbdh/TAACGYRMKWSBVHDtaacgyrmkwsbvhd/;
        $fasta1 = reverse($fasta1);
        $fasta2 =~ tr /ATUGCRYKMSWVBDHatugcrykmswvbdh/TAACGYRMKWSBVHDtaacgyrmkwsbvhd/;
        $fasta2 = reverse($fasta2);
        $middle = reverse($middle);
    }
    
    # 3) format selection
    if ($format eq "fasta") { 
        # 3.1) fasta but with "-" to keep alignments
        print $fh ">".$nbAlignments."a||".$name_q.": [".$pos_q_str."-".$pos_q_end ."]".$reverse."\n";
        print $fh $fasta1."\n";
        print $fh ">".$nbAlignments."b||".$name_s.": [".$pos_s_str."-".$pos_s_end ."]".$reverse."\n";
        print $fh $fasta2."\n";

    } elsif ($format eq "axt") { 
        # 3.2) axt (blastz) format
        print $fh $nbAlignments."\t".$name_q."\t".$pos_q_str."\t".$pos_q_end ."\t".$name_s."\t".$pos_s_str."\t".$pos_s_end ."\t".$reverse."\t".$score."\n";
        print $fh $fasta1."\n";
        print $fh $fasta2."\n";
        
    } elsif ($format eq "blast") { 
        # 3.3) blast like format
        my $nbMismatches = ($middle =~ tr /:\./  /);
        my $nbMatches    = ($middle =~ tr /\|/\|/);
        # a) print Sequence Header
        if ($newchunk_s){
	    if (!(exists $hash{$name_s})) {
		$hash{$name_s} = $hash_count++;
	    }
	    my $id_5 = sprintf("%-5d", ($hash{$name_s}+$timeMod)%65536);
            print $fh ">lcl|".$id_5." ".$name_s."\n";
            print $fh "Length=".$size_s."\n";
            print $fh "\n";
        }

        # b) print Alignment Header
	if ($evalue eq "0") {
	    printf $fh (" Score = %d bits (%d),  Expect = 0.0\n",$bitscore,$score);
	} else {
	    printf $fh (" Score = %d bits (%d),  Expect = %.2g\n",$bitscore,$score,$evalue);
	}
        print $fh " Identities = ".$nbMatches."/".(length($middle))." (".(int(100*$nbMatches/length($middle)))."%), Gaps = ".(length($middle)-$nbMatches-$nbMismatches)."/".(length($middle))." (".(int(100*(length($middle)-$nbMatches-$nbMismatches)/length($middle)))."%)\n";
        print $fh " Strand=Plus/";
        if ($reverse eq "-"){
            print $fh "Minus"."\n";
        }else{
            print $fh "Plus"."\n";
        }
        print $fh "\n";
        

        # c) print Alignment
        my $c = 0;
        my $i_fasta1 = $pos_q_str;
        my $i_fasta2 = 0;
        if ($reverse eq "-") {
            $i_fasta2 = $pos_s_end;
        } else {
            $i_fasta2 = $pos_s_str;
        }
        while ($c < length($middle)){
            my $subfasta1 = substr($fasta1,$c,60);
            my $submiddle = substr($middle,$c,60);
            my $subfasta2 = substr($fasta2,$c,60);
            my $letters_subfasta1 = $subfasta1; $letters_subfasta1 =~ s/-//g;
            my $letters_subfasta2 = $subfasta2; $letters_subfasta2 =~ s/-//g;
            my $nbletters_subfasta1 = length($letters_subfasta1);
            my $nbletters_subfasta2 = length($letters_subfasta2);
            printf $fh ("Query  %-9d ",$i_fasta1); $i_fasta1 += $nbletters_subfasta1-1;
            print  $fh ($subfasta1);
            printf $fh ("  %d\n",      $i_fasta1); $i_fasta1 += 1;

            printf $fh ("                 ");
            print  $fh ($submiddle);
            printf $fh ("\n");

            if ($reverse eq "-") {
                printf $fh ("Sbjct  %-9d ",$i_fasta2); $i_fasta2 -= $nbletters_subfasta2-1;
                print  $fh ($subfasta2);
                printf $fh ("  %d\n",      $i_fasta2); $i_fasta2 -= 1;
            } else {
                printf $fh ("Sbjct  %-9d ",$i_fasta2); $i_fasta2 += $nbletters_subfasta2-1;
                print  $fh ($subfasta2);
                printf $fh ("  %d\n",      $i_fasta2); $i_fasta2 += 1;
            }
            $c += 60;     
            print $fh "\n\n";
        }
    } elsif ($format eq "blastheader") { 
        # 3.4) just the header of blast
        if ($newchunk_s){
	    if (!(exists $hash{$name_s})) {
		$hash{$name_s} = $hash_count++;
	    }
	    my $id_5 = sprintf("%-5d", ($hash{$name_s}+$timeMod)%65536);

            if (length ($name_s) > 50) {
                print $fh "lcl|".$id_5." ".(substr($name_s,0,50))."... "
            } else {
                print $fh "lcl|".$id_5." ".$name_s;
                for(my $e = length ($name_s); $e <= 53 ;$e ++) {
                    print $fh " ";
                }
            }
	    if ($evalue eq "0") {
		printf $fh ("%7.1f   0.0\n",$bitscore);
	    } else {
		printf $fh ("%7.1f   %.2g\n",$bitscore,$evalue);
	    }
        }#1_0 embl|AF417609|AF417609 Colpidium campylum telomerase RNA gen...    32   0.063
    }
}


sub ScanFile($$) {

    my ($yassOutputFile,$format) = @_ ;

    open(FIC,$yassOutputFile) or die print "cant open ".$yassOutputFile."\n";

    my $nbAlignments = 0;
    my $nbSSequences = 0;
    my $selector     = 0;
    my $previous_name_s = "";
    my $previous_name_q = "";

    my @last_align = ();


    #A) alignments
    my ($pos_q_str,
        $pos_q_end,
        $pos_s_str,
        $pos_s_end,
        $name_q,
        $name_s,
        $size_q,
        $size_s,
        $reverse,
        $evalue,
        $score,
        $bitscore,
        $bias,
        $ts,
        $tv,
        $proba,
        $entropy,
        $fasta1,
        $middle,
        $fasta2,)
        =
       (
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "");

    if ($format eq "blast") {
	print $fh "BLASTN 2.2.25+\n";
	print $fh "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro\n";
	print $fh "A. Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and\n";
	print $fh "David J. Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new\n";
	print $fh "generation of protein database search programs\", Nucleic\n";
	print $fh "Acids Res. 25:3389-3402.\n";
	print $fh "\n";
	print $fh "\n";
	print $fh "RID: none\n\n\n";
    }


    while (my $line = <FIC>){

        # potential header (for alignemnts)
        if (($format eq "blast") && !($name_q eq $previous_name_q)) {
            # push the memorized blast alignments
            if ((scalar @last_align) > 0) {
		print $fh "\nALIGNMENTS\n";
                foreach my $m (@last_align) {
                    &PrintAlign($m->{nbAlignments},$m->{nbSSequences},
                                $m->{pos_q_str},$m->{pos_q_end},$m->{pos_s_str},$m->{pos_s_end},
                                $m->{name_q},$m->{name_s},$m->{size_q},$m->{size_s},
                                $m->{reverse},$m->{evalue},$m->{score},$m->{bitscore},
                                $m->{bias},$m->{ts},$m->{tv},$m->{proba},$m->{entropy},
                                $m->{fasta1},$m->{middle},$m->{fasta2},$format,
                                $m->{newchunk_q},
                                $m->{newchunk_s}
                    );
                }
                @last_align = ();
            }
            print $fh "Query= ".$name_q."\n\n";
            print $fh "Length=".$size_q."\n\n\n";
            print $fh "                                                                   Score    E\n";
            print $fh "Sequences producing significant alignments:                       (Bits) Value\n\n";
            $previous_name_q = $name_q;
	    $previous_name_s = "";
        }
        
        # lines begining with "*"
        if($line  =~ /^\*/ ){
            
            # lines begining with "*("
            if($line  =~ /^\*\(/  ) {

                # print previous alignments
                if (!($middle eq "")){
                    if ($format eq "blast") {                   
                        push (@last_align,({
                            nbAlignments => $nbAlignments,
                            nbSSequences => $nbSSequences,
                            pos_q_str => $pos_q_str,
                            pos_q_end => $pos_q_end,
                            pos_s_str => $pos_s_str,
                            pos_s_end => $pos_s_end,
                            name_q => $name_q,
                            name_s => $name_s,
                            size_q => $size_q,
                            size_s => $size_s,
                            reverse => $reverse,
                            evalue => $evalue,
                            score => $score,
                            bitscore => $bitscore,
                            bias => $bias,
                            ts => $ts,
                            tv => $tv,
                            proba => $proba,
                            entropy => $entropy,
                            fasta1 => $fasta1,
                            middle => $middle,
                            fasta2 => $fasta2,
                            newchunk_q => (($previous_name_q eq "") || !($previous_name_q eq $name_q)),
                            newchunk_s => (($previous_name_s eq "") || !($previous_name_s eq $name_s) || ($previous_name_q eq "") || !($previous_name_q eq $name_q))
                       }));
                    }
                    &PrintAlign($nbAlignments,$nbSSequences,
                                $pos_q_str,$pos_q_end,$pos_s_str,$pos_s_end,
                                $name_q,$name_s,$size_q,$size_s,
                                $reverse,$evalue,$score,$bitscore,
                                $bias,$ts,$tv,$proba,$entropy,
                                $fasta1,$middle,$fasta2,(($format eq "blast")?"blastheader":$format),
                                ($previous_name_q eq "") || !($previous_name_q eq $name_q),
				($previous_name_s eq "") || !($previous_name_s eq $name_s) || ($previous_name_q eq "") || !($previous_name_q eq $name_q)
			);
                    $fasta1 = "";
                    $middle = "";
                    $fasta2 = "";
                }
                
                
                # count number of new alignments (and new sequences)
                $nbAlignments++;
                if (($previous_name_s eq "") || 
                    !($previous_name_s eq $name_s)){
                    $nbSSequences++;
                }


                #first line :
                #*(546486-566813)(515659-536310) Ev: 0 s: 20328/20652 f
                $line  =~ m/\(([0-9]+)-([0-9]+)\)\(([0-9]+)-([0-9]+)\) Ev: (.+) s: ([0-9]+)\/([0-9]+) ([rf])$/;
                # get positions
                $pos_q_str  = $1;
                $pos_q_end  = $2;
                $pos_s_str  = $3;
                $pos_s_end  = $4;
                $evalue     = $5;
                $size_q     = $6;
                $size_s     = $7;
                $reverse    = $8;
                if ($reverse =~ /^f/) {
                    $reverse = "+";
                } else {
                    $reverse = "-";
                }
           } elsif ($line  =~ /^\* score/ ){
               # third line :
               $line  =~  m/([0-9]+).*bitscore = (.*)/;
               $score    = $1 + 0;
               $bitscore = $2 + 0;
           } elsif ($line =~ /^\* mutations/){
                # fourth line :
                #* mutations per triplet 129, 103, 124 (4.54e-04) | ts : 190 tv : 166 | entropy : 5.88743              
                $line  =~ m/^\* mutations per triplet ([0-9]+, [0-9]+, [0-9]+) \((.*)\) \| ts : ([0-9]+) tv : ([0-9]+) \| entropy : (.*)/;
                $bias    = $1;
                $proba   = $2;
                $ts      = $3;
                $tv      = $4;
                $entropy = $5;
            }else{
                # second line :
                $line  =~ m/\*[ ]*"(.+)" \(([0-9]+) bp\)[ ]*\/[ ]*"(.+)" \(([0-9]+) bp\)[ ]*$/;
                $previous_name_q = $name_q;
                $name_q = $1;
                $size_q = $2;
                $previous_name_s = $name_s;
                $name_s = $3;
                $size_s = $4;
            }
        } else {
            # alignment :
            if($line =~ /^[A-Za-z-]+$/){
                $line =~ s/\n//;
                $selector++;
                if ($selector % 2){
                    $fasta1 .= "".$line;
                } else {
                    $fasta2 .= "".$line;
                }
            } elsif ($line =~ /^[| :.]+$/ && (length($fasta1) > length($fasta2)) ){  
                $line =~ s/\n//;
                $middle .= $line;
            }
        }
    } # while <FIC>


    # print previous alignments
    if (!($middle eq "")){
       if ($format eq "blast") {                        
          push (@last_align,({
              nbAlignments => $nbAlignments,
              nbSSequences => $nbSSequences,
              pos_q_str => $pos_q_str,
              pos_q_end => $pos_q_end,
              pos_s_str => $pos_s_str,
              pos_s_end => $pos_s_end,
              name_q => $name_q,
              name_s => $name_s,
              size_q => $size_q,
              size_s => $size_s,
              reverse => $reverse,
              evalue => $evalue,
              score => $score,
              bitscore => $bitscore,
              bias => $bias,
              ts => $ts,
              tv => $tv,
              proba => $proba,
              entropy => $entropy,
              fasta1 => $fasta1,
              middle => $middle,
              fasta2 => $fasta2,
              format => $format,
              newchunk_q => (($previous_name_q eq "") || !($previous_name_q eq $name_q)),
              newchunk_s => (($previous_name_s eq "") || !($previous_name_s eq $name_s) || ($previous_name_q eq "") || !($previous_name_q eq $name_q))
          }));
       }
       &PrintAlign($nbAlignments,$nbSSequences,
                      $pos_q_str,$pos_q_end,$pos_s_str,$pos_s_end,
                      $name_q,$name_s,$size_q,$size_s,
                      $reverse,$evalue,$score,$bitscore,
                      $bias,$ts,$tv,$proba,$entropy,
                      $fasta1,$middle,$fasta2,(($format eq "blast")?"blastheader":$format),
                      ($previous_name_q eq "") || !($previous_name_q eq $name_q),
                      ($previous_name_s eq "") || !($previous_name_s eq $name_s) || ($previous_name_q eq "") || !($previous_name_q eq $name_q)
       );
       $fasta1 = "";
       $middle = "";
       $fasta2 = "";
    }

    if ($format eq "blast") {
       # push the memorized blast alignments
       if ((scalar @last_align) > 0) {
	   print $fh "\nALIGNMENTS\n";
           foreach my $m (@last_align) {
                &PrintAlign($m->{nbAlignments},$m->{nbSSequences},
                            $m->{pos_q_str},$m->{pos_q_end},$m->{pos_s_str},$m->{pos_s_end},
                            $m->{name_q},$m->{name_s},$m->{size_q},$m->{size_s},
                            $m->{reverse},$m->{evalue},$m->{score},$m->{bitscore},
                            $m->{bias},$m->{ts},$m->{tv},$m->{proba},$m->{entropy},
                            $m->{fasta1},$m->{middle},$m->{fasta2},$format,
                            $m->{newchunk_q},
                            $m->{newchunk_s}
                );
            }
            @last_align = ();
       }
    }
    close FIC;
    return $nbAlignments;
}





#------#
# Main #
#------#

my $outputFormat = "blast";
my $outputFile   = undef;

my @timeData = localtime(time);
my $timeDataJoin = join('', @timeData);
   $timeMod  = ($timeDataJoin % 65536);


($#ARGV >= 0) or die "yass2blast:\n".
                     "  a quick&dirty script to convert \"yass -d 1 \" output files into\n".
                     "  - axt format\n".
                     "  - fasta format\n".
                     "  - blast format\n".
                     "\n".
                     "  usage: yass2blast.pl {-o <outputFile>} [-blast|-fasta|-axt] <yassOutputFiles>\n".
                     "\n";


for (my $i = 0 ; $i <= $#ARGV ; $i++) {
    # a) select on the fly parameters ...
    if (($ARGV[$i]) eq "-blast") {
        $outputFormat = "blast";
    } elsif (($ARGV[$i]) eq "-fasta") {
        $outputFormat = "fasta";
    } elsif (($ARGV[$i]) eq "-axt") {
        $outputFormat = "axt";
    } elsif (($ARGV[$i]) eq "-o") {
        $i = $i+1;
        $i <= $#ARGV or die "-o found without parameter";
        $outputFile = $ARGV[$i];
	open($fh , '>>', $outputFile) or die "Unwritable \"$outputFile\" file";
    } else {
        # b) or scan files
        &ScanFile($ARGV[$i],$outputFormat);
	if (defined $outputFile) {
	    close $fh;
	    $fh = *STDOUT;
	    $outputFile = undef;
	}
    }
}
