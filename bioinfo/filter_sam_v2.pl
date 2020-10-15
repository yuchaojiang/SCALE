# Filter SAM file for pair-end sequencing. Require

# 1) Mapping quality score of each read is >= 30

# 2) Read pair mapped to the same chromosome with expected orientation; distance < 500,000 bp

# 3) Allow each read to map to <=10 positions in the genome

# 4) Keep reads mapped to MT

# 5) Keep soft clipped or hard clipped reads

# 6) CIGAR should have an "M" in it

 

 

use strict;

my $sam = $ARGV[0];

my $out = $ARGV[1];

 

my @temp;

 

open OUT, ">$out"||die "$!";

 

my $count = 0;

 

open IN, "samtools view -hS $sam | "||die "$!";

while (<IN>)

{

        chomp;

        my $line = $_;

 

        my $output;             # Whether to output this line

 

        if (m/^\@/)             # Keep the header line

        {

                print OUT "$line\n";

                next;

        }

 

#       $count ++;

#       if (int($count / 10000) == $count / 10000)

#       {

#               print "Processing $count read\r";

#       }

        my ($qname, $flag, $chr, $pos, $mapQ, $CIGAR, $rnext, $pnext, $length, $read, $qual, @tags) = split(/\t/, $line);

 

#       if ($chr =~m/M/i)       # Exclude reads mapped to MT

#       {

#               next;

#       }

 

        my ($strand_r1, $strand_r2);

        if ($flag & 16)

        {       $strand_r1 = "-";}

        else

        {       $strand_r1 = "+";}

 

        if ($flag & 32)

        {       $strand_r2 = "-";}

        else

        {       $strand_r2 = "+";}

 

        # Require proper pair; R1 and R2 on different strand; same chromorome; distance <= 500000 bp; mapQ >=30; CIGAR string contains "M"

        if ($flag & 2 && $strand_r1 ne $strand_r2 && $rnext eq "=" && abs($pos - $pnext) <= 500000 && $mapQ >= 30 && $CIGAR =~m/M/)

        {

                # Require opposite orientation and proper positions like this ---->  <----

                if ($strand_r1 eq "+" && $pos <= $pnext || $strand_r1 eq "-" && $pos >= $pnext)

                {

                        # Require <= 10 mapping

                        foreach my $tag (@tags)

                        {

                                if ($tag =~m/NH:i:(\d)/)

                                {

                                        my $n_hit = $1;

                                        if ($n_hit <= 10)

                                        {

                                                $output = 1;

                                        }

                                        else

                                        {

                                                $output = 0;

                                        }

                                }

                        }

#                       # Exclude soft clipped or hard clipped reads

#                       if ($CIGAR =~m/[SH]/)

#                       {

#                               $output = 0;

#                       }

                }

        }

 

        if ($output == 1)

        {

                print OUT "$line\n";

        }

}

close IN;

close OUT;

 

print STDERR "Filtering finished.\n";
