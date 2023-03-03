#!/usr/bin/perl -w
use strict;
use POSIX qw(ceil floor);
use Storable qw(dclone);
use Getopt::Long;
use Data::Dumper;

# Usage:
# filterfasta.pl --match '$chromosome' --start '$start' --end '$end' '$db'/Zm*

#uses:
#1) input fasta only; reformats to one line fasta
#2) pulls out one seq and makes it one line fasta: use --match
#3) pulls out one seq and give substr: use --match and --start --end
#4) only tell me sequence lengths: use --size Y (and --match string)
#5) sort by headers: use --sort (name=sort by name, size=sort by size in rev)
#6) rev-comp sequence: add -rev
#7) pull out multiple sequences: Use match with a file name
#8) directly write individual fastas to file named $ARGV[0].#instead of to screen: use --file Y, use --file head to write to a file with the header name plus '.fasta', header name must be clean

my $match = '';
my $start = '';
my $end = '';
my $size = 'N';
my $order = 'N';
my $rev = 'N';
my $match_type = 'C'; #(F)ull means the full match string matches the full header line; (P)artial means the match string matches part of the header line; (C)omplete means the match string matches between the '>' and a space/tab/eol
my $width = 'F';
my $verbose = 'N';
my $file = 'N';

GetOptions('match=s' => \$match, 'width=s' => \$width, 'start=i' => \$start, 'end=i' => \$end, 'size=s' => \$size, 'order=s' => \$order, 'rev=s' => \$rev, 'file=s' => \$file, 'match_type=s' => \$match_type, 'verbose=s' => \$verbose);

if($start ne '' && $end ne '') {
  $start -= 1;
  $end -= 1;
}
if($start) {
  if($start < 0) {
    $start = 1;
  }
}
my $fasta_ref = LOAD_FASTA($ARGV[0]);
my %fasta = %$fasta_ref;
if($verbose eq 'Y' || $verbose eq 'YY') {
  print "fasta loaded\n";
}

my %out = ();
if($match ne '') {
  my @match = ();
  if(-f $match) {
    open(MATCH, $match);
    while(<MATCH>) {
      chomp($_);
      $_ =~ s/\|/\_\_\_\_\_/g;
      push @match, $_;
    }
    close(MATCH);
  } else {
    $match =~ s/\|/\_\_\_\_\_/g;
    push @match, $match;
  }
  if($verbose eq 'Y' || $verbose eq 'YY') {
    print "match loaded\n";
  }
  foreach my $m (@match) {
    if($verbose eq 'YY') {
      print "searching $m\n";
    }
    foreach my $head (keys %fasta) {
      my $head_fix = $head;
      $head_fix =~ s/\_\_\_\_\_/\|/g;
      if($match_type eq 'F') {
        if($head eq $m) {
          if($start ne '' && $end ne '') {
            $out{$head_fix} = substr($fasta{$head},$start,($end-$start)+1);
          } else {
            $out{$head_fix} = $fasta{$head};
          }
        }
      } elsif($match_type eq 'P') {
        if($head =~ m/$m/) {
          if($start ne '' && $end ne '') {
            $out{$head_fix} = substr($fasta{$head},$start,($end-$start)+1);
          } else {
            $out{$head_fix} = $fasta{$head};
          }
        }
      } elsif($match_type eq 'C') {
        if($head =~ m/^$m / || $head =~ m/^$m\t/ || $head =~ m/^$m$/) { #space one works, not positive about rest.  The ^$m$ is the same as Full
          if($start ne '' && $end ne '') {
            $out{$head_fix} = substr($fasta{$head},$start,($end-$start)+1);
          } else {
            $out{$head_fix} = $fasta{$head};
          }
        }
      }
    }
  }
} else {
  %out = %fasta;
}

#write results
if($order eq 'N') {
  my $filecount = 0;
  foreach my $head (keys %out) {
    if($size eq 'Y') {
      my $len = length($out{$head});
      print "$head $len\n";
    } elsif($size eq 'N') {
      if($rev eq 'Y') {
        my $tmp = reverse($out{$head});
        $tmp =~ tr/ACGTacgt/TGCAtgca/;
        $out{$head} = $tmp;
      }
      if($width eq 'T') {
        $out{$head} =~ s/(.{1,100})/$1\n/gs;
      }
      chomp($out{$head});
      if($file eq 'N') {
        print ">$head\n";
        print "$out{$head}\n";
      } elsif($file eq 'Y') {
        open(OUT,">$ARGV[0].$filecount");
        print OUT ">$head\n";
        print OUT "$out{$head}\n";
        close(OUT);
        $filecount++;
      } elsif($file eq 'head') {
        open(OUT,">$head.fasta");
        print OUT ">$head\n";
        print OUT "$out{$head}\n";
        close(OUT);
        $filecount++;
      }
    }
  }
} elsif($order eq 'name') {
  foreach my $head (sort keys %out) {
    if($size eq 'Y') {
      my $len = length($out{$head});
      print "$head $len\n";
    } elsif($size eq 'N') {
      print ">$head\n";
      if($rev eq 'Y') {
        my $tmp = reverse($out{$head});
        $tmp =~ tr/ACGTacgt/TGCAtgca/;
        $out{$head} = $tmp;
      }
      if($width eq 'T') {
        $out{$head} =~ s/(.{1,100})/$1\n/gs;
      }
      chomp($out{$head});
      print "$out{$head}\n";
    }
  }
} elsif($order eq 'size') {
  foreach my $head (sort {length($out{$b}) <=> length($out{$a})} keys %out) {
    if($size eq 'Y') {
      my $len = length($out{$head});
      print "$head $len\n";
    } elsif($size eq 'N') {
      print ">$head\n";
      if($rev eq 'Y') {
        my $tmp = reverse($out{$head});
        $tmp =~ tr/ACGTacgt/TGCAtgca/;
        $out{$head} = $tmp;
      }
      if($width eq 'T') {
        $out{$head} =~ s/(.{1,100})/$1\n/gs;
      }
      chomp($out{$head});
      print "$out{$head}\n";
    }
  }
}


############################
sub LOAD_FASTA {
  my $dna = $_[0];
  my %fasta = ();
  open(DNA,$dna) or die "$dna not found\n";
  my $line = <DNA>;
  $line =~ s/\r//g;
  do {
    if($line =~ m/\>/) {
      chomp($line);
      my $head = $line;
      my $seq = '';
      while(<DNA>) {
        $line = $_;
        $line =~ s/\r//g;
        if($line =~ m/\>/) {
          last;
        }
        chomp($line);
        $seq .= $line;
      }
      $head =~ s/\>//g;
$head =~ s/\|/\_\_\_\_\_/g;
#one value
      $fasta{$head} = $seq;
#multiple values
#      $fasta{$head} = {seq => $seq, more => $something};
    }
  }until(eof);
  close(DNA);
  return(\%fasta);
}