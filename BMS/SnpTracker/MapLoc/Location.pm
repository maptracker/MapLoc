# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Location;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

=head1 SYNOPSIS

=head1 AUTHOR

Charles Tilford <podmail@biocode.fastmail.fm>

//Subject __must__ include 'Perl' to escape mail filters//

=head1 LICENSE

Copyright 2014 Charles Tilford

 http://mit-license.org/

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

=cut

use strict;
use BMS::SnpTracker::MapLoc::Tagged;
use BMS::SnpTracker::MapLoc::Accessioned;
use BMS::SnpTracker::MapLoc::TransientLocation;
use BMS::SnpTracker::MapLoc::LockedLocation;
use BMS::SnpTracker::MapLoc::Range;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Tagged
               BMS::SnpTracker::MapLoc::Accessioned);

our $tinyFakeFreq = 0.000000001234;

sub new {
    my $proto = shift;
    $proto->bench_start();
    my $class = ref($proto) || $proto;
    my $ml    = shift;
    
    my ($chr, $lft, $rgt, $bld) = @_;
    my $pkey;
    if ($#_ == 0) {
        $pkey = $chr;
        my $read = $ml->{STH}{READ_LOCATION} ||= $ml->dbh->prepare
            ( -name => "Read basic location data for PKEY",
              -sql  =>
              "SELECT chr, pos_to_left, pos_to_right, build ".
              "FROM location WHERE loc_id = ?");
        $read->execute($pkey);
        my $rows = $read->fetchall_arrayref();
        if ($#{$rows} == -1) {
            $proto->bench_end();
            return undef;
        }
        ($chr, $lft, $rgt, $bld) = @{$rows->[0]};
    }
    my $self = {
        LFT => $lft,
        RGT => $rgt,
        CHR => $chr,
        BLD => $bld,
        PKEY   => $pkey,
        MAPLOC => $ml,
        CATS   => {},
    };
    $proto->bench_end();
    bless ($self, $class);
}

sub make_transient {
    my $self = shift;
    bless ($self, 'BMS::SnpTracker::MapLoc::TransientLocation');
}

sub lock {
    my $self = shift;
    bless ($self, 'BMS::SnpTracker::MapLoc::LockedLocation');
}

*is_in_db = \&is_known;
sub is_known {
    my $self = shift;
    return ($self->{PKEY}) if ($self->{PKEY});
    $self->bench_start();
    my $rv   = 0;
    my $ml   = $self->maploc();
    my $sth  = $ml->{STH}{CHECK_FOR_LOCATION} ||= $ml->dbh->prepare
        ( -name => "Check for existing location in DB",
          -sql  => "SELECT loc_id FROM location WHERE chr = ? AND pos_to_left = ? AND pos_to_right = ? AND build = ?");
    my @ids = $sth->get_array_for_field
        ( $self->chr, $self->lft, $self->rgt, $self->build );
    $rv = $self->{PKEY} = $ids[0] unless ($#ids == -1);
    $self->bench_end();
    return $rv;
}

*chr_build = \&build;
sub build     { return shift->{BLD}; }
sub chr_builds {
    my $self = shift;
    return wantarray ? ($self->build) : $self->build;
}
sub chr       { return shift->{CHR}; }
sub chr_names {
    my $self = shift;
    return wantarray ? ($self->chr) : $self->chr;
}
*lft = \&left;
sub left     { return shift->{LFT}; }
*rgt = \&right;
sub right    { return shift->{RGT}; }
sub obj_type { return "Location"; }
*chr_start = \&start;
sub start    { return shift->{LFT} + 1; }
*chr_end   = \&end;
sub end      { return shift->{RGT} - 1; }

*left_anchor  = \&left;
*right_anchor = \&right;
*anchor_id    = \&loc_acc_id;

# ::Location
sub hsp {
    my $self = shift;
    unless ($self->{AS_HSP}) {
        $self->bench_start();
        my $hsp  = $self->{AS_HSP} =
            BMS::SnpTracker::MapLoc::Range->new( $self->maploc );
        $hsp->{STRAND} = [1];
        $hsp->hsps( [[ $self->lft, $self->rgt ]] );
        $hsp->seq_ids( $self->loc_acc_id() );
        $hsp->set_anchor( $self->anchor_id() );
        $self->bench_end();
    }
    return $self->{AS_HSP};
}

*len = \&width;
*length = \&width;
*chr_width = \&width;
sub width {
    my $self = shift;
    return $self->right - $self->left - 1;
}

# ::Location
sub to_text {
    my $self = shift;
    $self->bench_start();
    my $pre  = shift || "";
    my $txt  = $pre . $self->full_loc();
    if (my $pkey = $self->handle) {
        $txt .= " [loc_id = $pkey]";       
    } else {
        $txt .= " /ERR - No loc_id/";
    }
    $txt .= "\n";
    my @accs = $self->each_accession();
    my @cats = $self->each_category();
    $txt .= "$pre Categories: ".join(', ', @cats)."\n" unless ($#cats == -1);
    if ($#accs == -1) {
        $txt .= "$pre  -no external accessions-\n";
    } else {
        my $accH = $self->accessions_with_authorities();
        foreach my $acc (@accs) {
            $txt .= sprintf("%s  %15s : %s\n", $pre, $acc, join
                            ('+', sort {uc($a) cmp uc($b) } @{$accH->{$acc}}));
        }
    }
    my $ml      = $self->maploc;
    my $imps    = $self->impact( @_ );
    my @iData   = values %{$imps};
    if ($#iData == -1) {
        $txt .= "$pre  -no RNAs nearby-\n";
    } else {
        my @impacts;
        while (my ($rId, $idat) = each %{$imps}) {
            my $rna = $ml->get_rna_by_id($rId);
            my $sym = $rna->symbol() || '';
            push @impacts, [ $idat->{imp} || "UNK", $sym, $rna->acc(),
                             $idat->{nucNom} || "", $idat->{protNom}  || "" ];
        }
        @impacts = sort { $a->[0] cmp $b->[0] || $a->[1] cmp $b->[1]
                              || $a->[2] cmp $b->[2] } @impacts;
        $txt .= sprintf("%s%d Impact%s:\n", $pre, $#impacts + 1, $#impacts == 0
                        ? '' : 's');
        map { $txt .= sprintf("%s  %3s %5s %-18s %12s %12s\n", $pre, @{$_}) } @impacts;
    }
    my $alleles = $self->each_allele();
    my @aArr    = sort keys %{$alleles};
    if ($#aArr == -1) {
        $txt .= "$pre  -no alleles noted-\n";
    } else {
        my @head = ( "", map { length($_) <= 10 ? $_ : $ml->comma_number(length($_)).'bp' } @aArr);
        my %pops;
        for my $b (0..$#aArr) {
            my $bHash = $alleles->{$aArr[$b]};
            foreach my $pid (keys %{$bHash}) {
                my $targ = $pops{$pid} ||= [ map { "" } @head ];
                my ($freq, $count) = @{$bHash->{$pid}};
                $targ->[$b+1] = !defined $freq ? '' : $count ? 
                    sprintf("%d/%d", 0.5 + $count * $freq, $count) :
                    sprintf("%d%%", int(0.5 + 100 * $freq));
            }
        }
        while (my ($pid, $arr) = each %pops) {
            my $pop = $ml->get_population_fast($pid);
            $arr->[0] = $pop->name();
        }
        my @rows = (\@head);
        push @rows, sort { $a->[0] cmp $b->[0] } values %pops;
        my $frm = "    ";
        my $bar = $frm;
        for my $i (0..$#head) {
            my ($m) = sort { $b <=> $a } map { length($_->[$i]) } @rows;
            my $sfrm;
            if ($i) {
                $sfrm = sprintf(" %%-%ds", $m);
                $sfrm .= " |" unless ($i % 4);
            } else {
                # Population column
                $sfrm = sprintf("%%%ds |", $m);
            }
            $bar .= sprintf($sfrm, '-' x $m);
            $frm .= $sfrm;
        }
        $frm .= "\n";
        $bar .= "\n";
        $txt .= sprintf($frm, @{shift @rows});
        for my $i (0..$#rows) {
            $txt .= $bar unless ($i % 3);
            $txt .= sprintf($frm, @{$rows[$i]});
        }
    }
    
    $txt .= $self->tag_text($pre."  ");
    $self->bench_end();
    return $txt;
}

sub to_html {
    my $self = shift;
    return "<pre>".$self->to_text(@_)."</pre>";
}

# ::Location
*to_canvasXpress = \&cx_genome_part;
sub cx_genome_part {
    my $self = shift;
    $self->bench_start();
    my $args = $self->parseparams( @_ );
    my ($l, $r) = ($self->left, $self->right);
    my ($cl, $cr);
    if (defined $args->{LEFTPOS}) {
        ($cl, $cr) = ($l, $r);
        $l = $args->{LEFTPOS};
        $r = $args->{RIGHTPOS};
    }
    my $wid = $r - $l - 1;
    my ($type, $dat);
    if ($wid == 0) {
        ($type, $dat) = ("Deletion", [[$l, $r]]);
    } else {
        my ($s, $e) = ($l + 1, $r -1);
        $dat = [[$s, $e]];
        if ($s == $e) {
            # Single position
            $type = "Point";
        } else {
            # Range
            $type = "Range";
        }
    }
    my $imp  = $args->{IMPDATA} || $self->impact( @_ );
    my $pop  = $args->{POPDATA};
    unless ($pop) {
        my @pbits;
        foreach my $key ('TOSSPOP', 'KEEPPOP') {
            if (my $val = $args->{$key}) {
                push @pbits, ($key, $val);
            }
        }
        $pop = $self->each_population_id( @pbits );
    }
    my $MAF = 0;
    my %variance;
    while (my ($pid, $alH) = each %{$pop || {}}) {
        my @fs;
        while (my ($base, $fcnt) = each %{$alH}) {
            if (my $f = $fcnt->[0]) {
                push @fs, $f;
                $variance{$base} = $f if (!defined $variance{$base} ||
                                          $variance{$base} < $f);
            }
        }
        @fs = sort { $b <=> $a } @fs;
        if ($fs[1]) {
            # Consider MAF as "all non major alleles" within one population
            my $mf = 1 - $fs[0];
            # Capture the highest MAF across all populations
            $MAF = $mf if ($MAF < $mf);
        }
    }
    # Variance ignores the highest observed allele, and sums the highest
    # observed frequency for all other alleles.
    # I am sure that this is a butchery of a more formal meaning for
    # 'variance'
    my @fs  = sort { $b <=> $a } values %variance;
    my $VAR = 0; map { $VAR += $fs[$_] } (1..$#fs);
    $VAR    = 1 if ($VAR > 1);
    # warn "$VAR = ".join(',', @fs)."\n";

    # PopA = 100 G / PopB = 100A
    #   $MAF = 0, $VAR = 1
    # PopA = 50 A 50 G / PopB = 75 A 25 C
    #   $MAF = 0.5, $VAR = 0.75
    
    my ($scb, $scNorm) = $args->{SCALEMAF} ? ($MAF, 0.5) : ($VAR, 1);
    my $scale = int(0.5 + 100 * $scb / $scNorm) / 100;
    if ($scale > 1) {
        $scale = 1;
    } elsif ($scale < 0.1) {
        $scale = 0.1;
    }

    my $cats = $self->each_category();
    $self->read_tags();
    my $rv = {
        id       => $self->best_accession(),
        chrid    => $self->full_loc(),
        accs     => $self->accessions_with_authorities(),
        auth     => [ $self->each_authority() ],
        data     => $dat,
        cats     => $cats,
        type     => $type,
        MAF      => int(0.5 + 1000 * $MAF) / 10,
        # VAR      => int(0.5 + 1000 * $VAR) / 10,
        w        => $wid,
        hideName => 1,
        impact   => $imp,
        scaleHeight => $scale,
        freqs    => $pop,
        # dist     => $self->each_population_distribution(),
        tags     => $self->all_tag_values(),
    };
    if (my $pkey = $self->pkey()) {
        $rv->{loc_id} = $pkey;
    } else {
        $rv->{handle} = $self->handle();
    }
    if (defined $cl) {
        # We have remapped the original coordinates
        $rv->{chrdata} = [[$cl, $cr]];
    }
    $rv->{insertion} = 1 unless ($wid);
    my $oi = $args->{OVIDATA} || $self->overall_impact
        ( -impdata => $imp, -popdata => $pop, @_ );
    while (my ($k, $v) = each %{$oi}) {
        $rv->{$k} = $v;
    }
    $self->bench_end();
    return $rv;
}

sub overall_impact {
    my $self = shift;
    $self->bench_start();
    my $args    = $self->parseparams( @_ );
    my $imp     = $args->{IMPDATA} || $self->impact( @_ );
    my @impHash = values %{$imp};
    my $ml      = $self->maploc();
    my $impReq  = $args->{IMPACT};
    my $useBad  = $args->{USEBAD};
    my $pop     = $args->{POPDATA};
    unless ($pop) {
        my @pbits;
        foreach my $key ('TOSSPOP', 'KEEPPOP') {
            if (my $val = $args->{$key}) {
                push @pbits, ($key, $val);
            }
        }
        $pop = $self->each_population_id( @pbits );
    }
    my (%alleles, %imps);
    # If there are no RNAs observed set impact to Genomic:
    $imps{GEN} = 1 if ($#impHash == -1);
    my $rv = {};
    foreach my $ih (@impHash) {
        my $hb  = $ih->{rnaHB};
        # Do not consider suboptimal alignments
        if ($hb && !$useBad) {
            next;
        }
        my $racc = "Unk";
        if (my $rid = $ih->{rid}) {
            my $rna = $ml->get_rna_by_id( $rid );
            $racc   = $rna->acc();
        }
        next if ($impReq && $racc ne $impReq);
        if (my $im = $ih->{imp}) { 
            map { $imps{$_} = 1 } split(/\//, $im);
        }
        while (my ($fid, $srcs) = each %{$ih->{features} || {}}) {
            push @{$rv->{features}{$fid}}, @{$srcs};
        }
        if (exists $ih->{var}) { 
            map { $alleles{$_}++ } keys %{$ih->{var}};
        }
    }
    if ($pop) {
        foreach my $ih (values %{$pop}) {
            map { $alleles{$_}++ } keys %{$ih};
        }
    }
    my @als = sort keys %alleles;
    $rv->{alleles} = [ \@als ] unless ($#als == -1);

    my ($rd) = $self->maploc->heaviest_impact( keys %imps );
    $rv->{outline}  = $rd->{color};
    $rv->{impName}  = $rd->{name};
    $rv->{impToken} = $rd->{token};

    while (my ($fid, $srcs) = each %{$rv->{features} || {}}) {
        my %u = map { $_ => 1 } @{$srcs};
        $rv->{features}{$fid} = [ keys %u ];
    }
    $self->bench_end();
    return $rv;
}

# ::Location
sub impact {
    my $self = shift;
    $self->bench_start();
    my $args     = $self->parseparams( @_ );
    my $ml       = $self->maploc();
    my $flank    = 500;
    my $locHsp   = $self->hsp();
    my $anc      = $self->anchor_id;
    my $rnaFilt  = $args->{RNAFILTER};
    my $nearRNAs = $self->nearby_rnas( @_ );
    my $minDist  = $#{$nearRNAs} == -1 ? 0 : $nearRNAs->[0][2];
    my %imps;
    foreach my $rdat (@{$nearRNAs}) {
        my ($rId, $alnId, $dist) = @{$rdat};
        last if (abs($dist) > abs($minDist));
        my $rna  = $ml->get_rna_by_id($rId);
        my $aln  = $ml->get_alignment( $alnId );
        if ($rnaFilt && &{$rnaFilt}( $rna, $aln ) ) {
            next;
        }
        my $imp = $rna->impact( @_, -snp => $self, -align => $aln );
        push @{$imps{$rId}}, $imp;
    }
    my $nullTag = '-NULL-';
    # $ml->prebranch(\%imps) if ($self->pkey() == 44030822);
    foreach my $rId (keys %imps) {
        # What if there are two alignments for the same RNA that overlap?
        # This is possible for genes in regions of tandem duplication
        my @sets;
        foreach my $set (sort { $a->{rnaHB} <=> $b->{rnaHB} } @{$imps{$rId}}) {
            # Take the best alignment of an RNA at this position
            last if ($#sets != -1 && $set->{rnaHB} &&
                     $set->{rnaHB} > $sets[0]{rnaHB});
            push @sets, $set;
        }
           
        if ($#sets == 0) {
            # We expect this in the vast majority of cases
            $imps{$rId} = $sets[0];
            next;
        }
        # Need to combine sets
        my %comb;
        # map { $comb{$_} ||= $_ eq 'var' ? {} : [] } map { keys %{$_} } @sets;
        my %isArray;
        # my @cols = keys %comb;
        foreach my $set (@sets) {
            foreach my $col (keys %{$set}) {
                my $val = $set->{$col};
                my $r   = ref($val);
                if (!$r) {
                    $val = $nullTag if (!defined $val);
                    push @{$comb{$col}}, $val;
                } elsif ($r eq 'HASH') {
                    while (my ($v, $row) = each %{$val}) {
                        push @{$comb{$col}{$v}}, $row;
                    }
                } elsif ($r eq 'ARRAY') {
                    push @{$comb{$col}}, @{$val};
                    $isArray{$col} = 1;
                }
            }
        }
        foreach my $col (keys %comb) {
            my $val = $comb{$col};
            my $r   = ref($val);
            if ($col eq 'var') {
                while (my ($v, $rows) = each %{$val}) {
                    my @src = @{$comb{$col}{$v}};
                    my $row = $comb{$col}{$v} = [];
                    for my $i (0..3) {
                        my @d = map { defined $_ ? $_ : '?' } map { $_->[$i] } @src;
                        my %u = map { $_ => 1 } @d;
                        my @k = keys %u;
                        $row->[$i] = $#k == 0 ? $d[0] : join('/', @d);
                    }
                }
            } elsif ($col eq 'features') {
                my $t = $comb{$col} = {};
                while (my ($v, $hashArr) = each %{$val}) {
                    my %vals;
                    foreach my $h (@{$hashArr}) {
                        while (my ($hk, $hv) = each %{$h}) {
                            $hv = "" unless (defined $hv);
                            $vals{$hk}{$hv} = 1;
                        }
                    }
                    while (my ($hk, $hvs) = each %vals) {
                        delete $hvs->{''};
                        $t->{$v}{$hk} = join('/', sort keys %{$hvs});
                    }
                }
            } else {
                my %u = map { $_ => 1 } @{$val};
                my @k =  keys %u;
                if ($isArray{$col}) {
                    $comb{$col} = [];
                    foreach my $v (sort @k) {
                        push @{$comb{$col}}, $v eq $nullTag ? undef : $v;
                    }
                } else {
                    my $v = $#k == 0 ? $val->[0] : join('/', @{$val});
                    $comb{$col} = $v eq $nullTag ? undef : $v;
                }
            }
        }
        $imps{$rId} = \%comb;

    }
    $self->bench_end();
    return \%imps;
}

# ::Location
sub to_one_line {
    my $self = shift;
    my @accs = $self->each_accession();
    return sprintf("%s [%d] %s", $self->full_loc, $self->handle,
                   join(',', @accs) || "");
}

*position = \&loc_text;
*pos = \&loc_text;
sub loc_text {
    my $self = shift;
    unless ($self->{LOC_TEXT}) {
        $self->bench_start();
        my ($l, $r) = ($self->left, $self->right);
        if ($r == $l + 1) {
            # Insertion
            $self->{LOC_TEXT} = "$l^$r";
        } else {
            my ($s, $e) = ($l + 1, $r -1);
            if ($s == $e) {
                # Single position
                $self->{LOC_TEXT} = $s;
            } else {
                # Range
                $self->{LOC_TEXT} = "$s..$e";
            }
        }
        $self->bench_end();
    }
    return $self->{LOC_TEXT};
}

# :: Location
sub name {
    return shift->best_accession();
}

*description = \&desc;
sub desc {
    my $self = shift;
    my $desc = $self->full_loc();
    my @cats = $self->each_category();
    $desc .= " = ".join(', ', @cats) unless ($#cats == -1);
    return $desc;
}

# ::Location
*full_location = \&full_loc;
sub full_loc {
    my $self = shift;
    return sprintf("%s.%s:%s", $self->chr, $self->build, $self->loc_text);
}

our $locPkeyCols = [ qw(build chr pos_to_left pos_to_right) ];
*handle  = \&pkey;
sub pkey {
    my $self = shift;
    unless (defined $self->{PKEY}) {
        $self->bench_start();
        my @args = map { $self->{$_} } qw(BLD CHR LFT RGT);
        my $len = $args[3] - $args[2] - 1;
        if ($len < 0) {
            $self->err("Attempt to set negative width location",
                       join(" + ", @args). " = $len");
            $self->{PKEY} = 0;
        } else {
            my $getSet = $self->maploc->{STH}{"GET_LOC_PKEY_FUNC"} ||= 
                $self->maploc->dbh->prepare
                ( -name   => "Get or create location pkey",
                  -sql    => "SELECT location_pkey(?,?,?,?)", );
            # 410 us for existing data, 2 ms for new insertion:
            unless ($self->{PKEY} = $getSet->get_single_value( @args)) {
                $self->death("Failed to generate or recover pkey for location",
                             join(',', @args));
            }

            # 480 us for existing data:
            #$self->{PKEY} = $self->_generic_get_pkey
            #    (\@args, 'location', 'loc_id', $locPkeyCols );
        }
        $self->bench_end();
    }
    return $self->{PKEY};
}

sub loc_acc_id {
    my $self = shift;
    unless (defined $self->{LOC_ACC_ID}) {
        # Benchmark ~ 32.32us
        $self->{LOC_ACC_ID} = $self->maploc->loc_to_accid
            ($self->build, $self->chr);
    }
    return $self->{LOC_ACC_ID};
}

sub location_accession {
    my $self = shift;
    unless (defined $self->{LOC_ACC}) {
        $self->bench_start();
        $self->{LOC_ACC} = $self->maploc->loc_to_acc
            ($self->build, $self->chr);
        $self->bench_end();
    }
    return $self->{LOC_ACC};
}

sub distance_to_alignment {
    my $self = shift;
    my $aln  = shift;
    my $ind  = $aln->seqindex( $self->acc );
    return undef unless (defined $ind);

}

sub nearby_locations {
    my $self = shift;
    my $slop = shift || 0;
    my @rv = $self->maploc->overlapping_locations
        ($self->lft - $slop, $self->rgt + $slop, $self->chr, $self->build);
    return @rv;
}

sub nearby_alignment_ids {
    my $self = shift;
    my $slop = shift || 0;
    return $self->maploc->overlapping_alignments_by_id
        ($self->lft - $slop, $self->rgt + $slop, $self->loc_acc_id );
}

sub nearby_alignments {
    my $self = shift;
    $self->bench_start();
    my $arows = $self->nearby_alignment_ids( @_ );
    my $ml    = $self->maploc;
    my @rv = map { $ml->get_alignment( $_->[0] ) } @{$arows};
    my $accid = $self->location_accession();
    map { $_->set_anchor( $accid ) } @rv;
    $self->bench_end();
    return wantarray ? @rv : \@rv;
}

# ::Location
sub nearby_rnas {
    my $self     = shift;
    my $args     = $self->parseparams( @_ );
    my $ml       = $self->maploc();
    my $howbad   = $args->{HOWBAD};
    my $flank    = $args->{RANGE} || $args->{FLANK};
    my $defFlank = $ml->default_loc_to_rna_distance();
    my $defHB    = $ml->default_howbad();
    my $allVers  = 0; # $args->{ALLVERS} || $args->{ALLVERSIONS};
    my $isDef    = (defined $howbad && $howbad != $defHB) ||
        (defined $flank && $flank != $defFlank) ? 0 : 1;
    my $rv       = $isDef ? $self->{DEFAULT_NEAR_RNAS} : undef;
    unless ($rv) {
        $flank  = $defFlank unless (defined $flank);
        $howbad = $defHB unless (defined $howbad);
        my $anc = $self->anchor_id;
        my $rnaRows;
        if ($isDef) {
            my $bulkSrc = $ml->bulk_cache_source( 'nearby_rnas' );
            $rnaRows    = $bulkSrc->{$self->handle()};
            # $ml->prebranch($rnaRows) if ($self->pkey == 44244280);
        }
        my ($l, $r) = ($self->lft, $self->rgt);
        $rnaRows ||= $ml->overlapping_rnas_by_id
            ( $l - $flank, $r + $flank, $anc, $howbad );
        my %rnas;
        foreach my $rr (@{$rnaRows}) {
            my ($alnId, $seqId, $sc, $ptl, $ptr, $rId) = @{$rr};
            # Note how far the alignment is from the location
            # Zero indicates overlap
            # abs($dist) == 1 -> Adjacent
            my $dist = 0;
            if ($ptl - $r > 2) {
                # Location is 'to the left'
                # Positive distance
                $dist = $ptl - $r - 2;
            } elsif ($l - $ptr > 2)  {
                # Location is 'to the right'
                # Negative distance
                $dist = 0 - ($l - $ptr - 2);
            }
            my $accV = $ml->cached_pkey_to_text( $seqId );
            my $accU = $accV;
            my $vers = 0;
            if (!$allVers && $accV =~ /(.+)\.([^\.]+)$/) {
                ($accU, $vers) = ($1, $2);
            }
            push @{$rnas{$accU}{$vers}}, [ $rId, $alnId, $dist ];
        }
        my @found;
        while (my ($accU, $vH) = each %rnas) {
            my ($vers) = sort { $b <=> $a } keys %{$vH};
            # rna_id / aln_id / distance
            push @found, @{$rnas{$accU}{$vers}};
        }
        $rv = [ sort { abs($a->[2]) - abs($b->[2]) ||
                           $a->[0] <=> $b->[0] } @found ];
        $self->{DEFAULT_NEAR_RNAS} = $rv if ($isDef);
    }
    return $rv;
}

sub closest_rnas {
    my $self = shift;
    return $self->maploc->closest_rnas
        ( -lft   => $self->lft,
          -right => $self->rgt,
          -accid => $self->loc_acc_id,
          @_);
}

# ::Location
sub each_rna_id {
    my $self = shift;
    my $data = $self->nearby_rnas( @_ );
    # $self->maploc->prebranch(\@data);
    return map { $_->[0] } @{$data || []};
}

# ::Location
sub each_rna {
    my $self = shift;
    my @ids  = $self->each_rna_id( @_ );
    my $ml   = $self->maploc();
    return map { $ml->get_rna_by_id( $_ ) } @ids;
}

sub lock_alleles {

    # This just prevents alleles from being read from the database
    # This is useful if alleles are being manually provided via
    # extra_alleles(), and the user does not want other alleles to be
    # read in.

    my $self = shift;
    $self->{ALLELES} = {};
}

sub each_allele_id {
    my $self = shift;
    my $rv = $self->{ALLELE_IDS};
    unless ($rv) {
        $self->bench_start();
        $rv      = $self->{ALLELE_IDS} = {};
        my $ml   = $self->maploc();
        my $bcs  = $ml->bulk_cache_source( 'allele' );
        my $pk   = $self->pkey();
        my $rows = $bcs->{$pk};
        if (!$rows && $pk) {
            my $get = $ml->{STH}{ALLELES_FOR_LOC} ||= $ml->dbh->prepare
                ( -name => "Get alleles for location",
                  -sql  => "SELECT base_id, pop_id, freq, count ".
                  "FROM allele WHERE loc_id = ?");
            $get->execute( $pk );
            $rows = $get->fetchall_arrayref();
        }
        foreach my $row (@{$rows || []}) {
            my ($bid, $pid, $f, $c) = @{$row};
            $rv->{$bid}{$pid} = [$f,$c];
        }
        $self->bench_end();
    }
    return wantarray ? sort keys %{$rv} : $rv;
}

*alleles = \&each_allele;
sub each_allele {
    my $self = shift;
    my $rv = $self->{ALLELES};
    unless ($rv) {
        $self->bench_start();
        $rv      = $self->{ALLELES} = {};
        my $alli = $self->each_allele_id();
        my $ml   = $self->maploc();
        while (my ($bid, $pHash) = each %{$alli}) {
            my $base = $ml->cached_pkey_to_text( $bid );
            $rv->{$base} = $pHash;
        }
        $self->bench_end();
    }
    if ($#_ == -1 && !$self->{XTRA_ALLELE}) {
        # Just return everything
        return wantarray ? sort keys %{$rv} : $rv;
    }
    if (my $xa = $self->{XTRA_ALLELE}) {
        # We have to add in extra alleles
        my %xtra = %{$rv};
        while (my ($bid, $pHash) = each %{$xa}) {
            if ($xtra{$bid}) {
                my $targ = $xtra{$bid} = { %{$xtra{$bid}} };
                while (my ($pid, $dat) = each %{$pHash}) {
                    $targ->{$pid} = $dat;
                }
            } else {
                $xtra{$bid} = $pHash;
            }
        }
        $rv = \%xtra;
    }
    my $args = $self->parseparams( @_ );
    my $minFreq = $args->{MINFREQ};
    my $cid     = $args->{CATID};
    if ($minFreq || $cid) {
        # Filter alleles to a minimum frequency and/or category membership
        my ($cSth, %tossPid);
        if ($minFreq) {
            my $mafs     = $self->population_maf();
            my $keepNull = $args->{KEEPNULL};
            if (ref($minFreq)) {
                # The user has passed a hash reference for fine-grained
                # control of frequency filtering. It should be keyed
                # to PID and point to a frequency. Key 0 can be used for
                # a default minimum frequency
                while (my ($pid, $maf) = each %{$mafs}) {
                    if (!defined $maf) {
                        $tossPid{$pid} = 1 unless ($keepNull);
                    } else {
                        unless (defined $minFreq->{$pid}) {
                            # No explicit frequency for this specific pop
                            # See if we have a frequency defined for a parent
                            # population or category
                            my $useFreq = 999;
                            my $ml      = $self->maploc();
                            my @pkeys = $ml->all_parent_pop_ids($pid);
                            push @pkeys, $ml->all_parent_cat_ids($pid);
                            #warn "$pid -> ".join(',', @pkeys)."\n";
                            foreach my $id ( @pkeys ) {
                                # We will keep the most permisive freq we find
                                if (exists $minFreq->{$id} &&
                                    defined $minFreq->{$id} &&
                                    $useFreq > $minFreq->{$id}) {
                                    $useFreq = $minFreq->{$id};
                                }
                            }
                            # Use parent pop freq if found, or the default
                            # freq, or set as an always permissive -1
                            $minFreq->{$pid} = $useFreq < 999 ? $useFreq :
                                defined $minFreq->{0} ? $minFreq->{0} : -1;
                        }
                        if ($maf < $minFreq->{$pid}) {
                            # This population failed the frequency filters
                            $tossPid{$pid} = 1;
                        }
                    }
                }
            } else {
                # There is a global frequency filter defined
                # Simple, but can lead to confusion with users
                while (my ($pid, $maf) = each %{$mafs}) {
                    if (!defined $maf) {
                        $tossPid{$pid} = 1 unless ($keepNull);
                    } elsif ($maf < $minFreq) {
                        $tossPid{$pid} = 1;
                    }
                }
            }
        }
        if ($cid) {
            $cid  = [$cid] unless (ref($cid));
            $cSth = $self->_fast_category_filter( $#{$cid} );
        }
        
        my %loc;
        while (my ($base, $pH) = each %{$rv}) {
            while (my ($pid, $dat) = each %{$pH}) {
                next if ($tossPid{$pid});
                if ($cid) {
                    $cSth->execute($pid, @{$cid});
                    my $rows =  $cSth->fetchall_arrayref();
                    next if ($#{$rows} == -1);
                }
                $loc{$base}{$pid} = $dat;
            }
        }
        $rv = \%loc;
    }

    if (my $filter = $args->{TOSSPOP}) {
        # Request to exclude alleles only associated with one or more populations
        my %keep;
        while (my ($base, $pH) = each %{$rv}) {
            my %pLoc = %{$pH};
            map { delete $pLoc{$_} } @{$filter};
            my @left   = keys %pLoc;
            $keep{$base} = \%pLoc unless ($#left == -1);
        }
        $rv = \%keep;
    }
    if (my $filter = $args->{KEEPPOP}) {
        # Request to only keep alleles associated with one or more populations
        my %keep;
        while (my ($base, $pH) = each %{$rv}) {
            my $found = 0; map { $found++ if (exists $pH->{$_}) } @{$filter};
            $keep{$base} = $pH if ($found);
        }
        $rv = \%keep;
    }
    if ($args->{REVCOM} || $args->{REVCOMP}) {
        # Request to reverse complement the bases
        my %rc;
        my $su = $self->maploc->sequtils();
        while (my ($base, $pH) = each %{$rv}) {
            if (my $rev = $su->revcom($base)) {
                $rc{ $rev } = $pH;
            }
        }
        $rv = \%rc;
    }
    
    return wantarray ? sort keys %{$rv} : $rv;
}

sub extra_alleles {
    my $self = shift;
    if (my $nv = shift) {
        my @chk = ref($nv) ? @{$nv} : ($nv);
        my $ml   = $self->maploc();
        my $tpid = $ml->_temporary_pop_id();
        my $str  = shift || 0;
        foreach my $req (@chk) {
            next unless ($req);
            $req = uc($req);
            my @bases;
            if ($req =~ /^([^>]+)>([^>]+)$/) {
                # A>T, implies that A is going to T
                # We will put a small non-zero frequency on the first allele
                my $base = $1;
                $req = $2;
                if ($base =~ /^[ACTG]+$/ || $base eq '-') {
                    $base = $ml->revcom($base) if ($str < 0);
                    $self->{XTRA_ALLELE}{$base}{ $tpid } ||= [$tinyFakeFreq];
                } else {
                    $self->msg_once("[!]","Can not interpret allele '$base'");
                }
            }
            @bases = split(/[\s\,\/>]+/, $req);
            foreach my $base (split(/[\s\,\/>]+/, uc($req))) {
                if ($base =~ /^[ACTG]+$/ || $base eq '-') {
                    $base = $ml->revcom($base) if ($str < 0);
                    $self->{XTRA_ALLELE}{$base}{ $tpid } ||= [];
                } else {
                    $self->msg_once("[!]","Can not interpret allele '$base'");
                }
            }
        }
    }
    return wantarray ? keys %{$self->{XTRA_ALLELE} || {}} :
        $self->{XTRA_ALLELE};
}

sub clear_extra_alleles {
    my $self = shift;
    delete $self->{XTRA_ALLELE};
}

our $minorAlleleRound = 10000;
sub population_maf {
    my $self = shift;
    my $rv   = $self->{POP_MAF};
    unless ($rv) {
        my $ml   = $self->maploc();
        my $alls = $self->each_allele_id();
        my %pbHash;
        while (my ($bid, $pidH) = each %{$alls}) {
            while (my ($pid, $dat) = each %{$pidH}) {
                $pbHash{$pid}{$bid} = defined $dat->[0] ? $dat->[0] : -1;
            }
        }
        $rv = $self->{POP_MAF} = 
            $self->maploc->_maf_from_popBaseHash( \%pbHash );
        
        #my $pids = $self->each_population_id();
        #$self->maploc->prebranch($pids);
        #while (my ($pid, $allH) = each %{$pids}) {
        #    my @alls;
        #    while (my ($base, $dat) = each %{$allH}) {
        #        if (my $f = $dat->[0]) {
        #            push @alls, $f;
        #        }
        #    }
        #    if ($#alls == -1) {
        #        $rv->{$pid} = undef;
        #    } else {
        #        @alls = sort { $b <=> $a } @alls;
        #        $rv->{$pid} = int(0.5 + $minorAlleleRound * (1 - $alls[0])) 
        #            / $minorAlleleRound;
        #    }
        #}
    }
    return $rv if ($#_ == -1);
    # Filter requests have been passed ...
    my $args = $self->parseparams( @_ );
    # Make a local copy
    $rv = { %{$rv} };
    if (my $cid  = $args->{CATID}) {
        $cid     = [$cid] unless (ref($cid));
        my $cSth = $self->_fast_category_filter( $#{$cid} );
        foreach my $pid (keys %{$rv}) {
            $cSth->execute($pid, @{$cid});
            my $rows =  $cSth->fetchall_arrayref();
            delete $rv->{$pid} if ($#{$rows} == -1);
        }
    }
    return $rv;
}

sub _fast_category_filter {
    my $self = shift;
    my ($ind) = shift || 0;
    unless ($self->{MAPLOC}{"FastCategoryFilter$ind"}) {
        my $ml   = $self->maploc();
        my $dbh  = $ml->dbh();
        my $tid  = $ml->text_to_pkey("Category");
        $self->{MAPLOC}{"FastCategoryFilter$ind"} = $dbh->prepare
            ( -name => "Get 'Category' tag values for a location",
              -sql  => "SELECT tv.val_id FROM tagval tv WHERE tv.obj_id = ? AND tv.tag_id = $tid AND tv.val_id IN (".join(',', map {'?' } (0..$ind)).")",
              -limit => 1);
    }
    return $self->{MAPLOC}{"FastCategoryFilter$ind"};
}

sub OLD_minor_allele_representation {
    my $self = shift;
    my $pids = $self->each_population_id();
    my $args = $self->parseparams( @_ );
    my $cid  = $args->{CATID};
    my $cSth;
    
    if ($cid) {
        # Restrict 
        my $ml  = $self->maploc();
        $cSth   = $ml->{STH}{FastCategoryFilter};
        $cid    = [$cid] unless (ref($cid));
        unless ($cSth) {
            my $dbh = $ml->dbh();
            my $tid = $ml->text_to_pkey("Category");
            $cSth   = $ml->{STH}{FastCategoryFilter} = $dbh->prepare
                ( -name => "Get 'Category' tag values for a location",
                  -sql  => "SELECT tv.val_id FROM tagval tv WHERE tv.obj_id = ? AND tv.tag_id = $tid AND tv.val_id IN (".join(',', map {'?' } @{$cid}).")",
                  -limit => 1);
        }
    }
    my (%byAllele, %maxAllele);
    while (my ($pid, $allH) = each %{$pids}) {
        if ($cSth) {
            $cSth->execute($pid, @{$cid});
            my $rows =  $cSth->fetchall_arrayref();
            next if ($#{$rows} == -1);
        }
        my @alls;
        while (my ($base, $dat) = each %{$allH}) {
            my $f = $dat->[0];
            push @alls, [ $f, $base ] if ($f);
        }
        next if ($#alls == -1);
        @alls = sort { $b->[0] <=> $a->[0] } @alls;
        # There is at least one allele in this population with nonzero freq
        # We need to track these in case there are two or more populations
        # In such a case we could have:
        # Pop:1 100% A, Pop:2 100% C
        # Neither alone has a minor (secondary) allele, but combined they
        # indicate 'high' variability at the position
        my ($f, $base) = @{$alls[0]};
        if (!$maxAllele{$base} || $maxAllele{$base}[0] < $f) {
            $maxAllele{$base} = [ $f, $base, [$pid] ];
        } elsif ($maxAllele{$base}[0] == $f) {
            push @{$maxAllele{$base}[2]}, $pid;
        }
        next unless ($#alls > 0);
        # We have identified the second-most common allele for this pop
        # We will actually report the sum of all minor alleles, not just
        # the second-most common one
        $f = int(0.5 + $minorAlleleRound * (1 - $f)) / $minorAlleleRound;
        push @{$byAllele{$alls[1][1]}}, [$f, $pid] if ($f);
    }
    my @rv;
    while (my ($base, $datArr) = each %byAllele) {
        my @dats  = sort { $b->[0] <=> $a->[0] } @{$datArr};
        my $bestF = $dats[0][0];
        my @pids;
        while ($#dats > -1) {
            my $dat = shift @dats;
            my ($f, $pid) = @{$dat};
            last if ($f < $bestF);
            push @pids, $pid;
        }
        push @rv, [$bestF, $base, \@pids];
    }
    if ($#rv == -1) {
        # We will only consider the maxallele data if we do not find
        # a secondary allele the 'normal' way.
        @rv = sort { $b->[0] <=> $a->[0] } values %maxAllele;
        shift @rv if ($#rv <= 0 || $rv[0][0] > $rv[1][0]);
        map { $_->[0] = int(0.5 + $minorAlleleRound * $_->[0])
                  / $minorAlleleRound } @rv;
        # I am not sure this is the correct way to approach the problem.
        # If there are three populations:
        # Pop:1 95% A 5% T
        # Pop:2 100% C
        # Pop:3 100% G
        # ... ? 
    }
    # warn $self->branch({pids => $pids, ba => \%byAllele, rv => \@rv});
    @rv = sort { $b->[0] <=> $a->[0] } @rv;
    return wantarray ? @rv : \@rv;
}

sub each_population_id {
    my $self = shift;
    my $rv = $self->{POPULATION_IDS};
    unless ($rv) {
        # This is just a pivot of allele data
        $rv = $self->{POPULATION_IDS} = {};
        my $alls = $self->each_allele( );
        while (my ($base, $pidH) = each %{$alls}) {
            while (my ($pid, $dat) = each %{$pidH}) {
                $rv->{$pid}{$base} = $dat;
            }
        }
    }
    unless ($#_ == -1) {
        # User is passing arguments
        my $args = $self->parseparams( @_ );
        if (my $filter = $args->{TOSSPOP}) {
            # Request to exclude certain populations
            $rv = { %{$rv} };
            map { delete $rv->{$_} } @{$filter};
        }
        if (my $filter = $args->{KEEPPOP}) {
            # Request to only keep certain populations
            my %keep;
            foreach my $pid (@{$filter}) {
                if (exists $rv->{$pid}) {
                    $keep{$pid} = $rv->{$pid};
                }
            }
            $rv = \%keep;
        }
    }
    return wantarray ? sort keys %{$rv} : $rv;
}

sub each_population_distribution {
    my $self = shift;
    my $rv = $self->{POPULATION_DISTRIB};
    unless ($rv) {
        my $alls = $self->each_allele( );
        my %byPop;
        while (my ($base, $pidH) = each %{$alls}) {
            while (my ($pid, $dat) = each %{$pidH}) {
                my $f = $dat->[0];
                $f = 1 unless (defined $f);
                $byPop{$pid}{$base} = $f;
            }
        }
        $rv = $self->{POPULATION_DISTRIB} = {};
        while (my ($pid, $bdat) = each %byPop) {
            my @als = sort keys %{$bdat};
            my $key = join('/', @als);
            my $targ = $rv->{$key} ||= {
                key     => $key,
                alleles => \@als,
                pop     => [],
            };
            push @{$targ->{pop}}, [ $pid, map { $bdat->{$_} } @als ];
        }
    }
    return $rv;
}

sub each_population {
    my $self = shift;
    my $rv = $self->{POPULATIONS};
    unless ($rv) {
        # This is just a pivot of allele data
        my $pops = $self->each_population_id( );
        my $ml   = $self->maploc();
        my @unsort;
        while (my ($pid, $baseH) = each %{$pops}) {
            my $pop = $ml->get_population( $pid );
            my $json = $pop->json_data();
            push @unsort, $json->{name};
        }
        $self->{POPULATIONS} = [ sort @unsort ];
    }
    return sort @{$rv};
}

sub clear_alleles_for_population {
    my $self = shift;
    return if ($self->{READ_ONLY});
    my $pop  = shift;
    return unless ($pop);
    $self->bench_start();
    $pop = $pop->id if (ref($pop));
    my $ml    = $self->maploc();
    my $clear = $ml->{STH}{"CLEAR_"} ||= $ml->dbh->prepare
        ( -name => "Delete all al",
          -sql  => "FOO",
          -ignore => "duplicate key value" );

    die "WORKING";

    $self->bench_end();
}

# Check for misbehavior:
# SELECT * FROM allele WHERE freq > 1 OR freq < 0;
sub add_allele_to_db {
    my $self = shift;
    return "" if ($self->{READ_ONLY});
    # $self->bench_start();
    # Average time ~ 410 us
    # A lot of loading time will occur here, I have commented out benchmarking
    # since the benchmarks themselves were taking significant time

    my $base = shift;
    unless ($base) {
        #$self->bench_end();
        return 'No base provided';
    }
    my $pkey = $self->pkey();
    unless ($pkey) {
        # $self->bench_end();
        return 'No PKEY for location';
    }
    # Upcase all bases that are pure nucleotides:
    my $ml     = $self->maploc();
    $base      = uc($base) if ($base =~ /^[A-Z]$/i);
    my $bid    = $ml->cached_text_to_pkey( $base );
    my $pop_id = shift || $ml->_unknown_pop_id();
    my $rv     = 0;
    my ($freq, $denom, $defer) = @_;
    # If the base row already exists, it is faster to check first, then try to 
    # add the new row only if needed:
    # 290-300us
    if (!defined $freq) {
        $denom = undef;
    } elsif ($freq > 1 || $freq < 0) {
        $rv = "Frequency '$freq".($denom ? " / $denom" : "").
            "' is not between 0 and 1";
        $self->msg("[?]", $rv, $self->to_one_line());
        $freq = undef;
        $denom = undef;
    }
    # Testing query performance under different circumstances:
    # $freq = 0.222; $denom = 222; $self->msg_once("[!!]","MALFORMING DATA to ". join(" + ", map { defined $_ ? $_ : 'NULL'} ($freq, $denom)));

    # my $data = [$pkey, $bid, $pop_id, $freq, $denom];
    my $sth = $ml->{STH}{"CREATE_ALLELE_BY_FUNCTION"} ||= $ml->dbh->prepare
        ( -name => "Create an allele row with function call",
          -sql  => "SELECT add_allele_to_db(?,?,?,?,?,?)",
          -ignore => "duplicate key value" );

    #$self->bench_start('execute');
    my $addRv = $sth->get_single_value
        ($pkey, $bid, $pop_id, $freq, $denom, $defer ? 1 : 0);
    #$self->bench_end('execute');
    $ml->{ALLELE_UPDATE_STATUS}{$addRv || 'UNDEF'}++;
    unless ($addRv) {
        # The database has no data, and we are waiting to load in bulk
        $ml->defer_alleles_to_db( [$pkey, $bid, $pop_id, $freq, $denom] );
    }
    # $self->bench_end();
    return $rv;
}

sub clear_alleles_from_db {
    my $self = shift;
    return if ($self->{READ_ONLY});
    my $pkey = $self->pkey();
    $self->death("Can not clear allele unless location has pkey",
                 $self->to_one_line()) unless ($pkey);
    $self->bench_start();
    my $ml     = $self->maploc();
    my $pop_id = shift || $ml->_unknown_pop_id();
    my @rv;
    if (my $except = shift) {
        # We want to clear all alleles except specific ones
        # The code will return the alleles that were removed, so we
        # need to iterate them
        my %keepid;
        foreach my $req (ref($except) ? @{$except} : ($except)) {
            unless ($req =~ /^\d+$/) {
                $req = $ml->text_to_pkey( $req );
            }
            $keepid{$req} = 1;
        }
        my $find = $ml->{STH}{"FIND_ALL_POP_BASES"} ||= $ml->dbh->prepare
            ( -name => "Find all alleles for a population",
              -sql  => "SELECT base_id FROM allele WHERE ".
              "loc_id = ? AND pop_id = ?");
        my @found = $find->get_array_for_field($pkey, $pop_id);
        my $kill  = $ml->{STH}{"CLEAR_SINGLE_POP_BASE"} ||= $ml->dbh->prepare
            ( -name => "Clear specific allele for a population",
              -sql  => "DELETE FROM allele WHERE ".
              "loc_id = ? AND pop_id = ? AND base_id = ?");
        foreach my $bid (@found) {
            unless ($keepid{$bid}) {
                $kill->execute($pkey, $pop_id, $bid);
                push @rv, $bid;
            }
        }
    } else {
        # Just clear all alleles
        my $sth = $ml->{STH}{"CLEAR_ALL_POP_BASES"} ||= $ml->dbh->prepare
            ( -name => "Clear all alleles for a population",
              -sql  => "DELETE FROM allele WHERE ".
              "loc_id = ? AND pop_id = ?");
        $sth->execute($pkey, $pop_id);
    }
    return @rv;
}

sub read {
    my $self = shift;
    $self->read_tags();
    $self->read_categories();
}

# ::Location
*update = \&write;
sub write {
    my $self = shift;
    $self->write_tags();
    $self->write_categories();
}

*set_cat      = \&set_category;
*add_cat      = \&set_category;
*add_category = \&set_category;
sub set_category {
    my $self = shift;
    if (my $cat = shift) {
        if ($cat =~ /^\d+$/) {
            # Take this as a category ID
            $self->{CATS}{$cat} ||= undef;
        } else {
            my $id = $self->maploc->cached_text_to_pkey( $cat );
            $self->{CATS}{$id} ||= $cat;
        }
        delete $self->{CAT_NAMES};
    }
}

sub read_categories {
    my $self = shift;
    unless ($self->{DONE_READ_CATS}) {
        my $ml   = $self->maploc();
        my $sth  = $ml->{STH}{READ_LOCATION_CATEGORIES} ||= $ml->dbh->prepare
            ( -name => "Read categories for location",
              -sql  => "SELECT cat_id FROM loc_cat WHERE loc_id = ?");
        my @ids = $sth->get_array_for_field( $self->pkey() );
        map { $self->{CATS}{$_} = undef } @ids;
        $self->{DONE_READ_CATS} = 1;
        delete $self->{CAT_NAMES};
    }
}

sub clear_category {
    my $self = shift;
    if (my $cat = shift) {
        delete $self->{CATS}{$cat};
        delete $self->{CAT_NAMES};
        return if ($self->{READ_ONLY});
        my $ml  = $self->maploc();
        $cat    = $ml->cached_text_to_pkey( $cat ) unless ($cat =~ /^\d+$/);
        my $kill  = $ml->{STH}{"CLEAR_LOCATION_CATEGORY"} ||= $ml->dbh->prepare
            ( -name => "Delete a location category",
              -sql  => "DELETE FROM loc_cat WHERE loc_id = ? AND cat_id = ?" );
        $kill->execute( $self->pkey, $cat );
        return $cat;
    }
    return 0;
}

*all_cat_ids      = \&each_category_id;
*each_cat_id      = \&each_category_id;
*all_category_ids = \&each_category_id;
sub each_category_id {
    my $self = shift;
    return keys %{$self->{CATS}}
}

*all_cats       = \&each_category;
*each_cat       = \&each_category;
*all_categories = \&each_category;
sub each_category {
    my $self = shift;
    my $rv = $self->{CAT_NAMES};
    unless ($rv) {
        my @names;
        while (my ($cid, $cname) = each %{$self->{CATS}}) {
            push @names, $cname ||= $self->{CATS}{$cid} = 
                $self->maploc->cached_pkey_to_text( $cid );
        }
        $rv = $self->{CAT_NAMES} = [ sort @names ];
    }
    return wantarray ? @{$rv} : $rv;
}

sub write_categories {
    my $self = shift;
    return () if ($self->{READ_ONLY});
    my @ids  = $self->each_category_id();
    return () if ($#ids == -1);
    $self->bench_start();
    my $ml   = $self->maploc();
    # ~290 usec, but will spam log with 'violates unique constraint'
    #my $set  = $ml->{STH}{"WRITE_LOCATION_CATEGORIES"} ||= $ml->dbh->prepare
    #    ( -name => "Insert a location category",
    #      -sql  => "INSERT INTO loc_cat (loc_id, cat_id) VALUES (?,?)",
    #      -ignore => "duplicate key value" );

    # ~50 msec to write new category
    # 93 usec to ignore existing category entry
    my $set  =  $ml->{STH}{"ADD_CATEGORY_BY_FUNCTION"} ||= $ml->dbh->prepare
        ( -name   => "Add a category to database via function call",
          -sql    => "SELECT quiet_category_write(?,?)",
          -ignore => "duplicate key value" );
    # ~50 msec:
    # EXECUTE 'INSERT INTO loc_cat (loc_id, cat_id) VALUES ($1,$2)' 
    #   USING lid, cid;

    my $pkey = $self->pkey();
    map { $set->selectall_arrayref( $pkey, $_ ) } @ids;
    $self->bench_end();
    return @ids;
}

# ::Location
sub all_symbols {
    my $self = shift;
    unless ($self->{GENE_SYMS}) {
        my %syms;
        foreach my $gene ($self->each_gene('useCache')) {
            $gene->read();
            map { $syms{$_} = 1 } $gene->all_symbols();
        }
        $self->{GENE_SYMS} = [ sort keys %syms ];
    }
    my @rv = @{$self->{GENE_SYMS}};
    return wantarray ? @rv : \@rv;
}
