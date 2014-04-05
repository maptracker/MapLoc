#!/usr/bin/perl -w

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

my $isBeta;

use strict;
use BMS::ArgumentParser;
use BMS::SnpTracker::MapLoc;
use Bio::SeqIO;

my $args = BMS::ArgumentParser->new
    ( -limit    => 0,
      -nocgi    => 1,
      -fork     => 20,
      -verbose  => 1,
      -trial    => 1,
      -build    => 'GRCh37',
      -instance   => 'maploc',
      -progress => 120,
      -aspercent => 1,
      -binsize   => 10,
      );

my (%ignoredXREF, $only, $mapLoc, %textCache);

$args->shell_coloring();
$args->debug->max_any( 20 );


my $authority   = "Ensembl";
my $hasPtag     = "Has Protein";
my $sameAs      = "Same As";
my $dir         = $args->val(qw(directory dir));
my $isTrial     = $args->val(qw(testmode tm trial));
my $limit       = $args->val(qw(lim limit)) || 0;
unless ($dir) {
    &usage();
    $args->msg("[!]", "Please provide the directory path with -dir");
    exit;
}
$dir =~ s/\/$//;

# my $only;

if (my $oreq = $args->val(qw(only))) {
    my @reqs = ref($oreq) ? @{$oreq} : ($oreq);
    $only = {};
    foreach my $req (@reqs) {
        map { $only->{$_} = 1 } split(/\s*[\n\r\t\,]\s*/, uc($req || ""));
    }
    my @onlys = sort keys %{$only};
    if ($#onlys == -1) {
        $only = undef;
    } else {
        $args->msg('[-]',"Request to limit anlaysis to particular transcripts",
                   @onlys);
    }
}

my %tagMap = ( description   => 'Description',
               created_date  => 'Date Created',
               modified_date => 'Date Modified',
               status        => 'Ensembl Status',
               biotype       => 'Ensembl Biotype', );

$args->msg("[?]","FYI: Loading Human GRCh37 took ~3.5Gb of RAM");

my ($trialFile, $trialFH);
if ($isTrial) {
    $trialFile = "MapLoc-Ensembl_Load_Trial.txt";
    open($trialFH, ">$trialFile") || $args->death
        ("Failed to write trial file", $trialFile, $!);
    $args->msg("[!!]","----",
               "Trial mode, text summary of RNA and genes will be dumped to:",
               $trialFile,"To load database, pass the argument -trial 0",
               "----");
} else {
    $args->msg("[!!]","Live mode, database will be altered!");
}

my ($regLU);
my $ens_id_lookup = {};

my $gzFiles = &files();
my $schema  = &table_columns();
&get_build();

my %toDo = map { $_ => 1 } qw(REGION GENE RNA SEQ FEAT META);
foreach my $param (keys %toDo) {
    if ($args->val('no'.$param)) {
        # Do not do this task
        delete $toDo{$param};
    } elsif ($args->val('only'.$param, $param.'only')) {
        # Only do this task
        %toDo = ( $param => 1 );
    }
}

&region_lookup()    if ($toDo{REGION});
&genes()            if ($toDo{GENE});
&transcripts()      if ($toDo{RNA});
&rna_sequence()     if ($toDo{SEQ});
&protein_features() if ($toDo{FEAT});
&domain_meta()      if ($toDo{META});

$args->msg("Finished");
if ($isTrial) {
     $args->msg("[+]","Trial run complete. Human-readable reports at:",
                $trialFile,
                "To load database, re-run with -trial 0");
     close($trialFH);
}

sub pivot_translation {
    my $tab  = "translation";
    my $why  = "Assigns translations to transcripts, used to capture CDS";
    my $fh   = &gzfh( $tab, $why );
    my $cols = &cols_for_tab($tab, $why);
    my %lu;
    while (my $row = &table_row($fh, $cols)) {
        my $rnaid  = $row->{transcript_id} || 0;
        next unless ($rnaid);
        my $protid = $row->{translation_id};
        my $pens   = $row->{stable_id};
        # The CDS coordinates are in EXON space
        # We will need to alter them later when we know how long the exons
        # are and what order they occur in
        $lu{$rnaid} = [ $row->{seq_start},
                        $row->{seq_end},
                        $row->{start_exon_id},
                        $row->{end_exon_id}, 
                        $pens, $protid, ];
        # push @{$dat->{tags}{"Has Protein"}}, $pens;
    }
    close $fh;
    return \%lu;
}

sub pivot_exon {
    my $regLU = &region_lookup();
    my $rnaLU = &exonid_to_transid();
    my $tab   = "exon";
    my $why   = "This table anchors each exon to genomic coordinates";
    my $fh    = &gzfh( $tab, $why );
    my $cols  = &cols_for_tab($tab, $why);
    my (%exonData, %seen);
    my @grab = qw
        (stable_id seq_region_start seq_region_end seq_region_strand);
    while (my $row = &table_row($fh, $cols)) {
        my $exid  = $row->{exon_id};
        my $rDat  = $rnaLU->{$exid};
        my $srid  = $row->{seq_region_id};
        my $chr   = $regLU->{$srid};
        unless ($rDat) {
            $seen{"Failed to map ExonID to RnaID"}++;
            next;
 
        }
        unless ($chr) {
            $seen{ defined $chr ? 
                       "Ignored Accession" : "Irrelevant Accession"}++;
            next;
        }
        my ($ens, $gs, $ge, $str) = map { $row->{$_} } @grab;
        foreach my $info (@{$rDat}) {
            # One exon can be shared by many RNAs
            my ($dbid, $ind) = @{$info};
            my $store = $exonData{$dbid}{$chr}{$str} ||= [];
            $seen{"Anchored"}++;
            if ($store->[$ind - 1]) {
                &huh("Multiple exon locations defined: transcript_id = $dbid Exon $ind");
                next;
            }
            # To the best of my knowledge, Ensembl exon objects are always
            # built explicitly from the reference genome, and always match
            # perfectly without gaps.
            $store->[$ind - 1] = [$ens, $gs, $ge, $exid];
        }
    }
    close $fh;
    &summarize_tally(\%seen, "Exon Mapping Summary:");
    return \%exonData;
}

sub exons {
    my $token = 'RNA_Exons';
    my $pFile = &pivot_file($token);
    my $why   = "Denormalizes the assignments of exons to transcripts";
    unless (-s $pFile) {
        # Need to make the pivot file
        # Doing the pivot in stages to manage memory usage, which
        # was nearly 4 Gb when trying to load all files into RAM

        # This block consumes about 1 Gb to generate the pivot file
        # Loading the pivot file takes 208 Mb
        my $idlu = &transcript_lookup();
        my $elu  = &pivot_exon();
        my $trlu = &pivot_translation();
        my $errFile = "$pFile-Errors.txt";
        open(PIV, ">$pFile") || &polite_exit
            ("Failed to create $token pivot file", $pFile, $!, $why);
        open(PERR, ">$errFile") || &polite_exit
            ("Failed to create $token error file", $errFile, $!);
        print PIV join("\t", qw(rna_id chr strand prot_acc cds
                                hsps errors))."\n";
        foreach my $tid (sort {$a <=> $b} keys %{$elu}) {
            my $tacc = $idlu->{$tid};
            unless ($tacc) {
                my $err = "No accession found for ID";
                &huh("$tid : $err");
                print PERR join("\t", $tid, $err)."\n";
                next;
            }
            my @row = ($tacc, "", "", "", "", "", "");
            my (@data, @errs);
            my @chrs = sort keys %{$elu->{$tid}};
            my $pdat = $trlu->{$tid};
            $row[3] = $pdat->[4] || "" if ($pdat);
            if ($#chrs == 0) {
                my $chr  = $row[1] = $chrs[0];
                my @strs = sort keys %{$elu->{$tid}{$chr}};
                if ($#strs == 0) {
                    my $str = $row[2] = $strs[0];
                    @data   = @{$elu->{$tid}{$chr}{$str}};
                    my (@absent, %offsets, @hsps);
                    my $lastR = 1;
                    for my $i (0..$#data) {
                        unless ($data[$i]) {
                            push @absent, $i+1;
                            next;
                        }
                        my ($ens, $gs, $ge, $exid) = @{$data[$i]};
                        my $lft = $offsets{$exid} = $lastR - 1;
                        my $len = $ge - $gs + 1;
                        push @hsps, join(',', $lft, 
                                         $lastR = $lastR + $len,
                                         $gs - 1, $ge + 1 );
                        # push @{$tags->{"Has Exon"}}, $ens;
                    }
                    if ($#absent == -1) {
                        $row[5] = join(';', @hsps);
                        if ($pdat) {
                            # Update and set the CDS
                            my ($s, $e, $exS, $exE) = @{$pdat};
                            my ($os, $oe) = ($offsets{$exS}, $offsets{$exE});
                            if (defined $os && defined $oe) {
                                $row[4] = join(',', $s + $os, $e + $oe);
                                #$rna->set_cds( $atg, $e + $oe );
                                #if (my $feats = $rdat->{feat}) {
                                #    $rna->clear_features($authority);
                                #    foreach my $fdat (@{$feats}) {
                                #        my ($ps, $pe, $name) = @{$fdat};
                                #        # Adjust protein coordinates to RNA coordinates
                                #        $rna->set_feature
                                #            ( [[$atg + ($ps - 1) * 3, $atg - 1  + ($pe * 3) ]],
                                #              $name, $authority);
                                #    }
                                #}
                            } else {
                                my @miss;
                                push @miss, $exS unless (deinfed $exS);
                                push @miss, $exE unless (defined $exE);
                                push @errs, "Failed to set CDS, no offsets for exon_id = ".
                                     join(' + ', @miss);
                            }
                        }
                    } else {
                        push @errs, "Alignment data missing for exons: ".
                            join(',', @absent);
                    }
                } elsif ($#strs == -1) {
                    push @errs, "No strand defined for alignment";
                } else {
                    push @errs, "Aligned to multiple strands: ".
                        join(',', @strs);
                }
            } elsif ($#chrs == -1) {
                push @errs, "No chromosome alignments";
            } else {
                push @errs, "Aligned to multiple chromosomes: ".
                    join(',', @chrs);
            }
            unless ($#errs == -1) {
                $row[6] = join(';', @errs);
                print PERR join("\t", $tacc, @errs)."\n";
            }
            print PIV join("\t", @row)."\n";
        }
        close PIV;
        close PERR;
        $args->msg("[>]","Exon pivot file created", $pFile );
        if (-s $errFile) {
            $args->msg("[!]", "Errors encountered while pivotting exons",
                       $errFile);
            
        }
    } else {
        $args->msg("[<]","Using pre-generated exon pivot file", $pFile,
                   "Delete file to recalculate");
                   
    }
    my %rv;
    open(PIV, "<$pFile") || &polite_exit
        ("Failed to read $token pivot file", $pFile, $!, $why);
    my $head = <PIV>;
    my $num = 0;
    while (<PIV>) {
        s/[\n\r]+$//;
        my @data = split(/\t/);
        my $id = shift @data;
        if ($rv{$id}) {
            &huh("Multiple exon assignments for $id",
                 join(' + ', @data), join(' + ', @{$rv{$id}}));
        } else {
            $rv{$id} = \@data;
            $num++;
        }
    }
    close PIV;
    $args->msg("[+]","Exon assignments for $num transcripts loaded");
    return \%rv;
}

sub exonid_to_transid {
    my $tab  = "exon_transcript";
    my $why  = "Associates each exon with the relevant RNA, and in the proper order";
    my $fh   = &gzfh( $tab, $why );
    my $cols = &cols_for_tab($tab, $why);
    my %lu;
    while (my $row = &table_row($fh, $cols)) {
        my $dbid  = $row->{transcript_id};
        push @{$lu{$row->{exon_id}}}, [$dbid, $row->{rank}];
    }
    close $fh;
    return \%lu;
}

=head3 Protein featurename counts:

    Seg : 202442
    PF : 157586
    SSF : 129733
    SM : 129593
    PR : 128571
    PS : 116777
    Tmhmm : 57938
    ncoils : 35666
    Sigp : 12719
    TIGR : 4469
    PIRSF : 3977

=cut

sub prot_to_rna_lookup {
    my $luf  = &lookup_file('ProteinToRNA');
    my $what = "Protein to RNA lookup";
    my $why  = "Allow mapping of protein coordinates back to RNA";
    unless (-s $luf) {
        my $plu = &protein_lookup();
        # Reverse the hash
        my $ulp = {};
        while (my ($id, $acc) = each %{$plu}) {
            $ulp->{$acc} = $id;
        }
        undef $plu;
        my $vlu  = &transcript_version_lookup();
        my $exlu = &exons();
        open(LFILE, ">$luf") || &polite_exit
            ("Failed to write $what", $luf, $!, $why);
        foreach my $tacc (sort keys %{$exlu}) {
            my $exdat = $exlu->{$tacc};
            my ($pacc, $cdsTxt) = ($exdat->[2], $exdat->[3]);
            next unless ($pacc);
            my $pid = $ulp->{$pacc};
            unless ($pid) {
                &huh("Failed to find database key for translation", $pacc);
                next;
            }
           
            my $vnum = $vlu->{$tacc};
            unless ($vnum) {
                &huh("Failed to find versioned accession for RNA", $tacc,
                     "Can not assign sequence!");
                next;
            }
            my ($cs, $ce) = split(/\,/, $cdsTxt);
            unless ($cs && $ce) {
                &huh("CDS not defined for mRNA", "tacc = $pacc");
                next;
            }
            print LFILE join("\t", $pid, "$tacc.$vnum", $cs, $ce)."\n";
        }
        close LFILE;
        $args->msg("[>]","$what file created",
                   $luf);
    }
    my %rv;
    open(LFILE, "<$luf") || &polite_exit
        ("Failed to read $what", $luf, $!, $why);
    while (<LFILE>) {
        s/[\n\r]+$//;
        my @row = split(/\t/);
        my $pid = shift @row;
        $rv{$pid} = \@row;
    }
    # die $args->branch(\%rv);
    close LFILE;
    return \%rv;
}


sub _make_feature_pivot {
    my $pFile  = shift;
    my $prlu   = &prot_to_rna_lookup();
    my $tab    = "protein_feature";
    my $why    = "This table anchors each exon to genomic coordinates";
    my $fh     = &gzfh( $tab, $why );
    my $cols   = &cols_for_tab($tab, $why);
    my %ignore = map { $_ => 1 } qw(Seg);
    my (%seen, %byTrans);
    while (my $row = &table_row($fh, $cols)) {
        my $name  = $row->{hit_name};
        next unless ($name && !$ignore{$name});

        my $protid = $row->{translation_id};
        my $rnaDat = $prlu->{$protid};
        unless ($rnaDat) {
            &huh("Failed to find RNA information for protein id $protid");
            next;
        }
        my ($rnaacc, $cs, $ce) = @{$rnaDat};
        # Protein coordinates
        my ($ps, $pe) = ($row->{seq_start}, $row->{seq_end});
        # RNA coordinates
        my ($rs, $re) = ($cs + ($ps - 1) * 3, $cs - 1  + ($pe * 3));
        push @{$byTrans{$rnaacc}}, [$rs, $re, $name];
        my $prfx = $name;
        $prfx =~ s/\d+$//;
        $seen{$prfx}++;
    }
    close $fh;

    open(PIV, ">$pFile") || &polite_exit
        ("Failed to create protein feature pivot file", $pFile, $!);
    foreach my $rnaacc (sort keys %byTrans) {
        my @bits = map { join(' ', @{$_}) } sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[2] cmp $b->[2] } @{$byTrans{$rnaacc}};
        print PIV join("\t", $rnaacc, @bits)."\n";
    }
    close PIV;
    $args->msg("[>]","Feature pivot file created", $pFile);
    &summarize_tally(\%seen, "Observed protein feature prefices:");
}

sub protein_features {
    &location_dbi();
    my $token = 'Protein_Features';
    my $pFile = &pivot_file($token);
    my $why   = "Assigns protein features to transcript space";
    unless (-s $pFile) {
        &_make_feature_pivot( $pFile );
    } else {
        $args->msg("[<]","Using pre-generated feature pivot file", $pFile,
                   "Delete file to recalculate");
    }
    open(PIV, "<$pFile") || &polite_exit
        ("Failed to read $token pivot file", $pFile, $!, $why);
    # my $head = <PIV>;
    my $num = 0;
    my $fnum = 0;
    while (<PIV>) {
        s/[\n\r]+$//;
        my @data = split(/\t/);
        my $accv = shift @data;
        my $rna = BMS::SnpTracker::MapLoc::RNA->new
            ($mapLoc, $accv, $authority);
        $fnum += $#data + 1;
        foreach my $ftxt (@data) {
            if ($ftxt =~ /^(\d+) (\d+) (.+)/) {
                if ($isTrial) {
                    print $trialFH "$accv : [$1..$2] $3\n";
                } else {
                    $rna->set_feature( [[ $1, $2 ]], $3, $authority );
                }
            } else {
                &huh("Failed to interpret RNA feature information", $ftxt);
            }
        }
        $num++;
        if ($limit && $num >= $limit) {
            $args->msg("[LIMIT]","Only loading $limit entries");
            last;
        }
    }
    close PIV;
    $args->msg("[+]","$num RNAs have had $fnum features assigned to them");
}

sub set_exons {
    my ($rna, $exlu) = @_;
    my $acc = $rna->acc();
    $acc =~ s/\.\d+$//;
    my $exd = $exlu->{$acc};
    return unless ($exd);
    my ($chr, $str, $prot, $cdsSE, $hspTxt, $errs) = @{$exd};
    if ($prot) {
        $rna->set_tag_val_ids($textCache{$hasPtag} ||= 
                              $mapLoc->text_to_pkey( $hasPtag ),
                              $mapLoc->text_to_pkey($prot));
    }
    if ($cdsSE) {
        # 377,802
        my ($s, $e) = split(',', $cdsSE);
        $rna->set_cds($s, $e);
    }
    my @hsps;
    foreach my $hspt (split(/\;/, $hspTxt)) {
        # 0,50,85798319,85798369;49,214,85798772,85798937
        push @hsps, [split(/\,/, $hspt)];
    }
    if ($#hsps == -1) {
        return;
    }
    my $aln  = $isTrial ? 
        $mapLoc->get_transient_alignment() : $mapLoc->get_alignment();
    $aln->source($authority);
    $aln->score( 100 );
    $aln->howbad( 0 );
    $aln->strand( $str );
    $aln->seqs($rna->acc(), $chr);
    $aln->hsps(\@hsps);
    $aln->read_or_force_pkey( $chr, 10, 'useBoth', $authority); 
    $rna->add_alignment( $aln );
}

sub region_lookup {
    return $regLU if ($regLU);
    my $tab  = "seq_region";
    my $why  = "This table resolves sequence DB identifiers to names";
    my $fh   = &gzfh( $tab, $why );
    my $cols = &cols_for_tab($tab, $why);
    my %okCrd;
    my $numOk = 0;
    while (my ($ctyp, $crdid) = each %{$gzFiles->{_info}{coordid}}) {
        $okCrd{$crdid} = $ctyp;
        $numOk++;
    }
    &polite_exit
        ("Can not process genomic locations without the internal coord_system_id") unless ($numOk);
    
    &location_dbi();
    my %lu;
    my %seen;
    my $build = $gzFiles->{_info}{token};
    while (my $row = &table_row($fh, $cols)) {
        my $csid = $row->{coord_system_id} || 0;
        my $ctyp = $okCrd{$csid};
        next unless ($ctyp);
        my $frm  = $gzFiles->{_info}{fmt}{$ctyp};
        my $name = $row->{name};
        my ($nice, $type);
        if ($name =~  /_PATCH$/) {
            $type = "Ignored: Patch";
        } elsif ($ctyp eq 'scaffold') {
            if ($name =~ /^[A-Z0-9]+\.\d+$/) {
                $type = "Scaffold";
                $nice = sprintf($frm, $name);
            } else {
                &huh("Unexpected $ctyp name: '$name'");
            }
        } elsif ($ctyp eq 'chromosome') {
            if ($name =~ /^(\d{1,2}|[XYZW]|MT)$/) {
                $type = "Chromosome";
                $nice = sprintf($frm, $name);
            } else {
                &huh("Unexpected $ctyp name: '$name'");
            }
        } elsif ($name =~ /^.{1,2}_(NT_(\d{6}|\d{9}))$/) {
            $type = "Ignored: NT Contig";
            $nice = 0;
        } elsif ($name =~  /^(NC_(\d{6}|\d{9}))$/) {
            $type = "Ignored: NC Contig";
            $nice = 0;
        } elsif ($name =~  /_(MHC)_/) {
            $type = "Ignored: $1";
        } elsif ($name =~  /_CTG\d+(_\d+)?$/) {
            $type = "Ignored: CTG";
        } else {
            $type = "Ignored: Other";
        }
        next unless ($nice);
        $seen{$type}++;
        my $len = $row->{length};
        my $srid = $row->{seq_region_id};
        if ($nice && $lu{$srid}) {
            unless ($nice eq $lu{$srid}) {
                &huh("Multiple accessions for seq_region_id = $srid !!",
                     "$nice, $lu{$srid}");
            }
        } else {
            $lu{$srid} = $nice;
            my $txt = BMS::SnpTracker::MapLoc::Text->new($mapLoc, $nice);
            $txt->tag( 'Build', $build);
            # $txt->tag( 'Description', $bs->desc() );
            $txt->tag( 'Length', $len);
            $txt->tag( 'Name', $name );
            $txt->tag( 'Assembly Type', $ctyp );
            $txt->tag( 'Species', $gzFiles->{_info}{nicetx} );
            if ($isTrial) {
                print $trialFH $txt->to_text();
            } else {
                $txt->write_tags();
                $mapLoc->define_acc_for_loc( $nice, $build, $name );
            }
        }
    }
    close $fh;
    &summarize_tally(\%seen, "Sequence Types:");
    return $regLU = \%lu;
}


sub transcript_lookup {
    return &_generic_lookup( 'transcript' );
}

sub protein_lookup {
    return &_generic_lookup( 'translation' );
}

sub gene_lookup {
    return &_generic_lookup( 'gene' );
}

sub _generic_lookup {
    my $type = shift;
    my $luf  = &lookup_file($type);
    my $why  = "Lookup from ${type}_id to accession";
    my $idkey = $type . "_id";
    unless (-s $luf) {
        my $tab  = $type;
        my $fh   = &gzfh( $tab, $why );
        my $cols = &cols_for_tab($tab, $why);
        my %luh;
        while (my $row = &table_row($fh, $cols)) {
            my $unv  = $row->{stable_id};
            my $tid  = $row->{$idkey};
            next unless ($unv && $tid);
            if ($luh{$unv} && $luh{$unv} != $tid) {
                &huh("$type accession assigned to multiple DB Ids:",
                     "$unv : $luh{$unv} != $tid");
            } else {
                $luh{$unv} = $tid;
            }
        }
        close $fh;
        open(LFILE, ">$luf") || &polite_exit
            ("Failed to write $idkey lookup", $luf, $!, $why);
        my @accs  = sort keys %luh;
        foreach my $acc (@accs) {
            print LFILE "$acc\t$luh{$acc}\n";
        }
        close LFILE;
        $args->msg("[>]","$idkey to accession lookup created",
                   $luf);
    }
    my %rv;
    open(LFILE, "<$luf") || &polite_exit
        ("Failed to read $idkey lookup", $luf, $!, $why);
    while (<LFILE>) {
        s/[\n\r]+$//;
        my ($acc, $id) = split(/\t/);
        $rv{$id} = $acc;
    }
    close LFILE;
    return \%rv;
}

sub transcript_version_lookup {
    my $luf  = &lookup_file('RnaVersion');
    my $what = "RNA version lookup";
    my $why  = "Associate unversioned transcripts with their versioned ID";
    unless (-s $luf) {
        my $tab  = 'transcript';
        my $why  = "Finding transcript version numbers for this release";
        my $fh   = &gzfh( $tab, $why );
        my $cols = &cols_for_tab($tab, $why);

        my $rlu = &transcript_lookup;
        open(LFILE, ">$luf") || &polite_exit
            ("Failed to write $what", $luf, $!, $why);
        while (my $row = &table_row($fh, $cols)) {
            my $unv = $row->{stable_id};
            next unless ($unv);
            my $vnum = $row->{version};
            unless ($vnum) {
                &huh("No version number for $unv");
                next;
            }
            print LFILE "$unv\t$vnum\n";
        }
        close LFILE;
        $args->msg("[>]","$what file created",
                   $luf);
    }
    my %rv;
    open(LFILE, "<$luf") || &polite_exit
        ("Failed to read $what", $luf, $!, $why);
    while (<LFILE>) {
        s/[\n\r]+$//;
        my ($uacc, $vnum) = split(/\t/);
        $rv{$uacc} = $vnum;
    }
    close LFILE;
    return \%rv;
}

sub full_lookup {
    my %rv;
    foreach my $hash (&protein_lookup(),
                      &transcript_lookup(),
                      &gene_lookup()) {
        while (my ($id, $acc) = each %{$hash}) {
            if ($rv{$id}) {
                &huh("Database ID $id is shared by $acc and $rv{$id}");
            } else {
                $rv{$id} = $acc;
            }
        }
    }
    return \%rv;
}

sub transcripts {
    &location_dbi();
    my $glu  = &gene_lookup();
    my $xrlu = &xref_pivot();
    my $exlu = &exons();
    my $tab  = "transcript";
    my $why  = "This table captures basic information about each RNA";
    my $fh   = &gzfh( $tab, $why );
    my $cols = &cols_for_tab($tab, $why);
    my $num  = 0;
    my $start = time;
    my %seen;
    $args->msg("[>>]","Processing RNA objects");
    while (my $row = &table_row($fh, $cols)) {
        my $unv = $row->{stable_id};
        next unless ($unv);

        next if ($only && !$only->{$unv});
        if ($seen{$unv}++) {
            $args->death("Multiple entries for stable_id = $unv",
                         $args->branch($row));
        }
        my $vnum = $row->{version};
        unless ($vnum) {
            &huh("No version number for $unv");
            next;
        }
        my $vers = join('.', $unv, $vnum);
        my $tid  = $row->{transcript_id};
        my $rdat = {
            unv  => $unv,
            vers => $vers,
            vnum => $vnum,
            xref => $row->{display_xref_id},
            dbid => $tid,
            tags => { Error => [], },
        };
        $ens_id_lookup->{$tid} = $rdat;
        my $bt = lc($row->{biotype} || "");
        $rdat->{pseudo} = $bt =~ /pseudo\s?gene/ ? 1 : 0;
        
        &basic_tags($row, $rdat);
        &add_xrefs($rdat, $xrlu);
        if (my $gid = $row->{gene_id}) {
            if (my $gacc = $glu->{$gid}) {
                push @{$rdat->{tags}{"Has Gene"}}, $gacc;
            } else {
                &err("$vers fails to map gene_id = $gid to accession");
        }
        } else {
            &err("No gene_id for $vers");
        }
        my $rna = BMS::SnpTracker::MapLoc::RNA->new
            ($mapLoc, $vers, $authority);
        $rna->read_only() if ($isTrial);
        &set_tags( $rna, $rdat->{tags} );
        &set_exons( $rna, $exlu );
        if ($isTrial) {
            print $trialFH $rna->full_text()."\n";
        } else {
            foreach my $aln ($rna->each_alignment()) {
                $rna->delete_alignment( $aln )
                    unless ($aln->source() eq $authority);
            }
            $rna->write();
        }
        unless (++$num % 50) {
            my $msg = $vers;
            if (my $elaps = time - $start) {
                $msg .= sprintf(" (%.1f/min)", 60 * $num / $elaps);
            }
            $args->msg("[$num]", $msg);
        }
        if ($limit && $num >= $limit) {
            $args->msg("[LIMIT]","Only loading $limit entries");
            last;
        }
    }
    close $fh;
    $args->msg("[+]","Finished load of $num RNAs");
}

sub genes {
    &location_dbi();
    # First scan the RNAs to get list of each RNA associated with a gene
    my $rtab  = "transcript";
    my $rwhy  = "Extracting RNA-to-gene associations";
    my $rfh   = &gzfh( $rtab, $rwhy );
    my $rcols = &cols_for_tab($rtab, $rwhy);
    my %rnaLU;
    while (my $row = &table_row($rfh, $rcols)) {
        if (my $racc = $row->{stable_id}) {
            push @{$rnaLU{ $row->{gene_id} || 0 }}, $racc;
        }
    }
    close $rfh;

    # Then load basic Gene objects from the gene table
    my $xrlu   = &xref_pivot();
    my $tab    = "gene";
    my $why    = "This table is used to deconvolute DB IDs to stable gene identifiers";
    my $fh     = &gzfh( $tab, $why );
    my $cols   = &cols_for_tab($tab, $why);
    my $num    = 0;
    my $start  = time;
    $args->msg("[>>]","Processing Gene objects");
    while (my $row = &table_row($fh, $cols)) {
        my $gid = $row->{gene_id};
        unless ($gid) {
            &err("No gene_id??", $args->branch($row));
            next;
        }
        my $acc = $row->{stable_id};
        my $dat = {
            unv  => $acc,
            vnum => $row->{version},
            xref => $row->{display_xref_id},
            dbid => $gid,
            tags => {
                "Object Type" => [ 'Gene' ],
                "Namespace"   => [ 'Ensembl' ],
            },
        };
        if (my $rnas = $rnaLU{$gid}) {
            push @{$dat->{tags}{"Has RNA"}}, @{$rnas};
        }
        &basic_tags($row, $dat);
        &add_xrefs($dat, $xrlu);
        my $gene = BMS::SnpTracker::MapLoc::Gene->new($mapLoc, $acc);
        $gene->read_only() if ($isTrial);
        &set_tags( $gene, $dat->{tags} );
         if ($isTrial) {
            print $trialFH $gene->to_text()."\n";
        } else {
            $gene->write();
        }
        unless (++$num % 50) {
            my $msg = $acc;
            if (my $elaps = time - $start) {
                $msg .= sprintf(" (%.1f/min)", 60 * $num / $elaps);
            }
            $args->msg("[$num]", $msg);
        }
        if ($limit && $num >= $limit) {
            $args->msg("[LIMIT]","Only loading $limit entries");
            last;
        }
    }
    close $fh;
    $args->msg("[+]","Finished load of $num Genes");
}

sub basic_tags {
    my ($row, $dat) = @_;
    my $tags = $dat->{tags} ||= {};
    $tags->{Taxa}        = [ $gzFiles->{_info}{nicetx} ];
    if (my $unv = $dat->{unv}) {
        $tags->{Unversioned} = [ $unv ];
    }
    $tags->{Deprecated} = $row->{is_current} ? [] : ['Deprecated'];
    while (my ($key, $tag) = each %tagMap) {
        my $val = $row->{$key};
        next unless (defined $val && $val ne "");
        $val =~ s/\s+.+// if ($key =~ /date$/);
        # I find '[Source:RFAM;Acc:RF00402]' distracting:
        $val =~ s/\s*\[[^\]]+\]\s*$// if ($key eq 'description');
        push @{$tags->{$tag}}, $val;
    }
}

sub xref_pivot {
    my $token = 'XREF_Tags';
    my $pFile = &pivot_file($token);
    my $why   = "Collects useful tag/value pairs for objects";
    unless (-s $pFile) {
        &_make_xref_pivot( $pFile );
    } else {
        $args->msg("[<]","Using pre-generated XREF pivot file", $pFile,
                   "Delete file to recalculate");
    }
    my %rv;
    open(PIV, "<$pFile") || &polite_exit
        ("Failed to read $token pivot file", $pFile, $!, $why);
    # my $head = <PIV>;
    my $num = 0;
    while (<PIV>) {
        s/[\n\r]+$//;
        my @data = split(/\t/);
        my $id = shift @data;
        if ($rv{$id}) {
            &huh("Multiple exon assignments for $id",
                 join(' + ', @data), join(' + ', @{$rv{$id}}));
        } else {
            $rv{$id} = \@data;
            $num++;
        }
    }
    close PIV;
    $args->msg("[+]","XREF pairs for $num objects loaded");
    return \%rv;
}

sub add_xrefs {
    my ($dat, $xrlu) = @_;
    if (my $unv = $dat->{unv}) {
        if (my $xd = $xrlu->{$unv}) {
            foreach my $kv (@{$xd}) {
                if ($kv =~ /^([^=]+)=(.+)/) {
                    push @{$dat->{tags}{$1}}, $2;
                }
            }
        }
    }
}

sub _make_xref_pivot {
    my $pFile  = shift;
    my $db2acc = &full_lookup();
    my %xrLU;

    my ($tab, $why, $fh, $cols);
    $tab  = "object_xref";
    $why  = "Links gene and transcript objects to external references (XREFs)";
    $fh   = &gzfh( $tab, $why );
    $cols = &cols_for_tab($tab, $why);
    while (my $row = &table_row($fh, $cols)) {
        my $dbid  = $row->{ensembl_id};
        if (my $acc = $db2acc->{$dbid}) {
            if (my $xid = $row->{xref_id}) {
                push @{$xrLU{$xid}}, $acc;
            }
        } else {
            &huh("No accession for ensembl_id = $dbid");
        }
    }
    close $fh;

    my %tags;

    # Now scan the XREF table
    my $xdb  = &xdb();
    $tab     = "xref";
    $why     = "Provides details on the XREF entries";
    $fh      = &gzfh( $tab, $why );
    $cols    = &cols_for_tab($tab, $why);

    while (my $row = &table_row($fh, $cols)) {
        my $cb = $xdb->{ $row->{external_db_id} };
        next unless ($cb);
        my $xid  = $row->{xref_id};
        foreach my $acc (@{$xrLU{$xid} || []}) {
            my $dat = $tags{$acc} ||= {};
            &{$cb}( $row, $dat );
        }
    }

    # And finally the 'External Synonym' table, which is apparently
    # unofficial synonyms
    $tab  = "external_synonym";
    $why  = "Stores non-official symbols";
    $fh   = &gzfh( $tab, $why );
    $cols = &cols_for_tab($tab, $why);

    while (my $row = &table_row($fh, $cols)) {
        # I wish these were just alternative symbols
        # However, they are also things like SwissProt IDs
        if (my $syn = $row->{synonym}) {
            my $xid  = $row->{xref_id};
            foreach my $acc (@{$xrLU{$xid} || []}) {
                my $dat = $tags{$acc} ||= {};
                push @{$dat->{tags}{"Has Alias"}}, $syn;
            }
        }
    }
    close $fh;

    open(PIV, ">$pFile") || &polite_exit
        ("Failed to create XREF pivot file", $pFile, $!);
    my @accs = sort keys %tags;
    my $num = 0;
    foreach my $acc (@accs) {
        my @row = ($acc);
        foreach my $tag (sort keys %{$tags{$acc}{tags}}) {
            $tag =~ s/[\t=]//g;
            my %u = map { $_ => undef } @{$tags{$acc}{tags}{$tag} || []};
            delete $u{""};
            my @vals = sort keys %u;
            next if ($#vals == -1);
            foreach my $val (@vals) {
                $val =~ s/\P{IsASCII}//g; # Strip out non-ascii
                $val =~ s/[\t\s]+/ /g; # Whitespace runs to single space
                $val =~ s/^ //;
                $val =~ s/ $//;
                next if ($val eq '');
                push @row, "$tag=$val";
                $num++;
            }
        }
        print PIV join("\t", @row)."\n";
    }
    close PIV;
    $args->msg("[>]","XREF pivot file created for $num tag/val pairs in ".
               ($#accs + 1). " accessions", $pFile);
}

sub domain_meta {
    &location_dbi();
    my $xdb  = &xdb('Domain');
    my $tab  = "xref";
    my $why  = "Capturing descriptive text for protein domain accessions";
    my $fh   = &gzfh( $tab, $why );
    my $cols = &cols_for_tab($tab, $why);

    my %info;
    my $num = 0;
    while (my $row = &table_row($fh, $cols)) {
        my $exdb = $row->{external_db_id};
        my $cb   = $xdb->{ $exdb };
        next unless ($cb);
        my $dat  = $info{$exdb} ||= {};
        &{$cb}( $row, $dat );
        $num++;
        if ($limit && $num >= $limit) {
            $args->msg("[LIMIT]","Only loading $limit entries");
            last;
        }
    }
    close $fh;
    $tab  = "interpro";
    $why  = "Associate Interpro with other domain databases";
    $fh   = &gzfh( $tab, $why );
    $cols = &cols_for_tab($tab, $why);
    my %aliases;
    while (my $row = &table_row($fh, $cols)) {
        $aliases{$row->{interpro_ac}}{$row->{id}} = 1;
    }
    close $fh;

    foreach my $db (sort keys %info) {
        foreach my $acc (sort keys %{$info{$db}}) {
            my @objs;
            push @objs, $mapLoc->get_text( $acc );
            foreach my $ali (keys %{$aliases{$acc} || {}}) {
                push @objs, $mapLoc->get_text( $ali );
            }

            while ( my ($tag, $vH) = each %{$info{$db}{$acc}}) {
                my $tid = $textCache{$tag} ||= $mapLoc->text_to_pkey( $tag );
                my @vids = map { $mapLoc->text_to_pkey( $_ ) } keys %{$vH};
                $objs[0]->set_tag_val_ids($tid, @vids);
                if ($tag eq 'Description') {
                    my $sa = $textCache{$sameAs} ||=
                        $mapLoc->text_to_pkey( $sameAs );
                    for my $o (1..$#objs) {
                        $objs[$o]->set_tag_val_ids($tid, @vids);
                        $objs[$o]->set_tag_val_ids
                            ($sa, $objs[0]->pkey());
                    }
                }
            }
            foreach my $obj (@objs) {
                if ($isTrial) {
                    print $trialFH $obj->to_text();
                } else {
                    $obj->write_tags();
                }
            }
        }
    }
    $args->msg("[+]","Processed metadata for $num protein domains");
}

sub xdb {
    my $doDomain = shift;
    my $tab  = "external_db";
    my $why  = "Just some basic information on non-Ensembl data sources";
    my $fh   = &gzfh( $tab, $why );
    my $cols = &cols_for_tab($tab, $why);
    # Assigns aliases from the primary accession
    my $aliAcc = sub {
        my ($row, $dat) = @_;
        my $tags = $dat->{tags} ||= {};
        if (my $pacc = $row->{dbprimary_acc}) {
            push @{$tags->{"Has Alias"}}, $pacc;
        }
    };
    
    # Semi-cautiously assigns aliases from the label
    my $aliLabel = sub {
        my ($row, $dat) = @_;
        my $tags = $dat->{tags} ||= {};
        if (my $dl = $row->{display_label}) {
            # I want reasonably distinct IDs captured here
            # In particular avoiding symbol-like short IDs
            if ($dl =~ /(\.\d+|\-\d+)$/) {
                push @{$tags->{"Has Alias"}}, $dl;
            }
        }
    };

    # Assigns 'Related To' tags from the primary accession
    my $relAcc = sub {
        my ($row, $dat) = @_;
        my $tags = $dat->{tags} ||= {};
        if (my $pacc = $row->{dbprimary_acc}) {
            push @{$tags->{"Related To"}}, $pacc;
        }
    };

    my $capture = {
        Clone_based_vega_transcript    => $aliLabel,
        Clone_based_vega_gene          => $aliLabel,
        Clone_based_ensembl_gene       => $aliLabel,
        Clone_based_ensembl_transcript => $aliLabel,
        RFAM_transcript_name           => $aliLabel,

        RFAM                           => $aliAcc,
        LRG                            => $aliAcc,
        OTTG                           => $aliAcc,
        UCSC                           => $aliAcc,
        UniGene                        => $aliAcc,
        EMBL                           => $aliAcc,

        RefSeq_mRNA                    => $relAcc,
        RefSeq_mRNA_predicted          => $relAcc,
        RefSeq_ncRNA                   => $relAcc,
        RefSeq_ncRNA_predicted         => $relAcc,
        RefSeq_peptide                 => $relAcc,
        RefSeq_peptide_predicted       => $relAcc,
        'Uniprot/SWISSPROT'            => $relAcc,
        protein_id                     => $relAcc,

        Uniprot_gn => sub {
            my ($row, $dat) = @_;
            if (my $dl = $row->{display_label}) {
                my $tags = $dat->{tags} ||= {};
                push @{$tags->{"Symbol"}}, $dl;
            }
        },
        MGI => sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $dl = $row->{display_label}) {
                push @{$tags->{"Official Symbol"}}, $dl;
                push @{$tags->{"Symbol"}}, $dl;
            }
            if (my $pacc = $row->{dbprimary_acc}) {
                push @{$tags->{"MGI"}}, $pacc;
            }
            if (my $desc = $row->{description}) {
                push @{$tags->{Description}}, $desc;
            }
        },
        HGNC => sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $dl = $row->{display_label}) {
                push @{$tags->{"Official Symbol"}}, $dl;
                push @{$tags->{"Symbol"}}, $dl;
            }
            if (my $pacc = $row->{dbprimary_acc}) {
                push @{$tags->{"HGNC"}}, "HGNC:$pacc";
            }
            if (my $desc = $row->{description}) {
                push @{$tags->{Description}}, $desc;
            }

        },
        HGNC_transcript_name => sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $dl = $row->{display_label}) {
                push @{$tags->{"Official Symbol"}}, $dl;
            }
            if (my $desc = $row->{description}) {
                push @{$tags->{Description}}, $desc;
            }
        },
        miRBase_transcript_name => sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            &{$aliLabel}( @_ );
            if (my $dl = $row->{display_label}) {
                push @{$tags->{"miRBase"}}, $dl;
            }
        },
        miRBase =>  sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $pacc = $row->{dbprimary_acc}) {
                push @{$tags->{"miRBase"}}, $pacc;
                push @{$tags->{"Has Alias"}}, $pacc;
            }
            if (my $dl = $row->{display_label}) {
                push @{$tags->{"miRBase"}}, $dl;
                push @{$tags->{"Has Alias"}}, $dl;
            }
        },
        MIM_GENE =>  sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $pacc = $row->{dbprimary_acc}) {
                push @{$tags->{"Has Alias"}}, "MIM:$pacc";
            }
        },
        GO =>  sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $pacc = $row->{dbprimary_acc}) {
                push @{$tags->{"Gene Ontology"}}, $pacc;
            }
        },
        MIM_MORBID =>  sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $pacc = $row->{dbprimary_acc}) {
                my $tag = "MIM:$pacc";
                if (my $desc = $row->{description}) {
                    $tag .= ' : '. $desc;
                }
                push @{$tags->{"MIM Morbid"}}, $tag;
            }
        },   
        Orphanet =>  sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $pacc = $row->{dbprimary_acc}) {
                my $tag = "ORPHA$pacc";
                if (my $desc = $row->{description}) {
                    $tag .= ' : '. $desc;
                }
                push @{$tags->{"Orphanet"}}, $tag;
            }
        },   
        EntrezGene =>  sub {
            my ($row, $dat) = @_;
            my $tags = $dat->{tags} ||= {};
            if (my $pacc = $row->{dbprimary_acc}) {
                push @{$tags->{"Related To"}}, "LOC$pacc";
            }
            if (my $dl = $row->{display_label}) {
                push @{$tags->{"Entrez Symbol"}}, $dl;
                push @{$tags->{"Symbol"}}, $dl;
            }
            if (my $desc = $row->{description}) {
                push @{$tags->{Description}}, $desc;
            }
        },
    };

    # Just ignore these altogether:
    my @doNotWant =  qw
        (OTTT shares_CDS_with_OTTT shares_CDS_and_UTR_with_OTTT HPA PDB
         Ens_Hs_transcript Ens_Hs_translation ENS_LRG_transcript
         Uniprot_genename Uniprot_genename_transcript_name DBASS3 DBASS5 
         Vega_transcript Vega_translation Ens_Hs_gene MEROPS
         ArrayExpress WikiGene ENS_LRG_gene CCDS UniParc goslim_goa);
    $args->msg("[IgnoredDB]","The following databases are being ignored", sort @doNotWant);
    map { $capture->{$_} = "" } @doNotWant;
    $capture->{"Uniprot/SPTREMBL"} = "";

    if ($doDomain) {
        $capture = {};
        foreach my $db (qw(ProSite HAMAP PFAM PRINTS Prints PRODOM ProDom SMART
                           TIGRFAM TigrFam PIRSF Superfamily  CATH-Gene3D
                           PANTHER Panther InterPro Interpro)) {
            $capture->{$db} = sub {
                my ($row, $dat) = @_;
                if (my $pacc = $row->{dbprimary_acc}) {
                    # IPR003116
                    my $targ = $dat->{$pacc} ||= { Database => { $db => 1 } };
                    if (my $dl = $row->{display_label}) {
                        # Raf-like_ras-bd
                        $targ->{Token}{$dl} = 1;
                    }
                    if (my $desc = $row->{description}) {
                        # Raf-like Ras-binding
                        $targ->{Description}{$desc} = 1;
                    }
                }
            };
        }
    }

    my %xdb;
    my $num = 0;
    while (my $row = &table_row($fh, $cols)) {
        my $name = $row->{db_name};
        next unless ($name);
        my $id   = $row->{external_db_id};
        my $cb = $capture->{$name};
        if (defined $cb) {
            $xdb{$id} = $cb;
            $num++;
        } elsif (0) {
            # For testing if there are interesting bits we are not capturing
            $xdb{$id} = sub {
                my $num = ++$ignoredXREF{$name};
                return if ($num > 10);
                my ($row, $dat) = @_;
                my ($acc, $lab, $desc) = map { $_ || "" }
                    ($row->{dbprimary_acc}, $row->{display_label},
                     $row->{description});
                $args->msg("[IgnoredDB]", "$dat->{unv} / $name / $acc [$lab] $desc");
            };
        }
    }
    $args->msg("[+]","Recognized $num XREF database types");
    close $fh;
    return \%xdb;
}

sub set_tags {
    my ($obj, $tagHash) = @_;
    while (my ($tag, $valArr) = each %{$tagHash || {}}) {
        $obj->clear_tag($tag);
        my %u = map { $_ => 1 } @{$valArr};
        my $tid = $textCache{$tag} ||= $mapLoc->text_to_pkey( $tag );
        my @vids;
        foreach my $val (keys %u) {
            $val =~ s/\P{IsASCII}//g; # Strip out non-ascii
            $val =~ s/[\t\s]+/ /g; # Whitespace runs to single space
            push @vids, $mapLoc->text_to_pkey( $val ) unless ($val =~ /^\s*$/);
        }
        $obj->set_tag_val_ids($tid, @vids);
    }
}

sub summarize_tally {
    my ($hash, $title) = @_;
    my @msg;
    push @msg, $title if ($title);
    foreach my $key (sort {$hash->{$b} <=> $hash->{$a} } keys %{$hash}) {
        push @msg, "$key : $hash->{$key}";
    }
    $args->msg("[#]", @msg);
}


sub files {
    my %files;
    $args->msg("[+]", "Reading input from local directory",
               $dir);
    foreach my $file ($args->read_dir( -dir => $dir,
                                       -keep => '.gz', )) {
        my $short = $file;
        $short =~ s/.+\///;
        $short =~ s/\.gz$//;
        if ($short =~ /(.+)\.txt$/) {
            # This should be a table flat file
            $files{$1} = $file;
        } elsif ($short =~ /((.+)_core_(\d+)_(.+))\.sql$/) {
            # SQL used to build database
            my ($fdb, $tx, $dv, $gv) = ($1, lc($2), $3, $4);
            my $niceTx = join(" ", split(/_/, $tx));
            substr($niceTx, 0, 1) = uc(substr($niceTx, 0, 1));
            $files{_info} = {
                fulldb => $fdb,
                taxa   => $tx,
                nicetx => $niceTx,
                dbvers => $dv,
                gvers  => $gv,
            };
            push @{$files{schema}}, $file;
        } elsif ($short =~ /\.(cdna)\.all\.fa$/ ||
                 $short =~ /\.(ncrna)\.fa$/ ) {
            my $type = $1;
            push @{$files{$type."Fasta"}}, $file;
        } else {
            $args->msg("[-]","Ignoring $short");
        }
    }
    return \%files;
}

sub get_build {
    my $info = $gzFiles->{_info} ||= {};
    my $tab  = "coord_system";
    my $why  = "Tracks the different frames of reference that sequences may be scaffolded on. Allows build 'token' to be extracted.";
    my $fh   = &gzfh( $tab, $why );
    my $cols = &cols_for_tab($tab, $why);
    my $tok;
    while (my $row = &table_row($fh, $cols)) {
        next unless ($row->{attrib} && $row->{attrib} eq 'default_version');
        my $name = lc($row->{name} || "");
        &polite_exit("Multiple builds recorded for database",
                     "$tok + $row->{version}")
            if ($tok && $tok ne $row->{version});
        $tok = $row->{version};
        &polite_exit("Multiple coordinate_system_id values",
                     "$info->{coordid} + $row->{coord_system_id}")
            if ($info->{coordid}{$name} && 
                $info->{coordid}{$name} != $row->{coord_system_id});
       $info->{coordid}{$name} = $row->{coord_system_id};
    }
    close $fh;
    $info->{token}  = $tok;
    my @tmp;
    foreach my $ctyp (sort keys %{$info->{coordid}}) {
        push @tmp, "$ctyp template: ". ($info->{fmt}{$ctyp} = sprintf
            ("%s.%s.%%s.%s", $info->{taxa}, $ctyp, $tok || "ERROR"));
    }
    $args->msg("[^]","Database identified as $info->{fulldb}",
               "Species: $info->{taxa}",
               "DB Version: $info->{dbvers}",
               "Genome Version: $info->{gvers}",
               "Build Token: ".($tok || 'NOT FOUND'), @tmp );
    &polite_exit("Failed to identify the build token for '$info->{fulldb}'")
        unless ($tok);
}

sub err {
    $args->msg("[!]", @_);
}

sub huh {
    $args->msg("[?]", @_);
}

sub polite_exit {
    $args->msg("[!!]", @_);
    exit;
}

sub rna_sequence {
    &location_dbi();
    # my $rlu  = &transcript_lookup();
    my $vlu = &transcript_version_lookup();
    foreach my $type ('cdna','ncrna') {
        my $faf = $gzFiles->{$type."Fasta"};
        if (!$faf || $#{$faf} == -1) {
            &polite_exit("Failed to find $type fasta file - it should be of general format", "<species>.<Build>.##.$type.all.fa.gz - for example:", "Homo_sapiens.GRCh37.71.$type.all.fa.gz", "Please download and add to the local directory");
        }
        if ($#{$faf} != 0) {
            &polite_exit("Multiple $type files found:", @{$faf},"Please remove all but the relevant one");
        }
        my $path = $faf->[0];
        my $chk  = join('.', '', $gzFiles->{_info}{token} || '', $gzFiles->{_info}{dbvers} || '', '');
        unless ($path =~ /\Q$chk\E/i) {
            &polite_exit("I found a $type fasta file, but it does not match the apparent database version.", "I was expecting something with '$chk' but found:",$path);
        }
        my $fh = &gzfh($path, "Provides RNA sequence data for $type entries" );
        my $reader = Bio::SeqIO->new( -fh => $fh, -format=> 'fasta');
        my $num  = 0;
        my $start = time;
        my $bps   = 0;
        while (my $bs = $reader->next_seq()) {
            my $uacc = $bs->display_id();
            my $vnum = $vlu->{$uacc};
            unless ($vnum) {
                &huh("Failed to find versioned accession for RNA", $uacc,
                     "Can not assign sequence!");
                next;
            }
            my $vacc = "$uacc.$vnum";
            my $seq = $bs->seq();
            my $len = CORE::length($seq);
            $bps += $len;
            if ($isTrial) {
                print $trialFH "$vacc : $len bp\n";
            } else {
                my $rna = BMS::SnpTracker::MapLoc::RNA->new
                    ($mapLoc, $vacc, $authority);
                $rna->set_seq( $bs->seq() );
            }
            unless (++$num % 50) {
                my $msg = $vacc;
                if (my $elaps = time - $start) {
                    $msg .= sprintf(" (%.1f/min)", 60 * $num / $elaps);
                }
                $msg .= " $bps bp";
                $args->msg("[$num]", $msg);
            }
            if ($limit && $num >= $limit) {
                $args->msg("[LIMIT]","Only loading $limit entries");
                last;
            }
        }
        close $fh;
        $args->msg("[+]","Finished load of $num $type RNAs");
    }
}

sub set_fasta {
    my ($dat, $lines) = @_;
    return unless ($dat);
    my $seq = "";
    foreach my $line (@{$lines}) {
        $line = uc($line || "");
        $line =~ s/[^A-Z]//g;
        $seq .= $line;
    }
    $dat->{sequence} = $seq;
}

sub table_columns {
    my $sch = $gzFiles->{schema};
    if (!$sch || $#{$sch} == -1) {
        &polite_exit("Failed to find schema file - it should be of general format", "<species>_core_##_##.sql.gz - for example:", "homo_sapiens_core_71_37.sql.gz", "Please download and add to the local directory");
    }
    if ($#{$sch} != 0) {
        &polite_exit("Multiple schema files found:", @{$sch},"Please remove all but the relevant one");
    }

    # Kludgy. If the MySQL dump formatting is changed this likely breaks.
    my %tables;
    my $fh = &gzfh( $sch->[0], "Provides column headers for the other files" );
    my $targ;
    while (my $row = <$fh>) {
        if ($row =~ /^\s*CREATE\s+TABLE\s+\`(.+?)\`/) {
            $targ = $tables{$1} ||= [];
            next;
        } elsif ($row =~ /^\s*\`(.+?)\`/) {
            push @{$targ}, $1;
        }
        
    }
    close $fh;
    return \%tables;
}

sub table_row {
    my ($fh, $cols) = @_;
    while (<$fh>) {
        s/[\n\r]+$//;
        my @row = split(/\t/);
        my $rv = {};
        for my $i (0..$#{$cols}) {
            my $val = $row[$i];
            $val = "" if ($row[$i] && $row[$i] eq '\N');
            $rv->{$cols->[$i]} = $val;
        }
        return $rv;
    }
    return undef;
}

sub pivot_file {
    my $token = shift;
    return "$dir/Pivot-$token.tsv";
}

sub lookup_file {
    my $token = shift;
    return "$dir/Lookup-$token.tsv";
}

sub pivot_fh {
    my ($file, $why) = @_;
    my $fh;
    open($fh, "<$file") || &polite_exit
        ("Failed to open pivot file", $file, $!, $why);
    my $short = $file; $short =~ s/^\Q$dir\E/\./;
    $args->msg("[FILE]", "Parsing $short", $why);
    return $fh;
}

sub gzfh {
    my ($path, $why) = @_;
    unless ($path) {
        $args->death("Request for file, but no file provided",
                     "File is interesting because '".($why || 'NO REASON PROVIDED')."'");
    }
    unless (-e $path) {
        if (my $p2 = $gzFiles->{$path}) {
            # The file token was passed rather than the path itself
            $path = $p2;
        } else {
            
        }
    }
    my $fh;
    open($fh, "gunzip -c \"$path\" |") || &polite_exit
        ("Failed to open .gz file", $path, $!, $why);
    my $short = $path; $short =~ s/^\Q$dir\E/\./;
    $args->msg("[FILE]", "Parsing $short", $why);
    return $fh;
}

sub cols_for_tab {
    my ($tab, $why) = @_;
    my $cols = $schema->{$tab};
    return $cols || &polite_exit("Failed to find column names for table $tab",
                                 $why);
}

sub location_dbi {
    return $mapLoc if ($mapLoc);
    #my $makeDb = $args->val(qw(makedatabase makedb rebuild));
    my $mlInt  = $args->val(qw(instance));
    $args->msg("[DB]","Using MapLoc instance '$mlInt'");
    $mapLoc ||= BMS::SnpTracker::MapLoc->new
        ( -build    => $args->val(qw(build)),
          -noenv    => $args->val(qw(noenvironment noenv)),
          -instance => $mlInt,
          # -makedb   => $makeDb,
          );
    #if ($makeDb) {
    #    $args->msg("Built database");
    #    exit;
    #}
    return $mapLoc;
}

sub usage {
    my $txt = <<EOF;
This program will load a MapLoc instance with RNA-to-genome mappings from
Ensembl, as well as some related metadata like symbols and names.
It requires the path to a local directory, eg:
   -dir ~/ensemblFiles
That directory should contain a subset (you do not need all) of the
Ensembl MySQL files. The files can be downloaded from their website:

   ftp://ftp.ensembl.org/pub/current_mysql/
   ftp://ftp.ensembl.org/pub/current_fasta/<species>/cdna/
   ftp://ftp.ensembl.org/pub/current_fasta/<species>/ncrna/

Because Ensembl embeds both the DB version number and the genome build
number in each file path, it is difficult to do this automatically for
you. Presuming that you want human, and that the current version is
'71_1', then you would want to go to folder:

    ftp://ftp.ensembl.org/pub/current_mysql/homo_sapiens_core_71_37/

... and save the following files to your local computer:

    homo_sapiens_core_71_37.sql.gz

You do not need every file from the FTP site! You will need to retrieve:

Directory: <species>_core_<version>:
   <species>_core_<version>.sql.gz
   coord_system.txt.gz
   exon.txt.gz
   exon_transcript.txt.gz
   external_db.txt.gz
   external_synonym.txt.gz
   gene.txt.gz
   interpro.txt.gz
   mapping_session.txt.gz
   object_xref.txt.gz
   protein_feature.txt.gz
   seq_region.txt.gz
   transcript.txt.gz
   translation.txt.gz
   xref.txt.gz

Fasta:
   <species.Build.Version>.cdna.all.fa.gz
   <species.Build.Version>.ncrna.fa.gz

For GRCh37.73 the above files totalled to 180Mb

Once you have the above files stored locally, use the -dir parameter
to provide the path to the containing directory.

EOF

$args->msg("[#]","Usage instructions", split(/[\n\r]+/, $txt));
}
