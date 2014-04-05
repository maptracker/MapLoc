# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::Utilities::SequenceUtilities;
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
use BMS::Utilities;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Scalar::Util qw(weaken);

use vars qw(@ISA);
@ISA      = qw(BMS::Utilities);

our $complementHash = {
    A => 'T',
    C => 'G',
    G => 'C',
    T => 'A',
    U => 'A',
    Y => 'R',
    R => 'Y',
    W => 'W',
    S => 'S',
    K => 'M',
    M => 'K',
    B => 'V',
    V => 'B',
    D => 'H',
    H => 'D',
    X => 'X',
    N => 'N',
    a => 't',
    c => 'g',
    g => 'c',
    t => 'a',
    u => 'a',
    y => 'r',
    r => 'y',
    w => 'w',
    s => 's',
    k => 'm',
    m => 'k',
    b => 'v',
    v => 'b',
    d => 'h',
    h => 'd',
    x => 'x',
    n => 'n',
    ' ' => ' ',
    '-' => '-',
};

our $ambigBases = {
    "-" => [ '-' ],
    "A" => [ 'A' ],
    "C" => [ 'C' ],
    "G" => [ 'G' ],
    "T" => [ 'T' ],
    "M" => [ 'A','C' ],
    "R" => [ 'A','G' ],
    "W" => [ 'A','T' ],
    "S" => [ 'C','G' ],
    "Y" => [ 'C','T' ],
    "K" => [ 'G','T' ],
    "V" => [ 'A','C','G' ],
    "H" => [ 'A','C','T' ],
    "D" => [ 'A','G','T' ],
    "B" => [ 'C','G','T' ],
    "N" => [ 'A','C','G','T' ],
    "X" => [ 'A','C','G','T' ],
};

our $baseChars = uc(join('', sort keys %{$ambigBases}));
$baseChars =~ s/\-//;

our $aminoAcidData = <<EOF;
1:A;3:Ala;name:Alanine;type:common
1:C;3:Cys;name:Cysteine;type:common
1:D;3:Asp;name:Aspartic acid;type:common
1:E;3:Glu;name:Glutamic acid;type:common
1:F;3:Phe;name:Phenylalanine;type:common
1:G;3:Gly;name:Glycine;type:common
1:H;3:His;name:Histidine;type:common
1:I;3:Ile;name:Isoleucine;type:common
1:K;3:Lys;name:Lysine;type:common
1:L;3:Leu;name:Leucine;type:common
1:M;3:Met;name:Methionine;type:common
1:N;3:Asn;name:Asparagine;type:common
1:P;3:Pro;name:Proline;type:common
1:Q;3:Gln;name:Glutamine;type:common
1:R;3:Arg;name:Arginine;type:common
1:S;3:Ser;name:Serine;type:common
1:T;3:Thr;name:Threonine;type:common
1:V;3:Val;name:Valine;type:common
1:W;3:Trp;name:Tryptophan;type:common
1:Y;3:Tyr;name:Tyrosine;type:common

1:U;3:Sec;name:Selenocysteine;type:uncommon
1:O;3:Pyl;name:Pyrrolysine;type:uncommon

1:B;3:Asx;name:Asparagine or aspartic acid;type:ambiguous
1:J;3:Xle;name:Leucine or Isoleucine;type:ambiguous
1:Z;3:Glx;name:Glutamine or glutamic acid;type:ambiguous
1:X;3:Xaa;name:Unspecified or unknown amino acid;alias:Aaa;type:ambiguous

# Abbreviations and symbols in peptide science: a revised guide and commentary
# Jones JH, Journal of Peptide Science 2006; 12: 1-12
# DOI: 10.1002/psc.725

name:Penicillamine;3:Pen;type:variant
name:alpha-Aminosuberic acid;3:Asu;type:variant
name:alpha-Aminoadipic acid;3:Aad;type:variant
name:Azetidine-2-carboxylic acid;3:Aze;type:variant
name:alpha-Aminoisobutyric acid;3:Aib;type:variant
name:gamma-Carboxyglutamic Acid;3:Gla;type:variant
name:4-Hydroxyproline;3:Hyp;type:variant
name:Homophenylalanine;3:Hph;type:variant
name:Cysteic acid;3:Cya;type:variant
name:Homoserine;3:Hse;type:variant
name:t-Leucine;3:Tle;type:variant
name:Norleucine;3:Nle;type:variant
name:Norvaline;3:Nva;type:variant
name:Ornithine;3:Orn;type:variant
name:Pyroglutamic Acid;3:Pyr;alias:pGlu;type:variant
name:Phenylglycine;3:Phg;type:variant
name:Pipecolic acid;3:Pip;type:variant
name:Statine;3:Sta;type:variant
name:pyroglutamic acid;3:Glp;type:variant
name:Hippuric acid;3:Hip;type:variant
name:Sarcosine;3:Sar;type:variant
name:beta-Thienylalanine;3:Thi;type:variant
name:Thiaproline;3:Thz;type:variant
name:5-oxo-L-proline;3:Pyr;type:variant
name:L-alpha-Aminobutyric acid;3:Abu;type:variant
name:Norvaline;3:Nva;type:variant
name:Anthraniloyl;3:Abz;type:variant
name:Nitrobenzylamide;3:Nba;type:variant
name:Fluorescein-thiocarbamoyl;3:FTC;type:variant
name:beta-cyclohexylalanine;3:Cha;type:variant
name:3-(2,4-dinitrophenyl)-L-2,3-diaminopropionyl;3:Dpa;type:variant
name:1,2,3,4-Tetrahydroisoquinoline-3-carboxylic acid;3:Tic;type:variant
name:Agmatine;3:Agm;type:variant
name:Citrulline;3:Cit;type:variant
name:Dehydroalanine;3:Dha;type:variant
name:5H-dibenzo[a,d]cycloheptene-5-glycine;3:Bhg;type:variant
name:Morpholinoureidyl;3:Mu;type:variant
name:Dimethoxycinnamoyl;3:DMC;type:variant
name:2-naphthylalanine;3:NaI;type:variant
name:delta-Hydroxylysine;3:Hyl;type:variant
name:Lanthionine;3:Lan;type:variant
name:3,3'-diphenylalanine;3:Dip;type:variant
name:N-Methylalanine;3:MeAla;type:variant
name:allo-Threonine;3:alphaThr;type:variant
name:beta-aminoadipic acid;3:betaAad;type:variant
name:P-Alanine;3:PAla;type:variant
name:Isovaline;3:Iva;type:variant

token:Ac;name:Acetyl;type:group
token:Acm;name:Acetamidomethyl;type:group
token:Adoc;name:1-Adamantyloxycarbonyl;type:group
token:Alloc;name:Allyloxycarbonyl;type:group
token:Boc;name:t-Butoxycarbonyl;type:group
token:Bom;name:pi-Benzyloxymethyl;type:group
token:Bpoc;name:2-(4-Biphenylyl)isopropoxycarbonyl;type:group
token:Btm;name:Benzylthiomethyl;type:group
token:Bum;name:Pi-t-Butoxymethyl;type:group
token:Bui;name:i-Butyl;type:group
token:Bun;name:n-Butyl;type:group
token:But;name:t-Butyl;type:group
token:Bz;name:Benzoyl;type:group
token:Bzl;name:Benzyl;type:group
token:Cha;name:Cyclohexylammonium salt;type:group
token:Clt;name:2-Chlorotrityl;type:group
token:Dcha;name:Dicyclohexylammonium salt;type:group
token:Dde;name:1-(4,4-Dimethyl-2,6-dioxocyclohex-1-ylidene)ethyl;type:group
token:Ddz;name:2-(3,5-Dimethoxyphenyl)-isopropoxycarbonyl;type:group
token:Dnp;name:2,4-Dinitrophenyl;type:group
token:Dpm;name:Diphenylmethyl;type:group
token:Bzh;name:Benzhydryl;type:group
token:Dpp;name:Diphenylphosphinyl;type:group
token:Et;name:Ethyl;type:group
token:Fmoc;name:9-Fluorenylmethoxycarbonyl;type:group
# Need to decide how to cope with N-Formyl-Met variants
token:For;name:Formyl;type:group;alias:Formyl;alias:N-Formyl
token:fMet;name:N-Formylmethionine;type:group;alias:N-Formyl-Met
token:Mbh;name:4,4'-Dimethoxydiphenylmethyl,4,4'-Dimethoxybenzhydryl;type:group
token:Mbs;name:4-Methoxybenzenesulphonyl;type:group
token:Me;name:Methyl;type:group
token:Mob;name:4-Methoxybenzyl;type:group
token:Mtr;name:2,3,6-Trimethyl,4-methoxybenzenesulphonyl;type:group
token:Nps;name:2-Nitrophenylsulphenyl;type:group
token:OAll;name:Allyl ester;type:group
token:OBt;name:1-Benzotriazolyl ester;type:group
token:OcHx;name:Cyclohexyl ester;type:group
token:ONp;name:4-Nitrophenyl ester;type:group
token:OPcp;name:Pentachlorophenyl ester;type:group
token:OPfp;name:Pentafluorophenyl ester;type:group
token:OSu;name:Succinimido ester;type:group
token:OTce;name:2,2,2-Trichloroethyl ester;type:group
token:OTcp;name:2,4,5-Trichlorophenyl ester;type:group
token:Tmob;name:2,4,5-Trimethoxybenzyl;type:group
token:Mtt;name:4-Methyltrityl;type:group
token:Pac;name:Phenacyl;conflict:PhCH2CO;type:group
token:Ph;name:Phenyl;type:group
token:Pht;name:Phthaloyl;type:group
token:Scm;name:Methoxycarbonylsulphenyl;type:group
token:Pmc;name:2,2,5,7,8-Pentamethylchroman-6-sulphonyl;type:group
token:Pri;name:i-Propyl;type:group
token:Prn;name:n-Propyl;type:group
token:Tfa;name:Trifluoroacetyl;type:group
token:Tos;name:4-Toluenesulphonyl;type:group
token:Ts;name:4-Toluenesulphonyl;type:group
token:Troc;name:2,2,2-Trichloroethoxycarbonyl;type:group
token:Trt;name:Trityl, triphenylmethyl;type:group
token:Xan;name:9-Xanthydryl;type:group
token:Z;name:Benzyloxycarbonyl;type:group
token:NCA;name:N-Carboxyanhydride;type:group
token:PTH;name:Phenylthiohydantoin;type:group
token:UNCA;name:Urethane N-carboxyanhydride;type:group

# Assorted
token:NH2;name:Amine;type:group

EOF

# ' Make emacs happy

our ($ambigKey, $ambigMatchHash, $cachedTranslations, $cachedSafeDna);

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Returns the reverse complement of a DNA sequence
sub revcom {
    my $self  = shift;
    my ($input) = @_;
    if ($input) {
        if (my $r = ref($input)) {
            if ($r eq 'ARRAY') {
                my @rc = reverse $self->complement( @_ );
                return wantarray ? @rc : \@rc;
            } else {
                $self->err("Do not know how to revcom() query '$input'");
                return undef;
            }
        } else {
            # String
            my $rc = $self->revcom( [ split('', $input) ] );
            return wantarray ? @{$rc} : join('', @{$rc});
        }
    }
    return wantarray ? () : "";
}

sub complement {
    my $self  = shift;
    my ($input, $subInd) = @_;
    if ($input) {
        if (my $r = ref($input)) {
            if ($r eq 'ARRAY') {
                my $comp;
                if (defined $subInd) {
                    if ($subInd =~ /^\d+$/) {
                        $comp = [];
                        foreach my $arr (@{$input}) {
                            my @locArr = @{$arr};
                            $locArr[$subInd] = $complementHash->
                            {$locArr[$subInd]||""} || $locArr[$subInd];
                            push @{$comp}, \@locArr;
                        }
                    } else {
                        $self->err("Can not use '$subInd' as an array index");
                        return wantarray ? () : [];
                    }
                } else {
                    $comp = [ map { $complementHash->{$_} || $_ } @{$input}];
                }
                return wantarray ? @{$comp} : $comp;
            } else {
                $self->err("Do not know how to complement() query '$input'");
                return undef;
            }
        } else {
            # String
            my $comp = $self->complement( [ split('', $input) ] );
            return wantarray ? @{$comp} : join('', @{$comp});
        }
    }
    return wantarray ? () : "";
}

sub guess_moltype {
    my $self = shift;
    my ($seq) = @_;
    my $molt  = "";
    if (my $r = ref($seq)) {
        if ($r =~ /Bio/) {
            $seq = $seq->seq();
        } else {
            $self->err("Not sure how to handle passed sequence '$seq'");
            return $molt;
        }
    }
    $seq = uc($seq);
    $seq =~ s/[^A-Z]+//g;
    my $len = CORE::length($seq);
    return $molt unless ($len);
    $seq =~ s/[ACTGUN]+//g;
    return (CORE::length($seq)/$len > 0.5) ? 'PROTEIN' : 'DNA';
}
    
sub ambiguous_code {
    my $self = shift;
    my ($input) = @_;
    my @chars   = ref($input) ? @{$input} : split('', $input);
    my %uniq;
    map { $uniq{$_} = 1 } map { @{$ambigBases->{uc($_ || "")} || []} } @chars;
    my $code = join('', sort keys %uniq);
    unless ($ambigKey) {
        $ambigKey = {};
        while (my ($ambig, $arr) = each %{$ambigBases}) {
            if ($arr && $#{$arr} != -1) {
                $ambigKey->{ join('', sort @{$arr}) } = 
                    $ambig eq 'X' ? 'N' : $ambig;
            }
        }
    }
    return $ambigKey->{$code || ""} || "N";
}

sub ambiguous_possibilities {
    my $self = shift;
    my $char = uc(shift || "");
    my @rv;
    if (exists $ambigBases->{$char}) {
        @rv = @{$ambigBases->{$char} || [ $char ]};
    } else {
        @rv = ($char);
    }
    return wantarray ? @rv : \@rv;
}

sub ambiguous_translation {
    my $self = shift;
    my $req  = uc(shift || "");
    my @nucs;
    foreach my $bp (split('', $req)) {
        next if ($bp eq '-');
        if (my $ambig = $ambigBases->{$bp}) {
            # This is a recognized base
            push @nucs, $ambig;
        }
    }
    # Now go codon by codon
    my $aas = [["", 1]];
    my @n = ['A','C','T','G'];
    for (my $i = 0; $i <= $#nucs; $i += 3) {
        # Assemble all possible codons for this position
        my $codons = [""];
        foreach my $ambig ($nucs[$i], $nucs[$i+1] || \@n, $nucs[$i+2] || \@n) {
            my @newCodon;
            foreach my $cdn (@{$codons}) {
                foreach my $bp (@{$ambig}) {
                    push @newCodon, $cdn . $bp;
                }
            }
            $codons = \@newCodon;
        }
        my %aaH;
        foreach my $cdn (@{$codons}) {
            $aaH{ $self->cached_translation( $cdn ) }++;
        }
        my @newAA;
        foreach my $pd (@{$aas}) {
            my ($oldprt, $cnt) = @{$pd};
            while (my ($aa, $num) = each %aaH) {
                push @newAA, [ $oldprt.$aa, $cnt * $num ];
            }
        }
        $aas = \@newAA;
    }
    return $aas;
}

sub is_ambiguous_match {
    my $self = shift;
    my ($char1, $char2) = @_;
    return 0 unless ($char1 && $char2);
    unless ($ambigMatchHash) {
        my @allChar = keys (%{$ambigBases});
        # How much ambiguity does each base have (1 lowest, 4 highest)
        my %cSize   = map { $_ => $#{$ambigBases->{$_}} + 1 } @allChar;
        foreach my $c1 (@allChar) {
            my %cMem = map { $_ => 1 } @{$ambigBases->{$c1}};
            my $sz1  = $cSize{$c1};
            foreach my $c2 (@allChar) {
                my $matches = 0; 
                map { $matches += $cMem{$_} || 0 } @{$ambigBases->{$c2}};
                $ambigMatchHash->{$c1}{$c2} = $matches / ($sz1 * $cSize{$c2});
            }
        }
    }
    return 0 unless exists ($ambigMatchHash->{uc($char1)});
    return 0 unless exists ($ambigMatchHash->{uc($char1)}{uc($char2)});
    return $ambigMatchHash->{uc($char1)}{uc($char2)};
}

sub safe_dna {
    my $self = shift;
    my $dna  = shift || "";
    my $xtra = uc(shift || "");
    unless (uc($dna) =~ /^[$baseChars$xtra]+$/) {
        $dna =~ s/[^$baseChars$xtra]/N/ig;
    }
    return $dna;
}

sub cached_safe_dna {
    my $self = shift;
    my $dna  = shift || "";
    unless (defined $cachedSafeDna->{$dna}) {
        $cachedSafeDna->{$dna} = $self->safe_dna($dna, shift);
    }
    return $cachedSafeDna->{$dna};
}

sub translate_string {
    my $self = shift;
    my ($dna, $table) = @_;
    return "" unless ($dna);
    $table ||= 1;
    $self->death("If a second argument (the codon translation table id) is provided to translate_string(), it must be an integer, not '$table'") unless ($table =~ /^\d+$/);
    my $bs = Bio::PrimarySeq->new
        ( -seq        => $self->safe_dna($dna),
          -alphabet   => 'DNA',
          -display_id => "Seq_$dna", );
    my $prot = $bs->translate( -codontable_id => $table );
    return $prot->seq();
}

sub cached_translation {
    # Primarily useful in whacking out translations for ambiguous codons
    my $self = shift;
    my ($dna, $table) = @_;
    return "" unless ($dna);
    $table ||= 1;
    unless ($cachedTranslations->{$table}{$dna}) {
        my $prot = $self->translate_string($dna, $table);
        $cachedTranslations->{$table}{$dna} = uc($prot);
    }
    return $cachedTranslations->{$table}{$dna};
}

sub codon_table {
    my $self = shift;
    my $table = shift || 1;
    my $rv    = $self->{CODON_TABLE}{$table};
    unless ($rv) {
        $rv = $self->{CODON_TABLE}{$table} = {
            BYAA => { '-' => ['---'] },
        };
        my @nucs = qw(A C G T);
        for my $i (0..3) { for my $j (0..3) { for my $k (0..3) {
            my $cod = $nucs[$i] . $nucs[$j] . $nucs[$k];
            my $aa  = $self->cached_translation( $cod, $table );
            $rv->{CODON}{$cod} = $aa;
            push @{$rv->{BYAA}{$aa}}, $cod;
            push @{$rv->{BYAA}{'X'}}, $cod;
        } } }
        while (my ($aa, $clist) = each %{$rv->{BYAA}}) {
            my $cod = "";
            for my $i (0..2) {
                my %nucs = map { substr($_, $i, 1) => 1 } @{$clist};
                my $ambig = $self->ambiguous_code( [ keys %nucs ] );
                $cod .= $ambig;
            }
            $rv->{AMBIGAA}{$aa} = $cod;
        }
    }
    return wantarray ? %{$rv->{CODON}} : $rv;
}

sub aa_to_codon {
    my $self = shift;
    my ($aa, $table) = @_;
    return wantarray ? () : "" unless ($aa);
    return $self->protein_to_nucleotide( $aa, $table )
        unless (CORE::length($aa) == 1);
    my $ct = $self->codon_table($table);
    $aa = uc($aa);
    return wantarray ?
        @{$ct->{BYAA}{$aa} || []} : $ct->{AMBIGAA}{$aa} || 'NNN';
}

sub protein_to_nucleotide {
    my $self = shift;
    my ($prot, $table) = @_;
    my $nuc  = "";
    foreach my $aa (split('', $prot || "")) {
        $nuc .= $self->aa_to_codon( $aa, $table );
    }
    return $nuc;
}

sub protein_align_to_nucleotide {
    my $self = shift;
    my ($string, $nuc) = @_;
    my @errs;
    my $rv = "";
    if (!$string || !$nuc) {
        push @errs, "No alignment string provided" unless ($string);
        push @errs, "No nucleotide sequence provided" unless ($nuc);
    } else {
        my $check = $nuc;
        $check =~ s/[A-Z]//ig;
        if ($check) {
            $nuc =~ s/[^A-Z]//ig;
            push @errs, "Non-alphabetic characters in nucleotide sequence removed: '$check'";
        }
        my $nLen = CORE::length($nuc);
        my @codons;
        for (my $i = 0; $i < $nLen; $i += 3) {
            push @codons, substr($nuc, $i, 3);
        }
        my $llen = CORE::length($codons[-1]);
        unless ($llen == 3) {
            $codons[-1] .= 'N' x (3 - $llen);
            push @errs, "$llen bp terminal codon N-padded to '$codons[-1]'";
        }
        my @aas = split('', uc($string));
        my $tooFewNucs;
        while (my $aa = shift @aas) {
            if ($aa =~ /[A-Z\*]/) {
                if ($#codons == -1) {
                    $rv .= 'NNN';
                    unless ($tooFewNucs) {
                        my $num = 1;
                        map { $num++ if /[A-Z\*]/ } @aas;
                        $tooFewNucs ||= "Codons exhausted with $num AAs remaining";
                    }
                } else {
                    $rv .= shift @codons;
                }
            } else {
                $rv .= $aa x 3;
            }
        }
        unless ($#codons == -1) {
            push @errs, "Mapping leaves ".scalar(@codons)." codons left over";
        }
    }
    return wantarray ? ($rv, \@errs) : $rv;
}

sub gapped_codon_translation {
    my $self = shift;
    my $string = uc(shift || "");
    my $trans  = "";
    my %posErr;
    return wantarray ? ($trans, \%posErr) : $trans unless ($string);

    my ($lft, $rgt) = ("");
    if ($string =~ /^(.*?)([A-Z].+)/) {
        ($lft, $rgt) = ($1, $2);
    }
    if (my $padLft = CORE::length($lft) % 3) {
        # We need to n-pad a codon on the left side
        map { chop($lft) } (1..$padLft);
        $string = join('', $lft, "n" x $padLft, $rgt);
    }
    ($lft, $rgt) = ("","");
    if ($string =~ /(.+[A-Z])(.*?)$/) {
        ($lft, $rgt) = ($1, $2);
    }
    if (my $padRgt = CORE::length($lft) % 3) {
        # We need to n-pad a codon on the right side
        $padRgt = 3 - $padRgt;
        substr($rgt, 0, $padRgt) = "n" x $padRgt;
        $string = join('', $lft, $rgt);
    }

    # $string =~ s/ /-/g;

    # Length of string:
    my $len = CORE::length($string);
    for (my $i = 0; $i < $len; $i += 3) {
        my $codon = substr($string, $i, 3);
        if ($codon =~ /^[A-Zn]{1,3}$/) {
            # Pure letter codon (no gaps / weird characters)
            my $aa = $self->cached_translation( $codon );
            my $aL = CORE::length($aa);
            if ($aL == 1) {
                $trans .= $aa;
            } else {
                $trans .= "X";
                push @{$posErr{"CodonLength$aL"}}, ($i/3) + 1
                    unless ($aL == 0);
            }
        } elsif ($codon =~ /[A-Zn]/) {
            # Oopsies. Codon that is part letter and part other
            my @bits;
            push @bits, substr($string, $i - 6, 6) if ($i >= 6);
            push @bits, "[$codon]";
            push @bits, substr($string, $i+1, 6) if ($i + 6 < $len);
            push @{$posErr{GapInsideCodon}},
            (($i/3) + 1)." ".join('', @bits);
            $codon =~ s/[^A-Zn]//g;
            $trans .= $codon ? $self->cached_translation( $codon ) : 'X';
        } else {
            # Pure non-letter codon
            my $char = substr($codon, 0, 1);

            # WE NEED TO CHECK THIS FOR CONSISTENCY

            $trans .= $char;
        }
    }
    return wantarray ? ($trans, \%posErr) : $trans;
}

sub best_orf {
    my $self = shift;
    my $args  = $self->parseparams( $#_ == 0 ? (-seq => $_[0]) : @_ );
    my @bss   = $self->standardize_seqs( $args->{SEQ} );
    my @notes;
    my $rv    = { note => \@notes };
    if ($#bss == -1) {
        return $rv;
    } elsif ($#bss != 0) {
        push @notes,
        "Multiple sequences provided to best_orf(), using only first";
        @bss = ($bss[0]);
    }
    my $bs  = $rv->{src} = $bss[0];
    my $par = $bs->{SU_PARENT};
    if ($par && $par->can('top_SeqFeatures')) {
        my @cdss;
        map { push @cdss, $_ if ($_->primary_tag eq 'CDS')
              } $par->top_SeqFeatures();
        my %uniq;
        foreach my $cds (@cdss) {
            my $loc  = $cds->location;
            my @locs = $loc->isa('Bio::Location::Split') ? 
                $loc->sub_Location : ($loc);
            my $key  = join("\t", $loc->strand, map { 
                $loc->start, $loc->end} @locs);
            $uniq{$key} ||= $cds;
        }
        my @locs = keys %uniq;
        if ($#locs == 0) {
            # There is a unique CDS feature on this sequence
            my $feat   = $uniq{$locs[0]};
            my @se     = split("\t", $locs[0]);
            $rv->{str} = shift @se;
            $rv->{orf} = \@se;
            $rv->{seq} = uc($feat->seq()->seq());
        } elsif ($#locs > 0) {
            push @notes, "Multiple differing CDS entries found, ignoring all";
        }
    }
    my $len = $bs->length();
    unless ($rv->{orf}) {
        my @frames;
        if (my $f = $args->{FRAME}) {
            if ($f =~ /^[\+\-]?[123]$/) {
                push @frames, $f;
            } else {
                push @notes, "Failed to understand frame request '$f'";
            }
        } elsif (my $str = $args->{STRAND}) {
            if ($str =~ /-1|rev/i) {
                @frames = (-3..-1);
            } elsif ($str =~ /\+1|fwd/i) {
                @frames = (1..3);
            } else {
                push @notes, "Failed to understand strand request '$str'";
            }
        }
        if ($#frames == -1) {
            @frames = (-3..-1,1..3);
        }
        my $forceATG = $args->{FORCEATG};
        my @results;
        foreach my $frame (@frames) {
            my ($seq, $str) = ($bs, 1);
            if ($frame < 0) {
                $frame = abs($frame);
                $str   = -1;
                $seq   = $seq->revcom();
            }
            my $pbs  = $seq->translate(undef, undef, $frame - 1);
            my $prot = uc($pbs->seq);
            my @pept = split(/\*/, $prot);
            map { s/^[^M]+// } @pept if ($forceATG);
            my ($best) = sort { CORE::length($b) <=> CORE::length($a) ||
                                    $a cmp $b } @pept;
            my $ind = index($prot, $best);
            unless (defined $ind) {
                push @notes, "Error finding substring '$best' in '$prot'";
                next;
            }
            my $olen = CORE::length($best) * 3;
            my $s    = $frame + $ind * 3;
            my $e    = $s - 1 + $olen;
            my $oseq = substr($seq->seq(), $s - 1, $olen);
            if ($str < 0) {
                ($s, $e) = ($len - $e + 1, $len - $s + 1);
            }
            push @results, [$olen, $frame, $s, $e, $str, $oseq ];
        }
        my ($br) = sort { $b->[0] <=> $a->[0] || $b->[1] <=> $a->[1]} @results;
        my ($olen, $frame, $s, $e, $str, $oseq) = @{$br};
        $rv->{str} = $str;
        $rv->{orf} = [ $s, $e ];
        $rv->{seq} = $oseq;
    }
    if (my $rng = $rv->{orf}) {
        my @srt = sort { $a <=> $b } @{$rng};
        my ($min, $max, $str) = ($srt[0] || 0, $srt[-1] || 0, $rv->{str} || 0);
        my ($lft, $rgt) = $str < 0 ? ('utr3','utr5') : ('utr5','utr3');
        $rv->{$lft} = [1, $min - 1]  if ($min && $min > 1);
        $rv->{$rgt} = [$max+1, $len] if ($max && $max < $len);
    }
    return $rv;
}

=head2 standardize_seqs

 Title   : standardize_seqs
 Usage   : my @seqs = $fcw->standardize_seqs( $seq1, $seq2, ... )
 Function: Standardizes sequenecs to a list of BioPerl sequence objects
 Returns : An array of Bio::PrimarySeq objects
 Args    : One or more sequence requests

This method is designed to standardize user-supplied sequences. It
will consider each member of the input, recognizing:

 Any BioPerl compliant sequence object
 Raw sequence text
 Raw fasta file text (not tested yet...)

=cut

sub standardize_seqs {
    my $self = shift;
    my @reqlist = @_;
    my @seqs;
    while ($#reqlist != -1) {
        my $seq = shift @reqlist;
        next unless ($seq);
        my $ref = ref($seq) || "";
        if ($ref eq 'ARRAY') {
            # Presume a list of sequences
            push @reqlist, @{$seq};
            next;
        }
        if ($ref eq 'Bio::PrimarySeq') {
            # Already a primary sequence, keep as is
            # Added 5 Apr 2011 - may cause trouble somewhere?
        } elsif ($ref) {
            # Get the PrimarySeq. This has been made more difficult by BioPerl
            my $ps;
            if ($seq->can('primary_seq')) {
                $ps = $seq->primary_seq();
                # Bio::Seq::Quality seems to have primary_seq() as a method
                # but it does not get populated. GRRR.
            }
            if (!$ps && $seq->can('seq')) {
                my $d = $seq->can('desc') ? $seq->desc() : "";
                my $n = $seq->can('display_id') ? $seq->display_id() : "";
                $ps = Bio::PrimarySeq->new
                    ( -seq         => $seq->seq(),
                      -alphabet    => $seq->alphabet(),
                      -desc        => $d || "Explicit user sequence",
                      -display_id  => $n || "SEQ". ++$self->{COUNTER});
            }
            if ($ps) {
                $ps->{SU_PARENT} = $seq;
                $seq = $ps;
            } else {
                $self->err("Failed to recover PrimarySeq", $seq);
            }
        } elsif ($seq =~ /^>/) {
            # Fasta format
            my $obj;
            foreach my $line (split(/[\n\r]+/, $seq)) {
                if ($line =~ /^>(\S+)(.*?)$/) {
                    if ($obj && $obj->{CAT_SEQ}) {
                        push @reqlist, $obj;
                        $obj->seq( $obj->{CAT_SEQ} );
                        delete $obj->{CAT_SEQ};
                    }
                    my ($name, $desc) = ($1, defined $2 ? $2 : "");
                    $desc =~ s/^\s+//; $desc =~ s/\s+$//;
                    $obj = Bio::PrimarySeq->new
                        ( -desc       => $desc,
                          -display_id  => $name );
                    $obj->{CAT_SEQ} = "";
                } else {
                    $obj->{CAT_SEQ} .= $line;
                }
            }
            if ($obj && $obj->{CAT_SEQ}) {
                push @reqlist, $obj;
                $obj->seq( $obj->{CAT_SEQ} );
                delete $obj->{CAT_SEQ};
            }
            next;
        } elsif (-e $seq && (length($seq) > 8 || $seq =~ /.+\..{2,3}/)) {
            # File of sequences
            # There were problems when the sequence was 'X' and a file called
            # X existed in the run path. So a regexp check was added to see
            # that it had something that looked like a file name
            # this will be imperfect, of course
            my $reader;
            eval {
                $reader = Bio::SeqIO->new
                    ( -file   => $seq );
            };
            $self->death("Failed to read sequence file '$seq'")
                unless ($reader);
            my $count = 0;
            my $limit = 50;
            while (my $bs = $reader->next_seq()) {
                # Just grab one sequence for the sake of it
                push @reqlist, $bs;
                if (++$count >= $limit) {
                    $self->err("Halting read of sequences from $seq at $limit");
                    last;
                }
            }
            $reader = undef;
            next;
        } else {
            $seq = Bio::PrimarySeq->new
                ( -seq         => $seq,
                  -desc        => "Literal sequence provided by user",
                  -display_id  => "SEQ". ++$self->{COUNTER});
        }

        my $update = 0;
        my $sd = $seq->seq();
        my $newDesc;
        unless ($sd) {
            $sd = "X";
            $newDesc = "EMPTY SEQUENCE!! ". ($seq->desc() || "");
            $update++;
        }

        # This should be handled automatically by character remapping
        # when auto_remap() is turned on!
        # if ($sd =~ /\*/) {
        #    # Replace stop "amino acids" with Xs
        #    $sd =~ s/\*/x/g;
        #    $update++;
        # }
        if ($sd =~ /[ \-]/) {
            # Remove gap tokens;
            $sd =~ s/[ \-]+//g;
            $update++;
        }
        if ($update) {
            # We don't want to modify the original sequence; it may be
            # used by another program
            $seq = Bio::PrimarySeq->new
                ( -seq         => $sd,
                  -desc        => $newDesc || $seq->desc(),
                  -display_id  => $seq->display_id() );
        }
        push @seqs, $seq;
    }
    return @seqs;
}

sub standardize_peptide {
    my $self = shift;
    my $args  = $self->parseparams( $#_ == 0 ? (-seq => $_[0]) : @_ );
    my $pep   = $args->{SEQ} || $args->{PEPTIDE} ||
        $args->{PROT} || $args->{PROTEIN};
    return wantarray ? ( "", ["No input provided"]) : "" unless ($pep);
    unless ($pep =~ /\-/) {
        # No dashes - the sequence does not appear to be in 
        # Cys-Tyr-Phe-Gln
        # Nomenclature
        my $err;
        my $simp = $self->_clean_aa_input( uc($pep) );
        $simp =~ s/\s+//g;
        my $junk = $simp;
        $junk =~ s/[A-Z]+//g;
        if ($junk) {
            $simp =~ s/[^A-Z]+//g;
            $err = [ "Unrecognized characters '$junk'" ];
        }
        return wantarray ? ($simp, $err) : $simp;
    }
    # Ok, we have dashes. We need to pick apart the sequence
    my $info = $self->aa_info();
    my $hack = $pep;
    $hack =~ s/\s+\-\s+/-/g; # remove spaces around dashes
    my @errs;
    # First we need to map over the aliases,
    # since they might themselves have dashes
    while (my ($ali, $prefer) = each %{$info->{ALIAS}}) {
        my $iloop = 0;
        while ($hack =~ /([\<\-])(\Q$ali\E)([\-\>])/i) {
            my ($l, $found, $r) = ($1 || "", $2, $3 || "");
            my $out = $l . $found . $r;
            my $in  = $l . $prefer . $r;
            $hack =~ s/\Q$out\E/$in/g;
            if (++$iloop > 100) {
                # Just to be safe. Not much sucks more than an infinite regexp
                push @errs, "Likely i-loop substituting alias '$ali' for '$prefer'";
                last;
            }
        }
    }
    # Ok, now break up the string into pieces and see if we can recognize them
    my @words = split(/\-/, $hack);
    my ($tok, $one) = ([],[]);
    for my $w (0..$#words) {
        
    }
}

sub map_amino_acid {
    my $self = shift;
    my ($in, $trg, $src) = shift;
    return "" unless ($in);
    my $info = $self->aa_info();
    $in = $self->_clean_aa_input( lc($in) );
    return "" unless (exists $info->{MAP}{ $in });
    my $hash = $info->{MAP}{ $in };
    my @keys = keys %{$info->{MAP}{$in}};
    if ($src) {
        # The user is explicitly defining the type of object
        $src =~ s/\s+//g;
        @keys = (lc($src));
    }
    my %out;
    foreach my $key (@keys) {
        map { $out{$_} ||= $_ } @{$hash->{$key}} if (exists $hash->{$key});
    }
    my @rv;
    if ($trg) {
        # Request for a particular part of the protein
        my %u;
        $trg = lc($trg);
        foreach my $dat(values %out) {
            $u{$dat->{$trg} || ""} = 1 if (exists $dat->{$trg});
        }
        delete $u{""};
        @rv = sort keys %u;
    } else {
        # Just return the hash data if no specific keys are provided
        @rv = values %out;
    }
    return @rv;
}

sub _clean_aa_input {
    my $self = shift;
    my $text = shift;
    $text =~ s/\P{IsASCII}//g; # Remove non-ascii characters
    $text =~ s/^\s+//; $text =~ s/\s+$//; # Lead / trail whitespace
    $text =~ s/\s+/ /g; # Whitespace runs to individual spaces
    return $text;
}

sub aa_info {
    my $self = shift;
    unless ($self->{AAI}) {
        my @errs;
        my $info = $self->{AAI} = { ALIAS => {} };
        foreach my $line (split(/[\n\r]/, $aminoAcidData)) {
            $line = $self->_clean_aa_input( $line );
            next if (!$line || $line =~ /^\#/);
            my $entry = {};
            foreach my $data (split(/\s*\;\s*/, $line)) {
                if ($data =~ /^([^\:]+)\:(.+)$/) {
                    my ($tag, $val) = (lc($1), $2);
                    $tag =~ s/\s+//g;
                    $info->{TAGS}{$tag}++;
                    push @{$info->{MAP}{lc($val)}{$tag}}, $entry;
                    if ($tag eq 'alias') {
                        # Allow multiple aliases per entry
                        push @{$entry->{$tag}}, $val;
                    } elsif ($entry->{$tag}) {
                        push @errs, "Multiple tags for '$tag' in '$data'";
                    } else {
                        $entry->{$tag} = $val;
                    }
                } else {
                    push @errs, "Tag:Value not seen in '$data'";
                }
            }
            if (exists $entry->{alias}) {
                my $val =
                    exists $entry->{token} ? $entry->{token} :
                    exists $entry->{3} ? $entry->{3} : "";
                map { $info->{ALIAS}{lc($_)}{$val} = 1 } @{$entry->{alias}}
                if ($val);
            }
        }
        while ( my ($ali, $prefer) = each %{$info->{ALIAS}}) {
            my @u = sort keys %{$prefer};
            if ($#u == 0) {
                $info->{ALIAS}{$ali} = $prefer;
            } else {
                delete $info->{ALIAS}{$ali};
                push @errs, "Multiple values for alias '$ali' : ".
                    join(", ", @u);
            }
        }
        $self->err("Errors parsing amino acid data", @errs)
            unless ($#errs == -1);
    }
    return $self->{AAI};
}

1;
