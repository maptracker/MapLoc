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
use BMS::ForkCritter;
use BMS::TableReader;
# use Statistics::Distributions;
use Bio::PrimarySeq;

my ($lookup, %populations, $localML, $header, $fc, %stuff, %mafValidate,
    $currentRow, $locCategory, $locCatName, %barcodeCache);

my %roman = ( I => 1, II => 2, III => 3, IV => 4, V => 5, VI => 6, VII => 7);
my $lastTime = 0;
my %builds = ( 36 => 'NCBI36',
               37 => 'GRCh37' );

my %wgetExcl = map { $_ => 1 } ('bcr','cgcc','.','..','lost+found','','','',);

=head1 MAF Column Normalization

The program utilizes a hash reference called $standardColumns. This
hash is critical, as it normalizes, as best as possible, the column
headers across the different MAF files. The hash is keyed to column
names as they appear in actual files, and has values that represent
the standard keys used by this program. Any of these names are
arbitrary, but standardization is required to allow the program to
reliably find required data. The most critical columns are:

  Location: Build Chr ChrStart ChrEnd ChrStr
   Alleles: AlleleRef AlleleNorm1 AlleleNorm2 AlleleTum1 AlleleTum2 Type
 Frequency: CountTum1 CountTum2 VafTum CountNorm1 CountNorm2 VafNorm
    Sample: BarcodeNorm BarcodeTum

Additionally, the get_standard_row() method is used to normalize
values. For example VAF values are sometimes reported as a fraction,
sometimes as a percentage. Null values may be empty strings, -, --,
---, ., x, and NA.

When new MAF files are encountered, it is possible that either the key
hash or the normalization method will need to be updated. You can
check each file by using the -check option:

   load_TCGA_maf.pl -check some_maf_file_I_just_downloaded.maf -getrow 10

The above command will check the specified MAF file for unknown column
names. It will also look for rudimentary consistency problems, and
indicate if read counts are available. the -getrow parameter is
optional, it specifies the number of sample rows to recover. They will
be written to a summary file.

There is also the $ignoreCol hash reference. This stores columns that
are being ignored for the purposes of MapLoc, preventing their being
flagged as non-standard when running checks.

=cut

my $standardColumns = {
    "Chromosome"                    => "Chr",
    "CHROM"                         => "Chr",
    "Start_position"                => "ChrStart",
    "Start_Position"                => "ChrStart",
    "START_POSITION"                => "ChrStart",
    "End_position"                  => "ChrEnd",
    "End_Position"                  => "ChrEnd",
    "END_POSITION"                  => "ChrEnd",
    "Strand"                        => "ChrStr",
    "STRAND"                        => "ChrStr",
    "NCBI_Build"                    => "Build",
    "NCBI_BUILD"                    => "Build",
    "Variant_Type"                  => "Type",
    "VARIANT_TYPE"                  => "Type",

    "Match_Norm_Seq_Allele1"        => "AlleleNorm1",
    "Match_Norm_Seq_Allele2"        => "AlleleNorm2",
    "MATCH_NORM_SEQ_ALLELE1"        => "AlleleNorm1",
    "MATCH_NORM_SEQ_ALLELE2"        => "AlleleNorm2",
    "Tumor_Seq_Allele1"             => "AlleleTum1",
    "Tumor_Seq_Allele2"             => "AlleleTum2",
    "TUMOR_SEQ_ALLELE1"             => "AlleleTum1",
    "TUMOR_SEQ_ALLELE2"             => "AlleleTum2",
    "Reference_Allele"              => "AlleleRef",
    "REFERENCE_ALLELES"             => "AlleleRef",
    "REFERENCE_ALLELE"              => "AlleleRef",
    "tumor_depth"                   => "DepthTum",
    "TTotCov"                       => "DepthTum",
    "TVarCov"                       => "CountTum2",
    "normal_depth"                  => "DepthNorm",
    "NTotCov"                       => "DepthNorm",
    "NVarCov"                       => "CountNorm1",
    "rna_depth"                     => "DepthRNA",

    "t_ref_count"                   => "CountTum1",
    "t_alt_count"                   => "CountTum2",
    "i_t_ref_count"                 => "CountTum1",
    "i_t_alt_count"                 => "CountTum2",
    "n_ref_count"                   => "CountNorm1",
    "n_alt_count"                   => "CountNorm2",
    "i_n_ref_count"                 => "CountNorm1",
    "i_n_alt_count"                 => "CountNorm2",
    "i_tumor_f"                     => "VafTum",
    "TumorRefReads_WU"              => "CountTum1",
    "TumorVarReads_WU"              => "CountTum2",
    "NormalRefReads_WU"             => "CountNorm1",
    "NormalVarReads_WU"             => "CountNorm2",
    "TumorVAF_WU"                   => "VafTum",
    "NormalVAF_WU"                  => "VafNorm",
    "RNAVAF_WU"                     => "VafRNA",
    "RNARefReads_WU"                => "CountRna1",
    "RNAVarReads_WU"                => "CountRna2",

    "normal_vaf"                    => "VafNorm",
    "tumor_vaf"                     => "VafTum",
    "rna_vaf"                       => "VafRNA",
    "transcript_name"               => "Transcript",
    "TranscriptID"                  => "Transcript",
    "Transcript_Id"                 => "Transcript",
    "Transcript"                    => "Transcript",
    "amino_acid_change"             => "Token-Protein",
    "Protein_Change"                => "Token-Protein",
    "Prot_String_Short"             => "Token-Protein",
    "AAChange"                      => "Token-Protein",
    "cDNA_Change"                   => "Token-cDNA",
    "c_position"                    => "Token-cDNA",
    "Genome_Change"                 => "Token-Genome",
    "ChromChange"                   => "Token-Genome",
    "Refseq_mRNA_Id"                => "Transcript",
    "transcript_name_WU"            => "TranscriptWU",

    "Validation_Status"             => "StatusValid",
    "VALIDATION_STATUS"             => "StatusValid",
    "Verification_Status"           => "StatusVerify",
    "VERIFICATION_STATUS"           => "StatusVerify",
    "condel"                        => "PredictCondel",
    "polyphen"                      => "PredictPolyPhen",
    "sift"                          => "PredictSIFT",
    "Mutation_Status"               => "PredictTCGA Status",
    "MUTATION_STATUS"               => "PredictTCGA Status",
    "MA_Func.Impact"                => "PredictMutationAssessor",
    "MA_FI.score"                   => "PredictMA Functional Impact",
    "MA_VC.score"                   => "PredictMA Variant Conservation",
    "MA_VS.score"                   => "PredictMA Variant Specificity",
    ""                   => "",
    "i_judgement"                   => "Judgement",

    "Matched_Norm_Sample_Barcode"   => "BarcodeNorm",
    "Tumor_Sample_Barcode"          => "BarcodeTum",
    "MATCH_NORM_SAMPLE_BARCODE"     => "BarcodeNorm",
    "TUMOR_SAMPLE_BARCODE"          => "BarcodeTum",
    "MATCH_NORMAL_SAMPLE_ID"        => "BarcodeNorm",
    "MATCH_NORM_SAMPLE_ID"          => "BarcodeNorm",
    "TUMOR_SAMPLE_ID"               => "BarcodeTum",
    "dbSNP_RS"                      => "dbSNP",
    "DBSNP_RS"                      => "dbSNP",

    "Sequence_Source"               => "SeqSource",
    "Sequencing_Phase"              => "SeqPhase",
    "Sequencer"                     => "SeqTech",
    "Center"                        => "SeqGroup",
    "GSC_Center"                    => "SeqGroup",
    "CENTER"                        => "SeqGroup",
    "filter"                        => "Filter",

    "Hugo_Symbol"                   => "Symbol",
    "HUGO_SYMBOL"                   => "Symbol",
    ""                          => "",
    ""                          => "",
};

# Malformed header in OV_genome.wustl.edu_000.026.000.000_ABI-26-0-0-somatic.maf
$standardColumns->{"Hugo_Symbol     Entrez_Gene_Id  GSC Center      NCBI_Build      Chromosome      Start_position  End_position    Strand  Variant_Classification  Variant_TypeReference_Allele Tumor_Seq_Allele1       Tumor_Seq_Allele2       dbSNP_RS        dbSNP_Val_Status  "} = "Symbol";

my $args = BMS::ArgumentParser->new
    ( -limit     => 0,
      -nocgi     => 1,
      -fork      => 20,
      -verbose   => 1,
      -mirrordir => "/tmp/TCGA",
      -base      => "MapLoc-TCGA",
      -build     => 'GRCh37',
      -progress  => 120,
      -aspercent => 1,
      -binsize   => 10,
      -paramfile  => 'BMS::SnpTracker::MapLoc',
      );

$args->shell_coloring( );

my $vb        = $args->val(qw(vb verbose)) || 0;
my $tm        = $args->val(qw(tm testmode trial)) || 0;

my $forkNum   = $args->val(qw(forknum fork)) || 1;
my $prog      = $args->val(qw(prog progress)) || 300;
my $limit     = $args->val(qw(limit));
my $deferCache = $args->val(qw(defer defercache defersize)) || 300;
my $pfile;

my $clobber   = $args->val(qw(clobber));
my $killOld   = $args->val(qw(clearold killold));
my $mafDir    = $args->val(qw(mafdir)) || ""; $mafDir =~ s/\/\s*$//;
my $remapSize = 120;
my $blastSz   = 120;

my $trgDir    = $args->val(qw(mirrordir)) || "/tmp/TCGA";
$trgDir =~ s/\/$//;
my $popRepDir = "$trgDir/popReports";
my $checkFh;
my $checkFile = "$trgDir/ColumnCheck.txt";
my $wgetLog   = "wget_log.txt";
my $srcUrl    = $args->val(qw(source)) || "https://tcga-data.nci.nih.gov".
    "/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor";

my $noDB     = $args->val(qw(nodb));
my $dbInst   = $args->val(qw(instance));
my $dumpCols = $args->val(qw(dumpcol dumpcols));
$dumpCols    = [split(/\s*[\t\r\n\,]+\s*/, $dumpCols)] 
    if ($dumpCols && !ref($dumpCols));
my $doDump   = $args->val(qw(dump)) || ($dumpCols ? 1 : 0);

my $clearPop = $args->val(qw(clearpop));
my $loadPops = $args->val(qw(loadpop loadpops dopop));
if ($clearPop) {
    $clearPop = {};
    $forkNum = 1;
}
my $studyFilter  = &parameter_to_hash('study','studies');
my $centerFilter = &parameter_to_hash('center','centers');
my $geneFilter   = &parameter_to_hash('gene','genes');

my $shellPreCol = "\033["; #1;33m
my $shellProCol = "\033[0;49m";

my @ignoring = qw
    (BAM_file BAM_File Entrez_Gene_Id ENTREZ_GENE_ID GENE_LIST
     Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2
     Tumor_Validation_Allele1 Tumor_Validation_Allele2
     TUMOR_VALIDATION_ALLELE1 TUMOR_VALIDATION_ALLELE2
     MATCH_NORM_VALIDATION_ALLELE1 MATCH_NORM_VALIDATION_ALLELE2
     transcript_status dbSNP_Val_Status DBSNP_VAL_STATUS
     Validation_Method VALIDATION_METHOD Algorithm
     Variant_Classification VARIANT_CLASSIFICATION 
     all_domains_WU deletion_substructures_WU categ_ignoring_null_categ 
     domain_WU context65 context_orig gene_name gene_name_WU 
     transcript_source_WU transcript_species_WU transcript_status_WU 
     transcript_version_WU start_WU stop_WU x type_WU domain
     ucsc_cons ucsc_cons_WU variant_WU reference_WU type patient_name
     gene chromosome_name_WU categ start end chr patient
     tum_allele1 tum_allele2 ref_allele newbase classification
     MA_MSA MA_PDB dataset Score Comments ASSAY

     Exon COSMIC_Gene dbSNPPopFreq
     Genome_Plus_Minus_10_Bp Genome_Plus_Minus_10_Bp

     pox qox pox_cutoff isArtifactMode oxoGCut transcript_error

     GO_Biological_Process GO_Cellular_Component GO_Molecular_Function

     SwissProt_acc_Id SwissProt_entry_Id Transcript_Exon Transcript_Position Transcript_Strand UniProt_AApos UniProt_Experimental_Info UniProt_Natural_Variations UniProt_Region UniProt_Site gc_content ref_context Other_Transcripts FamilialCancerDatabase_Syndromes DNARepairGenes_Role Description DrugBank  CCLE_ONCOMAP_overlapping_mutations CCLE_ONCOMAP_total_mutations_in_gene CGC_Mutation_Type CGC_Other_Diseases CGC_Translocation_Partner CGC_Tumor_Types_Germline CGC_Tumor_Types_Somatic Annotation_Transcript Drug_Target

     COSMIC_Codon COSMIC_fusion_genes COSMIC_overlapping_mutations COSMIC_tissue_types_affected COSMIC_total_alterations_in_gene


     TCGAscape_Amplification_Peaks TCGAscape_Deletion_Peaks Tumorscape_Amplification_Peaks Tumorscape_Deletion_Peaks Refseq_prot_Id MUTSIG_Published_Results
     OREGANNO_ID OREGANNO_Values Codon_Change MUTSIG_Significant_Genes

     i_ACHILLES_Top_Genes i_CCLE_ONCOMAP_overlapping_mutations i_CCLE_ONCOMAP_total_mutations_in_gene i_CCLE_SEQ_overlapping_mutations i_CCLE_SEQ_total_mutations_in_gene i_TCGAscape_Amplification_Peaks i_TCGAscape_Deletion_Peaks i_failure_reasons i_normal_best_gt i_init_t_lod i_t_lod_fstar

     amino_acid_change_WU annotation_errors_WU c_position_WU trv_type_WU strand_WU transcript_error_WU

     TCGAscape_Deletion_Peaks TCGAscape_Amplification_Peaks

     CCLE_ONCOMAP_overlapping_mutations CCLE_ONCOMAP_total_mutations_in_gene CGC_Mutation_Type CGC_Other_Diseases CGC_Translocation_Partner CGC_Tumor_Types_Germline CGC_Tumor_Types_Somatic COSMIC_fusion_genes COSMIC_overlapping_mutations COSMIC_tissue_types_affected COSMIC_total_alterations_in_gene FamilialCancerDatabase_Syndromes GO_Biological_Process GO_Cellular_Component  GO_Molecular_Function

     Other_Transcripts Tumorscape_Deletion_Peaks Tumorscape_Amplification_Peaks TCGAscape_Amplification_Peaks TCGAscape_Deletion_Peaks

     validation_alt_allele validation_method validation_tumor_sample

     Tumor_Sample_UUID Matched_Norm_Sample_UUID

     transcript_species transcript_source transcript_version trv_type all_domains deletion_substructures

     );

my $ignoreCol = { map { $_ => 1 } @ignoring };


# These appear to be duplicate column headers?
# OV_genome.wustl.edu_003.001.000.000_IlluminaGA-DNASeq-Level-3-1.maf
map { $ignoreCol->{$_} = 1 } qw(chromosome_name stop reference variant );

# Duplicated columns
map { $ignoreCol->{$_} = 1 } qw( validation_status strand );



&init_notes();


# &lookup_tables(); die;
my $doneSomething = 0;
if ($args->val(qw(init))) {
    &maploc();
    $doneSomething++;
}
if ($args->val(qw(mirror))) {
    &mirror();
    $doneSomething++;
}
if (my $dir = $args->val(qw(directory folder dir)) ) {
    &load_dir( $dir );
    $doneSomething++;
}
if (my $file = $args->val(qw(maf file input load)) ) {
    &load_maf( $file );
    $doneSomething++;
}
if (my $file = $args->val(qw(check)) ) {
    if (-d $file) {
        $file =~ s/\/$//;
        my $prob = 0;
        foreach my $f (&maf_files_in_dir( $file )) {
            $prob += &check_maf_file( $f );
        }
        $args->msg("[!]", "A total of $prob files had issues") if ($prob);
    } else {
        &check_maf_file( $file );
    }
       
    $doneSomething++;
}

if ($checkFh) {
    close $checkFh;
    $args->msg("[>]", "Column report generated", $checkFile);
}
&finish();

if ($doneSomething) {
} else {
    &usage();
}

sub benchmark_table {
    my $obj = shift;
    return "" unless ($obj);
    return $obj->showbench( -shell => 1, -minfrac => 0.0001 );
}

sub checkfh {
    unless ($checkFh) {
        $args->assure_dir($checkFile, 'isFile');
        open($checkFh, ">$checkFile") || $args->death
            ("Failed to write column check file", $checkFile, $!);
    }
    return $checkFh;
}

sub table_reader {
    my $file = shift || "";
    my $scan = shift || "NCBI";
    my $tr     = BMS::TableReader->new
        ( -limit     => $limit,
          -colmap    => $standardColumns,
          -scan      => $scan,
          -format    => 'tsv',
          -hasheader => 1, );
    

    map { $tr->ignore_column_name( $_ ) } @ignoring 
        unless ($args->val(qw(noignore)));

    $tr->format_from_file_name($file);
    $tr->input($file);
    return $tr;
}
sub init_notes {
    $args->msg("[!]","-testmode is on", "Most data will not be written") 
        if ($tm);
    $args->msg("[INFO]", "Using instance $dbInst",
               "Mirrored data will be read from / written to:", "   $trgDir");
}

sub finish {
    return unless ($pfile);
    if (open(PFILE, "<$pfile")) {
        my %u;
        while (<PFILE>) {
            s/[\n\r]+$//;
            my ($ptxt, $ov, $num) = split(/\t/);
            next unless ($ptxt);
            $u{$ptxt} ||= [ 0, $ov];
            $u{$ptxt}[0] += $num if ($num);
        }
        close PFILE;
        if (open(PFILE, ">$pfile")) {
            foreach my $ptxt (sort { $u{$b}[0] <=> $u{$a}[0] ||
                                     $a cmp $b } keys %u) {
                my ($num, $ov) = @{$u{$ptxt}};
                $num ||= "Parent";
                print PFILE "[$num] $ptxt\n";
                print PFILE "  $ov\n" if ($ov);
                print PFILE "\n";
            }
            close PFILE;
        }
        $args->msg("[Populations]", $pfile);
    } else {
        $args->msg("[!]","Failed to read file", $pfile, $!);
    }
}

sub dump_row {
    my $data = shift;
    my @k    = $dumpCols ? @{$dumpCols} : sort keys %{$data};
    my ($maxLen) = sort { $b <=> $a } map { length($_ || "") } @k;
    my $txt = ("- " x 20)."\n";
    my $fmt = sprintf(" %%%ds | %%s\n", $maxLen);
    foreach my $key (@k) {
        next if ($ignoreCol->{$key});
        my $v = $data->{$key};
        $txt .= sprintf($fmt, $key, !defined $v ? '-UNDEF-' : $v eq '' ? '-EMPTY STRING-' : $v);
    }
    return $txt;
}

sub get_standard_row {
    my ($data, $nonFatal) = @_;
    while (my ($k, $v) = each %{$data}) {
        $v = $data->{$k} = "" if (!defined $v);
        if ($v eq '--' || $v eq '---') {
            $data->{$k} = "";
        } elsif ($v eq '.') {
            $data->{$k} = "";
        } elsif ($v eq '-') {
            $data->{$k} = "" unless ($k =~ /^Allele/);
        } elsif ($v eq 'x') {
            $data->{$k} = "" if ($k =~ /^Seq/);
        } elsif ($v eq 'Unknown') {
            $data->{$k} = "" if ($k =~ /Allele/);
        } elsif ($v eq 'NA') {
            $data->{$k} = "" if ($k =~ /^(Count|Vaf)/);
        } elsif ($k =~ /^Vaf/) {
            if ($v =~ /^(.+)%$/) {
                $data->{$k} = $1 / 100;
            } elsif ($v > 1) {
                # Reported as a percentage
                # Will not adjust small frequencies less than one percent!
                $data->{$k} = $v / 100;
            }
        }
    }
    # No standard for reporting mutations...
    # Force standardization here:
    if (my $ref = $data->{AlleleRef}) {
        unless ($data->{AlleleNorm1}) {
            $data->{AlleleNorm1} = $ref;
            $data->{AtypicalData}++;
        }
        unless ($data->{AlleleNorm2}) {
            $data->{AlleleNorm2} = $ref;
            $data->{AtypicalData}++;
        }
        if ($data->{AlleleTum1} && $data->{AlleleTum2} &&
            $data->{AlleleTum1} eq $data->{AlleleTum2}) {
            # I am not sure why this was done
            # Fix so that Allele1 is consistent with other formatting:
            $data->{AlleleTum1} = $ref;
            $data->{AtypicalData}++;
        }
    }
    
    if ($data->{DepthTum} && $data->{CountTum2} && !$data->{CountNorm1}) {
        # Some files define the total read depth and just the variant depth
        # COAD_hgsc.bcm.edu_002.001.004.000_IlluminaGA-DNASeq-1-somatic.maf
        $data->{CountTum1}  = 
            ($data->{DepthTum}  || 0) - ($data->{CountTum2}  || 0);
        $data->{CountNorm1} = 
            ($data->{DepthNorm} || 0) - ($data->{CountNorm2} || 0);
    }
    foreach my $typ ('Tum', 'Norm') {
        my $vaf = $data->{"Vaf$typ"};
        my ($dk, $ck1, $ck2) = ("Depth$typ","Count${typ}1","Count${typ}2");
        if (!defined $data->{$ck1} && $data->{$ck2} && $data->{$dk}) {
            # Count1 not defined, but Depth and Count2 are
            $data->{$ck1} = $data->{$dk} - $data->{$ck2};
            $data->{AtypicalData}++;
        } elsif (!defined $data->{$ck2} && $data->{$ck1} && $data->{$dk}) {
            # Count2 not defined, but Depth and Count1 are
            $data->{$ck2} = $data->{$dk} - $data->{$ck1};
            $data->{AtypicalData}++;
        } elsif (!defined $data->{$dk} && $data->{$ck1} && $data->{$ck2}) {
            # Depth not defined, but both counts are
            $data->{$dk} = $data->{$ck1} + $data->{$ck2};
            $data->{AtypicalData}++;
        }
        my $chk = $data->{$dk} || 0;
        my ($c1, $c2) = ($data->{$ck1} || 0, $data->{$ck2} || 0);
        my $tot = $c1 + $c2;
        if ($chk && $tot && $chk != $tot) {
            # Both counts were provided, but do not add up to the total
            # The counts do not add up to the total depth!!
            my $msg = "$typ counts do not add up to reported depth";
            if ( $nonFatal) {
                return ($data, $msg);
            } else {
                warn &dump_row($data);
                $args->death($msg);
            }
        }
        # warn "$typ : $tot = $c1 + $c2 = $chk\n";
        if ($tot) {
            # We can calculate VAF from counts
            my ($calc, $problem) = &report_vs_calc($vaf, $c2, $tot);
            if ($problem) {

                # Most MAF files are fine, but there are some
                # with many inexplicable differences between the
                # reported counts and the reported VAFs. This used to
                # be a fatal error, but for now I think I will simply
                # clear the VAF calculation and just report the
                # alleles without frequencies
                
                delete $data->{"Vaf$typ"};
                $data->{AtypicalData}++;
                my $bin = int(10 * $problem) || 1;
                $mafValidate{$bin}++;
            } else {
                $data->{"Vaf$typ"} = $calc;
                $mafValidate{$vaf ? "" : 0}++;
            }
        }
    }
    if (my $str = $data->{ChrStr}) {
        if ($str eq '-' || $str eq '-1') {
            my $msg = "Unanticipated -1 strand assignment";
            if ( $nonFatal) {
                return ($data, $msg);
            } else {
                warn &dump_row($data);
                $args->death($msg);
            }
        }
    } elsif ($nonFatal) {
        return ($data, "Explicit strand designation not found");
    }
    if ($doDump) {
        print &dump_row($data)
            if ($doDump eq '1' ||
                (exists $data->{$doDump} && $data->{$doDump}));
    }
    return $currentRow = $data;
}

sub report_vs_calc {
    my ($vaf, $c2, $tot, $recurse) = @_;
    my $calc = int(0.5 + 1000 * $c2 / $tot) / 1000;
    # If we only have counts, return the calculated value
    return ($calc) unless (defined $vaf && $vaf ne '');

    # Validate that the counts match the reported frequency We will
    # consider up to "half a count" on either side of the actual
    # counts to accomodate rounding errors. If our

    return ($calc) unless ($vaf < ($c2-0.5) / $tot || 
                           $vaf > ($c2+0.5) / $tot);

    # Ok, we still are not agreeing. We will be generous and take a VAF
    # within 1% of what we calculate

    return $calc if (abs($vaf - $calc) < 0.01);

    # If a file is reporting VAFs as percentages, they will be
    # automatically corrected for values over 1, but not for smaller
    # values


    if ($calc < 0.01 && !$recurse) {
        my ($calc2, $prob2) = &report_vs_calc($vaf/100, $c2, $tot, 'RECURSE');
        if (!$prob2) {
            # warn "FIXED\n";
            return ($calc);
        }
    }
    
    # We STILL can not get the counts and the VAF to agree.

    my $delta = $calc - $vaf;
    my $c1 = $tot - $c2; warn sprintf("%d + %d = %d. %d / %d = %.4f != %.4f [%s]\n", $c1, $c2, $tot, $c2, $tot, $calc, $vaf, $delta);
    
    return ($calc, $delta);
}

sub load_dir {
    my $dir = shift;
    return unless ($dir);
    $dir   =~ s/\/$//;
    $pfile = "$popRepDir/BulkPopulationReport.txt";
    # $args->msg("Finding MAF files in directory", $dir);
    my @files = &maf_files_in_dir( $dir );
    foreach my $file (@files) {
        &load_maf($file);
    }
}

sub load_maf {
    my $file = &maf_file(shift);
    return unless ($file);
    if ($mafDir && $file !~ /\//) {
        $file = "$mafDir/$file";
    }
    $args->death("Failed to find requested file", $file) unless (-s $file);
    my $info = 
    $args->msg("[+]","Reading MAF file: ".&maf_details($file), $file,`date`);
    my $rmf = $file."-REMAP";
    $rmf .= "-LIMIT" if ($limit);
    $rmf .= ".tsv";
    my $short = $file;
    $short =~ s/.+\///;
    &lookup_tables();
    $fc ||= BMS::ForkCritter->new
        ( -init_meth   => \&init,
          -finish_meth => \&finish_fc,
          -method      => \&process_row,
          -limit       => $limit,
          -progress    => $prog,
          -verbose     => $vb );
    $fc->reset();
    $fc->set_colmap( $standardColumns );
    $fc->input_type('maf header hash clean');
    $fc->input($file);
    $fc->output_file('REMAP', ">$rmf");

    $args->assure_dir($popRepDir);
    $pfile ||= "$popRepDir/$short-PopulationReport.txt";
    unlink($pfile) unless ($stuff{PopFileCleared}++);
    $fc->output_file('POPS', ">>$pfile");

    # my $tr = &table_reader( $file );
    # &init();
    if (my $failed = $fc->execute( $forkNum )) {
        $args->err("$failed children failed to finish execution");
    }
    #$lastTime = 0;
    #while (my $data = &get_standard_row($tr)) {
    #    &process_row( $data );
    #    # print $loc->to_text();
    #    
    #    # die;
    #    if (time > $lastTime + 60) {
    #        $args->msg("[-]", "$chr:${l}_${r}" );
    #        $lastTime = time;
    #    }
    #}
    #close REMAP;
    return if ($clearPop || $loadPops);
    if (-s $rmf) {
        my @bits = map { split(/\n+/) } ( `wc -l $rmf`, `egrep -o '(High confidence|ERROR|MISS|Size Difference|Different Chromosome)' $rmf | sort | uniq -c`);
        map { s/^\s+// } @bits;
        $args->msg("[REMAP]","Build remapping status", @bits,
                   "egrep -v '(ERROR|MISS|100%)' $rmf | less");
    } else {
        unlink( $rmf );
    }
    
    $args->msg("Load finished", `date`);
}

sub init {
    $localML = &maploc();
}

sub finish_fc {
    my $ptxt = "";
    $localML->add_deferred_alleles() if ($localML);
    foreach my $pop (values %populations) {
        $ptxt .= join("\t", $pop->to_one_line(), 
                      $pop->simple_value("Overview"),
                      $pop->{RECORD_COUNT} || 0)."\n";
    }
    $fc->write_output('POPS', $ptxt);
    print STDERR &benchmark_table( $localML ) if
        ($fc->child() == 1 &&
         $args->val(qw(benchmark showbench dobench bench)) );
    if ($clearPop) {
        die "SAFETY - not clearing populations!!";
        my @popids = sort { $a <=> $b } keys %{$clearPop};
        my $num = $#popids + 1;
        $args->msg("Purging population data [$num]", join(', ', @popids));
        my $sth = $localML->dbh->prepare
            ( -name => "Delete population data",
              -sql  => "DELETE FROM allele WHERE pop_id = ?");
        foreach my $pid (@popids) {
            if ($tm) {
                warn $sth->pretty_print($pid);
            } else {
                $sth->execute($pid);
            }
        }
    }
    
}

sub fc_msg {
    return if ($fc && $fc->child != 1);
    $args->msg_once(@_);
}

sub process_row {
    my ($row) = @_;
    my $data = &get_standard_row($row);
    return if ($geneFilter && !$geneFilter->{uc($data->{Symbol} || "")});
    # warn &dump_row( $data); return;
    # die $args->branch($fc->{COL_MAP});
    my ($chr, $s, $e, $b, $type) = map
    { $data->{$_} } qw(Chr ChrStart ChrEnd Build Type);
    return unless ($chr && $s && $e);
    $b ||= "-UNDEFINED-";
    my $build = $builds{$b};
    unless ($build) {
        warn &dump_row( $data);
        $args->death("Failed to find build token for build '$b'");
    }
    if (my $stat = $data->{'PredictTCGA Status'}) {
        if ($stat =~ /^(variant|germline|unknown|none)$/i) {
            # Do not want these - right?
            return;
        } elsif ($stat =~ /^(loh|somatic)$/i) {
            # What we want
        } else {
            $args->death("Unrecognized status code", "'$stat'");
        }
    }
    my ($l, $r) = ($s - 1, $e + 1);
    if ($type eq 'INS') {
        ($l, $r) = ($s, $e);
        unless ($e == $s + 1) {
            $args->msg("[!]", "$type location does not properly define flanks", "$chr:$s..$e");
            return;
        }
    }
    if ($clearPop || $loadPops) {
        $args->msg_once("Clearing population data!!!!") if ($clearPop);
        &set_freqs(undef, $data);
        return;
    }
    my $loc = $localML->flanks_to_location( $chr, $l, $r, $build );
    unless ($build eq 'GRCh37') {
        #     &update_build($loc, $data);
    }
    # &fc_msg("Only mapping build over"); return;
    &set_freqs($loc, $data);
    &set_meta($loc, $data);
    $loc->add_category( $locCategory );
    if (!$tm) {
        $loc->write();
    } elsif ($tm ne '1') {
        print $loc->to_one_line()."\n";
    }
}

sub set_meta {
    my ($loc, $data) = @_;
    my %params;
    while (my ($t, $v) = each %{$data}) {
        next if (!defined $v || $v eq '');
        if ($t =~ /^Predict(.+)$/) {
            $params{"$1 Prediction"} = $v;
        } elsif ($t =~ /^Seq(Group|Phase|Tech)/) {
            my $typ = $1;
            if ($typ eq 'Phase') {
                if ($v =~ /([IV]+)$/) {
                    $v = $roman{$1} || $args->death
                        ("Failed to interpret $t for '$data->{$t}'");
                }
            }
            $params{"TCGA $t"} = $v;
        }
    }
    while (my ($t, $v) = each %params) {
        $loc->tag_values($t, $v);
    }
}

sub set_freqs {
    my ($loc, $data, $invert) = @_;
    my @errs;
    my $normSamp = $data->{"BarcodeNorm"};
    my $ml = $loc ? $loc->maploc() : $localML;
    my $su = $ml->sequtils();
    # &dump_row($data);
    foreach my $typ ('Norm', 'Tum') {
        my $ntype = $typ eq 'Norm' ? 'Normal' : 'Tumor';
        my $bc = $data->{"Barcode$typ"};
        unless ($bc) {
            push @errs, "No $ntype barcode";
            next;
        }
        my $bcdat = &parse_barcode( $bc, $ntype, $normSamp );
        next if ($loadPops);
        unless ($bcdat && $bcdat->{Pop}) {
            push @errs, "No $ntype population generated";
            next;
        }
        my $popid = $bcdat->{Pop}->pkey();
        if ($clearPop) {
            $clearPop->{$popid}++;
            next;
        }

        my ($a1, $a2) = ($data->{"Allele${typ}1"}, $data->{"Allele${typ}2"});
        my $f1  = $data->{"Vaf$typ"};
        $f1     = undef unless (defined $f1 && $f1 ne "");
        my @adat = ([$a1, $f1]);
        if ($a2 && $a2 ne $a1) {
            my $f2 = defined $f1 ? 1 - $f1 : undef;
            push @adat, [$a2, $f2];
        }
        my $depth = $data->{"Depth$typ"} || undef;
        next if ($tm);
        my @keep;
        foreach my $ad (@adat) {
            my ($base, $otherFrac) = @{$ad};
            if ($base =~ /^\-+([A-Z]+)\-*$/ || $base =~ /([A-Z]+)\-+$/) {
                # eg "-G" at 1.GRCh37:215408280..215408281 TCGA-A8-A08L-01
                $loc->tag_values("Caution", "Corrected allele '$base'");
                $base = $1;
                $args->msg_once("Correcting Allele '$ad->[0]' to $base");
            }
            unless ($base eq '-' || 
                    $base eq $su->cached_safe_dna( $base )) {
                push @errs, "Illegal allele base '$base'";
                print &dump_row($currentRow);
                next;
            }
            $base = $ml->sequtils->revcom($base) if ($invert);
            if ($depth && defined $otherFrac) {
                my $frac = 1 - $otherFrac;
                $loc->add_allele_to_db($base, $popid, $frac, $depth, 'copy');
            } else {
                $loc->add_allele_to_db($base, $popid, undef, undef, 'copy');
            }
            push @keep, $base;
        }
        if ($killOld) {
            my @rem = $loc->clear_alleles_from_db( $popid, \@keep );
            unless ($#rem == -1) {
                my @bases = map { $ml->cached_pkey_to_text( $_ ) } @rem;
                $args->msg("Removed old alleles: ".join(',', @bases),
                           $loc->to_one_line(), $bcdat->{Pop}->to_one_line());
            }
        }
    }
    $localML->add_deferred_alleles() 
        if ($localML->deferred_allele_count() > $deferCache);
    $args->msg("[FreqErr]", $loc->to_one_line(), 
               $normSamp || "Unknown Barcode", @errs) unless ($#errs == -1);
}

sub parse_barcode {
    my $bc = shift || "";
    my ($type, $normSamp) = @_;
    # https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode#TCGAbarcode-ReadingBarcodes
    my $rv = $barcodeCache{"$bc\t$type"};
    return $rv if ($rv);
    # TCGA-AB-2804-03B-01W-0728-08 / TCGA-AB-2804-11B-01W-0728-08
    # Last bit is Plate/Center

    #            [ PopName                      ]
    #                  [ TSS  ] [ PAT  ] [ SAMP] [ VIAL] [ PORT] [ ANA ]
    if ($bc =~ /^(TCGA-([^\-]+)-([^\-]+)-(\d{2}))([A-Z])?(\d{2})?([A-Z])?/) {
        $localML->bench_start();
        my ($popN, $tss, $pat, $samp, $vial, $portion, $analyte) = 
            ($1, $2, $3, $4, $5, $6, $7);
        $rv = $barcodeCache{"$bc\t$type"} = {
            TSS     => $tss,  # Tissue source site
            Patient => $pat,  # Study participant
            Sample  => $samp, # Sample type
            PopName => $popN, # MapLoc population
        };
        
        if (my $luTab = &lookup_tables) {
            if (my $s = $luTab->{Sample}{$samp || ""}) {
                $rv->{SampleType} = $s->{Desc};
            } else {
                $args->msg_once("Failed to parse Sample code '$samp'");
            }
            if (my $c = $luTab->{Source}{$tss || ""}) {
                $rv->{StudySite}  = $c->{Site};
                my $name = $rv->{StudyName} = $c->{Study};
                if (my $sc = $luTab->{Study}{$name}) {
                    $rv->{StudyCode} = $sc->{Code};
                } else {
                    $args->death("Unrecognized Study name $name");
                }
            } else {
                $args->msg_once("Unrecognized TSS $tss");
            }
        }
        # warn $args->branch($rv);
        my $pop = &barcode_population( $rv, $type, $normSamp );
        $localML->bench_end();
    } else {
        $rv = $barcodeCache{"$bc\t$type"} = {};
        $args->msg_once("Unrecognized barcode $bc");
    }
    return $rv;
}

sub barcode_population {
    my ($dat, $type, $normSamp) = @_;
    return undef unless ($dat);
    my $pn = $dat->{PopName};
    return undef unless ($pn && $type);
    my $pname = "$pn $type";
    unless ($populations{$pname}) {
        my $pop = $populations{$pname} = $localML->get_population
            ( $pname, "TCGA Samples");
        my $sc  = $dat->{StudyCode} || "Unknown";
        my $par = "TCGA $sc $type";
        my $parPop = $populations{$par};
        unless ($parPop) {
            $parPop = $populations{$par} = $localML->get_population
                ( $par, "TCGA Samples");
            map {$parPop->tag_values($_, $dat->{$_})} qw(StudyName StudySite);
            my $gpar = "TCGA $type";
            unless ($populations{$gpar}) {
                my $gpop = $populations{$gpar} = $localML->get_population
                    ( $gpar, "TCGA Samples");
                my $ggpar = "TCGA Populations";
                unless ($populations{$ggpar}) {
                    $populations{$ggpar} = $localML->get_population
                        ( $ggpar, "TCGA Samples");
                }
                $gpop->parent($populations{$ggpar});
                $gpop->update();
            }
            $parPop->parent($populations{$gpar});
            $parPop->update();
            # $args->msg("[PAR]", $parPop->to_text());
        }
        $pop->parent($parPop);
        my %params = %{$dat};
        my $nm = $params{StudyName};
        $nm = $nm ? "the $nm" : "an unknown";
        my $tp = "$type sample";
        if (my $st = $params{SampleType}) {
            $tp = "$st sample";
            $tp = "$type sample matched against $st"
                if ($type eq 'Normal' && $st !~ /normal/i);
        }
        $params{Category} = $locCatName;
        $params{Overview} = "This is a $tp from $nm TCGA study. Frequency counts represent read depth covering this location.";
        if ($type ne 'Normal') {
            if ($normSamp) {
                my $nPop = &parse_barcode($normSamp, 'Normal');
                $params{"Normal Population"} = $nPop->{PopName} . " Normal";
                $pop->clear_tag('Normal Population');
            }
            $params{Tally} = "TCGA";
        } else {
            $pop->clear_tag('Tally');
            #$params{Tally} = "TCGA";
        }
        while (my ($t, $v) = each %params) {
            $pop->tag_values($t, $v) if (defined $v && $v ne '');
        }
        $pop->update();
        # $args->msg("[POP]", $pop->to_text());
    }
    my $pop = $dat->{Pop} = $populations{$pname};
    $pop->{RECORD_COUNT}++;
    return $pop;
}


sub lookup_tables {
    return $lookup if ($lookup);
    $lookup = {};
    $args->assure_dir("$trgDir");
    my $url = "https://tcga-data.nci.nih.gov/datareports/codeTablesExport.htm";

    my $tabMap = {
        sampleType       => "Sample",
        tissueSourceSite => "Source",
        diseaseStudy     => "Study",
        centerCode       => "Center",
        # bcrBatchCode     => "BCR",
        # dataLevel        => "Level",
        # platformCode     => "Platform",
        # portionAnalyte   => "Analyte",
        # tissue           => "Tissue",
    };
    my $colMap = {
        Sample => {
            'Code'               => 'Code',
            'Definition'         => 'Desc',
            'Short Letter Code'  => 'Short',
        },
        Source => {
            'TSS Code'           => 'Code',
            'Source Site'        => 'Site',
            'Study Name'         => 'Study',
            'BCR'                => 'BCR',
        },
        Study => {
            'Study Abbreviation' => 'Code',
            'Study Name'         => 'Name',
        },
        BCR => {
            'BCR Batch'          => 'Code',
            'Study Name'         => 'Study',
        },
        Center => {
            'Code'               => 'Code',
            'Center Name'        => 'URL',
            'Center Type'        => 'Type',
            'Display Name'       => 'Name',
            'Short Name'         => 'Short',
        },
    };
    my $luKeys = {
        Study  => ['Name', 'Code'],
    };
    my @msg;
    foreach my $ctr (sort keys %{$tabMap}) {
        my $file = "$trgDir/$ctr.tsv";
        my $post = "exportType=tab&codeTablesReport=$ctr";
        &wget_location($url, $trgDir, $file, $post);
        my $tab = $tabMap->{$ctr} || "";
        my $tr = BMS::TableReader->new
        ( -file      => $file,
          -colmap    => $colMap->{$tab},
          -format    => 'tsv',
          -hasheader => 1, );
        my $dtrg = $lookup->{$tab} ||= {};
        my @keyz = @{$luKeys->{$tab} || ['Code']};
        while (my $data = $tr->next_clean_hash()) {
            map { $dtrg->{$data->{$_}} = $data } @keyz;
        }
        my @k = keys %{$dtrg};
        push @msg, sprintf("%s : %d", $tab, $#k + 1);
    }
    $args->msg("[<]", "Code Lookup entries: ".join(', ', @msg));
    
    return $lookup;
}

sub maploc {
    my $ml = BMS::SnpTracker::MapLoc->new
        ( -instance => $dbInst,
          -noenv  => $args->val(qw(noenvironment noenv)),
          -makedb => $args->val(qw(makedatabase makedb)) );
    $locCatName = "TCGA Mutations";
    my $cat = $ml->get_text( $locCatName );
    $cat->tag("Description", 
              "Sites observed in at least one TCGA cancer sample");

    $cat->tag("SNP Class Tag", "StudyCode");

    $cat->tag("Useful Tag", "Condel Prediction", "PolyPhen Prediction", "SIFT Prediction", "TCGA Status Prediction", "Mutation Assessor Prediction", "MA Functional Impact Prediction", "MA Variant Conservation Prediction", "MA Variant Specificity Prediction", "TCGA SeqSource", "TCGA SeqPhase", "TCGA SeqTech", "TCGA SeqGroup", "MutationAssessor", "Mapped From", "TCGA Build Remapping");

    $cat->tag("URL Tag", 'tag="TCGA SeqGroup" url="http://TAGVAL"');
    $cat->tag("URL Tag", 'tag="TSS" url="https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=Tissue%20Source%20Site"');
    $cat->clear_tag("URL Tag", 'tag="Study Code" url="https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=Disease%20Study"');
    $cat->tag("URL Tag", 'tag="StudyCode" url="https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=Disease%20Study"');

    $locCategory = $cat->pkey();
    $cat->write();

    my $catPar = $ml->get_text( "Mutation Categories" );
    $catPar->tag("Description", 
                 "Variant locations attributed to somatic mutations");
    $catPar->tag("Member", $cat->text());
    $catPar->write();
    return $ml;
}

sub mirror {

    # I tried to do this all just using wget. It was a mess. I could
    # not get wget to honor directory exclusion wildcards (even using
    # explicit full paths). So instead I am using it to recover
    # specific files, and to read directory listings.
    $wgetLog = "$trgDir/wget_log.txt";
    unlink($wgetLog) if (-e $wgetLog);

    $args->assure_dir("$trgDir/readmes");
    $args->msg("[>]", "Mirroring TCGA MAF files to local computer",
               " -source $srcUrl", " -mirror $trgDir", "LOG: $wgetLog",
               "Files already present will be ".
               ($clobber ? "Re-downloaded" : "skipped"));
    unlink($wgetLog) if (-e $wgetLog);
    my @toRead = ($srcUrl);
    while ($#toRead > -1) {
        my $dir = shift @toRead;
        push @toRead, &wget_location($dir, $trgDir);
    }
    $args->msg("Finished");
    exit;
}

sub wget_location {
    my ($src, $tmpDir, $out, $post) = @_;
    my $isDir;
    if (!$out) {
        # Index listing
        $out = "$tmpDir/wgetIndex.html";
        unlink($out) if (-e $out);
        $isDir = 1;
    } elsif (!$clobber && -s $out) {
        $args->msg("[.]","Keeping existing file", $out);
        return ();
    }
    $post = $post ? "--post-data '$post' " : "";
    my $cmd = "wget -O $out -a $wgetLog $post$src";
    $src =~ s/\/+$//;
    my $shortSrc = $src;
    $shortSrc =~ s/\Q$srcUrl\E//;
    $shortSrc =~ s/^\///;
    $shortSrc ||= '.';
    system($cmd);
    unless (-s $out) {
        $args->msg("[!!]", "Failed to read remote site", 
                   "From: $shortSrc/", "  To: $out");
        return ();
    }
    if ($isDir) {
        open(IND, "<$out") || $args->death
            ("Failed to read local wget listing", $out, $!);
        my @rv;
        while (<IND>) {
            if (/href=\"([^\"]+)\"/) {
                my $href = $1;
                if ($href =~ /(.+)\/$/) {
                    # This is a directory
                    my $dir = $1;
                    if ($wgetExcl{$dir}) {
                        # Do not want
                        $args->msg_once("[-]", "Excluding all '$dir' directories");
                    } else {
                        # recurse
                        push @rv, "$src/$dir";
                    }
                } else {
                    # A file
                    if ($href =~ /(maf|readme)/i) {
                        # What we are looking for
                        # The file names are highly non-standard, so clean up:
                        my $name = &stnd_maf_name($shortSrc, $href);
                        my $dest = $href =~ /readme/i ? 
                            "$tmpDir/readmes/$name" : "$tmpDir/$name";
                        if (0) {
                            $args->msg("[DEBUG]", $name, $href, "$shortSrc/");
                        } else {
                            &wget_location("$src/$href", $tmpDir, $dest);
                        }
                    }
                }
            }
        }
        close IND;
        return @rv;
    } elsif ($out =~ /maf/) {
        &check_maf_file( $out );
    } else {
        $args->msg("[>]", "Downloaded: $out");
        return ();
    }
    return ();
}

sub parameter_to_hash {
    my $rv;
    foreach my $val ($args->each_split_val
                      ('/\s*[\n\r\,]+\s*/', @_)) {
        if (defined $val && $val ne '') {
            $rv ||= {};
            $rv->{uc($val)} = 1;
        }
    }
    return $rv;
}

sub maf_files_in_dir {
    my $dir = shift;
    if ($studyFilter) {
        $args->msg("[#]","Only capturing requested studies:",
                   join(' ', sort keys %{$studyFilter}));
    }
    if ($centerFilter) {
        $args->msg("[#]","Only capturing requested centers:",
                   join(' ', sort keys %{$centerFilter}));
    }
    my %found;
    my $noSrc   = $args->val(qw(nosource nocenter));
    my $oldToo  = $args->val(qw(useold));

    foreach my $file ($args->read_dir( -dir => $dir,
                                       -recurse => 1,
                                       -keep => '.maf$', )) {
        my $info = &parse_name($file);
        next if ($studyFilter  && !$studyFilter->{$info->[0]});
        next if ($centerFilter && !$centerFilter->{$info->[1]});
        my @kbits;
        push @kbits, $info->[0];
        push @kbits, $info->[1] unless ($noSrc);
        my $key  = join("\t", @kbits) || '';
        push @{$found{$key}}, $file;
    }

    # Organize results by ascending age:
    foreach my $key (sort keys %found) {
        my $arr = $found{$key};
        my @sorted;
        foreach my $file (@{$arr}) {
            my $info = &parse_name($file);
            push @sorted, [$info->[4], $file];
        }
        @sorted = map { $_->[1] } sort { $a->[0] <=> $b->[0] } @sorted;
        $found{$key} = \@sorted;
    }

    if (!$oldToo) {
        my $grpBy = $noSrc ? "study" : "study/center pair";
        $args->msg("[#]", "Only keeping newest MAF file for each $grpBy");
    } else {
        $args->msg("[#]", "Keeping all MAF files, including older ones");
    }
    my @rv;
    foreach my $key (sort keys %found) {
        my @report = @{$found{$key}};
        @report = ($report[0]) unless ($oldToo);
        push @rv, @report;
    }
    my $num = $#rv + 1;
    if ($num) {
        $args->msg("[>]","$num file".($num > 1 ? 's' : '')." recovered");
    } else {
        $args->msg("[!]","No maf files found matching your criteria!");
    }
    return @rv;
}

sub stnd_maf_name {
    my ($shortSrc, $name) = @_;
    my @sbits = split(/\//, $shortSrc);
    my $src = &get_domain_name( $sbits[-1] ) || &get_domain_name( $name );
    my $isNotMaf = ($name =~ /readme/i || $name =~ /html?$/i) ? 1 : 0;

    my $type;
    if ($shortSrc =~ /^([a-z]+)\//) {
        $type = uc($1);
    } elsif ($isNotMaf) {
        $type = "zzGlobal";
    } else {
        $args->msg("[!]", "Failed to identify tumor type in directory path",
                   $shortSrc, $name);
        return "";
    }
    my @bits;
    if ($shortSrc =~ /(\d+)\.(\d+)\.(\d+)(\.\d+)?/) {
        @bits = ($1, $2, $3, $4 || "");
        map { s/\.// } @bits;
        while (!defined $bits[-1] || $bits[-1] eq "") { pop @bits; }
    } elsif ($isNotMaf) {
        @bits = (0,0,0);
    } elsif ($shortSrc =~ /^[a-z]+\/gsc$/) {
        # This is at the root level of a project
        $src ||= 'cancergenome.nih.gov';
        @bits = (0,0,1);
    } else {
        $args->msg("[!]", "Failed to identify level / version in path",
                   $shortSrc, $name);
        return "";
    }
    $name =~ s/\Q$src\E//gi;
    $name =~ s/\Q$type\E//gi;
    $name =~ s/(tcga)//gi;
    $name =~ s/[\.\-_]+/-/g; $name =~ s/^\-//; $name =~ s/\-$//;
    $name =~ s/\-maf$/\.maf/;
    $name =~ s/\-txt$/\.txt/;
    while ($#bits < 3) { unshift @bits, 0; }
    my $vers  = join('.', map { sprintf("%03d", $_) } @bits);
    $src    ||= "www.example.com";
    return sprintf("%s_%s_%s_%s", $type, $src, $vers, $name);
}

sub maf_file {
    my $req = shift;
    return $req if (-s $req);
    return "$trgDir/$req" if (-s "$trgDir/$req");
    return $req;
}

sub check_maf_file {
    my $file  = &maf_file(shift);
    my $info  = &parse_name($file);
    my $short = $file;
    $short =~ s/.+\///;
    open (CHK, "<$file") || $args->death
        ("Failed to read MAF file", $file, $!);
    my $bar = ("#" x 70)."\n";
    my $age = sprintf("%.1f days", $info->[4]);
    my $sz  = sprintf("%.3f MB", (-s $file) / (1024 * 1024));
    &checkfh();
    print $checkFh "\n$bar" . $file."\n";
    print $checkFh "  AGE: $age SIZE: $sz  STUDY: $info->[0]\n";
    print $checkFh "$bar";
    # $args->msg("[CHECK]", $short);
    my ($head, @msg, @errs);
    my $numRows = 0;
    my $cntr = 0;
    my $getRows = $args->val(qw(numrow numrows getrow getrows)) || 1;
    my (%ignored, %errs, $noTumCount, %usedMaps);
    while (<CHK>) {
        next if (/^\#/);
        s/[\n\r]+//;
        if ($head) {
            # Data row
            $numRows++;
            last if ($numRows > $getRows);
            my @row = split(/\t/);
            print $checkFh "\n[ROW \#$numRows]\n";
            my %hash;
            for my $i (0..$#{$head}) {
                my $col  = $head->[$i];
                my $val  = $row[$i];
                my $stnd = $standardColumns->{$col};
                $usedMaps{$stnd}{$col} ||= ++$cntr if ($stnd);
                $hash{$stnd || ""} = $val;
                $val     = '-UNDEF-' if (!defined $val);
                my $vlen = CORE::length($val);
                if ($vlen > 50) {
                    $val = substr($val, 0, 50) . " ... $vlen chars";
                }
                if ($ignoreCol->{$col}) {
                    $ignored{$col} ||= $val;
                    next;
                }
                my $tok  = $stnd ? '   ' : '???';
                printf($checkFh "%40s [%s] %s\n", $stnd || $col, $tok, $val);
            }
            my ($h, $err) = &get_standard_row( \%hash, 'nonFatal' );
            $errs{$err}++ if ($err);
            $noTumCount = 1 unless (($h->{CountTum1} || 0) +
                                    ($h->{CountTum2} || 0));
            foreach my $type ('Norm', 'Tum') {
                $errs{"No $type Barcode"}++ unless ($h->{"Barcode$type"});
            }
        } else {
            $head = [split(/\t/)];
            foreach my $col (@{$head}) {
                next if ($standardColumns->{$col});
                next if ($ignoreCol->{$col});
                push @msg, $col;
            }
        }
    }
    close CHK;

    foreach my $needed (qw(Chr ChrStart ChrEnd Build Type
                           BarcodeNorm BarcodeTum 
                           AlleleNorm1 AlleleTum1 AlleleTum2)) {
        unless (exists $usedMaps{$needed}) {
            $errs{"Required standard column $needed not found"}++;
        }
    }
    print $checkFh "[NO TUMOR READ COUNTS]\n" if ($noTumCount);
    my $umTxt = "\n[COLUMN MAPPINGS USED]:\n";
    foreach my $stnd (sort keys %usedMaps) {
        my %um  = %{$usedMaps{$stnd}};
        my @raw = sort { $um{$a} <=> $um{$b} } keys %um;
        $umTxt .= sprintf("   %20s <= %s\n", $stnd, join(' & ', @raw));
        $errs{"Multiple columns assigned to $stnd"}++ if ($#raw != 0);
    }
    my @err = sort keys %errs;
    print $checkFh "$umTxt\n" unless
        ($#msg == -1 && $#err == -1 && !$noTumCount);

    my @igc = sort keys %ignored;
    if ($#igc != -1) {
        print $checkFh "\n[IGNORED COLUMNS] : sample value\n";
        foreach my $col (@igc) {
            printf($checkFh "%40s [---] %s\n", $col, $ignored{$col});
        }
    }
    if ($#msg == -1 && $#err == -1) {
        my $tag = "OK {$info->[0] $age / $sz \@$info->[1]}";
        if ($noTumCount) {
            $args->msg("[#]","$tag (but no read counts)", $short);
        } else {
            $args->msg("[+]","$tag", $short);
        }
        return 0;
    }
    if ($#err != -1) {
        $args->msg("[!]","$short - Data processing errors", @err);
        print $checkFh "\n[PROCESSING ERRORS]:\n".join
            ("", map { "      [$errs{$_}] $_\n" } @err);
    }
    if ($#msg != -1) {
        $args->msg("[-]","$short - some unrecognized columns", @msg);
        print $checkFh "\n[UNRECOGNIZED COLUMNS]:\n".join
            ("", map { "      $_\n" } @msg);
    }
    return 1;
}

sub get_domain_name {
    my $txt = shift;
    if ($txt =~ /^([a-z]{2,12}\.[a-z]{2,12}\.[a-z]{2,12}[a-z\.]+)_/) {
        my $dom = lc($1);
        while ($dom && $dom !~ /\.(com|org|edu|uk|ca|gov)$/) {
            $dom =~ s/\.[^\.]+$//;
        }
        return $dom;
    }
    return "";
}

sub parse_name {
    my $file = shift;
    my $sz   = ( -M $file );
    $file    =~ s/.+\///;
    # Study, Center, Level, FileStuff
    my @bits = split(/_/, $file);
    # DaysOld
    $bits[4] = $sz;
    map { $_ = "" unless (defined $_) } @bits;
    return \@bits;
}

sub maf_details {
    my $file = shift;
    my $info = &parse_name($file);$info->[1]/
    return sprintf("%s %.1f days / %.3f MB \@%s", $info->[0],
                   $info->[4], (-s $file) / (1024 * 1024), $info->[1]);
}

sub markup_notes {
    $args->_stop_intercept();
    foreach my $txt (@_) {
        foreach my $line (split(/\n/, $txt)) {
            my $tca;
            if ($line =~ /^\#/) {
                $tca = 31;
            } elsif ($line =~ /^\s*http/) {
                $tca = 34;
            }
            my $msg = $line;
            if ($tca) {
                $msg = sprintf("%s%sm%s%s", $shellPreCol, $tca,
                               $msg, $shellProCol);
            }
            warn "$msg\n";
        }
    }
}

sub usage {

    my $txt = <<INFO;

# This script is designed to help download TCGA MAF files, and to load
# them into a MapLoc variant database. If the prior sentance does not
# make sense to you, you should not use this program.

# Before downloading TCGA data, you should check their terms of
# use. Generally, the data are freely available, but the consortium
# has some "first publication" restrictions.

  http://cancergenome.nih.gov/publications/publicationguidelines

# To download the data, use the -mirror argument

  load_TCGA_maf.pl -mirror -mirrordir /some/place/to/store/theData

# The -mirrordir argument can be left out, and will default to:

  /tmp/TCGA

# The program will recover all MAF files from the TCGA and store them
# locally. The files will be renamed in an attempt to standardize
# them, and the directory structure will be flattened. You can specify
# the root URL using the -source parameter; The default URL for
# retrieving data is:

  $srcUrl

# As each file is recovered, it will be checked for any unexpected
# column headers and basic problems with the data. New files
# frequently have new columns, and occasionally change column
# names. If this occurs, you will need to alter the code to help the
# program recognize the location of essential data in the
# file. Comments within the top of this script provide more
# information.

# Once you have one or more MAF files, it is time to load them. You
# will of course need MapLoc to be set up; please run setupMapLoc.pl
# if you have not done so already. You then will likely want to test
# load a few entries to verify that things are working all right:

  load_TCGA_maf.pl


INFO
    
    &markup_notes( $txt );
}
