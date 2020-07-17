use strict;
use vars;
use warnings;
use Getopt::Long;
use Math::CDF qw(qnorm);

my $chr_header = "CHR";
my $pos_header = "POS";
my $ref_header = "REF";
my $alt_header = "ALT";
my $var_id_header = "VAR_ID";
my $effect_allele_header = "Effect_Allele_PH";
my $eaf_header = "EAF_PH";
my $beta_header = "BETA";
my $or_header = "ODDS_RATIO";
my $p_header = "P_VALUE";
my $n_header = "Neff";
my $neff = undef;
my $ncase = undef;
my $nctrl = undef;
my $info_header = "Info";

GetOptions(
		'var-id-header=s' => \$var_id_header,
		'chr-header=s' => \$chr_header,
		'pos-header=s' => \$pos_header,
		'ref-allele-header=s' => \$ref_header,
		'alt-allele-header=s' => \$alt_header,
		'effect-allele-header=s' => \$effect_allele_header,
		'eaf-header=s' => \$eaf_header,
		'or-header=s' => \$or_header,
		'beta-header=s' => \$beta_header,
		'p-header=s' => \$p_header,
		'n-header=s' => \$n_header,
		'info-header=s' => \$info_header,
		'n-eff=f' => \$neff,
		'n-case=f' => \$ncase,
		'n-ctrl=f' => \$nctrl,
		);

my $header = <>; 
chomp($header);
my @header = split(/\s+/, $header);
my %header = ();
map {$header{$header[$_]}=$_} 0..$#header;

print "SNP\tCHR\tPOS\tA1\tA2\tREF\tEAF\tBeta\tse\tP\tN\tINFO\n";

while (<>)
{
		chomp;
		my @line = split;
		my $var_id = $line[$header{$var_id_header}];
		my @var_id = split("_", $var_id);
		my $snp = undef;
		my $chr = undef;
		my $pos = undef;
		my $ref = undef;
		my $alt = undef;
		if (scalar(@var_id) == 4)
		{
				$snp = $var_id;
				$chr = $var_id[0];
				$pos = $var_id[1];
				$ref = $var_id[2];
				$alt = $var_id[3];
		}
		elsif (defined $header{$var_id_header} && defined $header{$var_id_header} && defined $header{$var_id_header} && defined $header{$var_id_header})
		{
				$snp = $var_id;
				$chr = $line[$header{$chr_header}];
				$pos = $line[$header{$pos_header}];
				$ref = $line[$header{$ref_header}];
				$alt = $line[$header{$alt_header}];
		} 
		else
		{
				next;
		}
		my $effect_allele = $line[$header{$effect_allele_header}];
		my $a1 = $effect_allele;
		my $a2 = $ref;
		if ($ref eq $effect_allele)
		{
				$a2 = $alt
		}
		elsif ($alt ne $effect_allele)
		{
				next;
		}
		my $eaf = defined $header{$eaf_header} ? $line[$header{$eaf_header}] : "NA";
		my $beta = undef;
		if ($header{$beta_header})
		{
				$beta = $line[$header{$beta_header}];
		}
		elsif ($header{$or_header})
		{
				$beta = log($line[$header{$or_header}]);
		}
		my $p = $line[$header{$p_header}];
		my $zscore = -qnorm($p / 2);
		my $se = abs($beta / $zscore);
		my $n = defined $header{$n_header} ? $line[$header{$n_header}] : (defined $neff ? $neff : ($ncase && $nctrl ? 4 / (1/$ncase + 1/$nctrl) : "NA"));
		my $info = defined $header{$info_header} ? $line[$header{$info_header}] : 1;

		print "$snp\t$chr\t$pos\t$a1\t$a2\t$ref\t$eaf\t$beta\t$se\t$p\t$n\t$info\n";
}

