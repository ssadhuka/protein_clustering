BEGIN{
    #push @INC,'/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.10.1/lib/5.10.1/';
    #push @INC,'/broad/software/free/Linux/redhat_6_x86_64/pkgs/perl_5.10.1/lib/site_perl/5.10.1/';
    push @INC,'/home/unix/tgreen/perl5/lib/site_perl/5.10.1/';
    #push @INC,'/home/unix/tgreen/dig-jenkins/lib/';
}


use strict;
use Getopt::Long;
use Pod::Usage;
use Config::Simple;
use DBI;
use Log::Log4perl;
#use mysql_db;
#use KeyArgs;


my $help;
my $man;
my $dbi_configuration_file;
my $debug;
my $trace;
my $autocommit = 1;
my $log_file_name='dump_property.log';
my @param = @ARGV;
my $no_phenotype;
my $parent;
my $dataset;
my $dataset_id;
my $phenotype;
my $min_pvalue;


GetOptions(
    'trace'                         => \$trace,
    'debug'                         => \$debug,
    'help'                          => \$help,
    'man'                           => \$man,
    'no-phenotype'                  => \$no_phenotype,
    'parent=s'                      => \$parent,
    'dataset=s'                     => \$dataset,
    'dataset-id=s'                  => \$dataset_id,
    'phenotype=s'                   => \$phenotype,
    'min-pvalue=f'                  => \$min_pvalue,
    'dbi-configuration-file=s'      => \$dbi_configuration_file,
    ) or pod2usage(-verbose         => 1) && exit;

pod2usage({-noperldoc => 1, -verbose => 1,-exitval=>1}) if defined $help;
pod2usage({-noperldoc => 1, -verbose => 2,-exitval=>1}) if defined $man;
pod2usage({-message=>"\n\n ERROR: --dbi-configuration-file not configured\n\n",-noperldoc => 1, -verbose => 2}) if (! -e $dbi_configuration_file );

my $log4perlconf = qq(
    log4perl.category.Foo.Bar        = INFO, Screen , Log
    log4perl.category.debugger       = DEBUG, Screen
    log4perl.category.tracer         = TRACE, Screen
    log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr  = 1
    log4perl.appender.Screen.layout= PatternLayout
    log4perl.appender.Screen.layout.ConversionPattern=  %d %p %m%n
    log4perl.appender.Log            = Log::Log4perl::Appender::File
    log4perl.appender.Log.filename   = $log_file_name
    log4perl.appender.Log.mode       = append
    log4perl.appender.Log.layout= PatternLayout
    log4perl.appender.Log.layout.ConversionPattern=  %d %p %m%n
  );

Log::Log4perl::init( \$log4perlconf );
my $logger;
if($debug){
     $logger = Log::Log4perl->get_logger('debugger');
}elsif($trace){
     $logger = Log::Log4perl->get_logger('tracer');
}else{
     $logger = Log::Log4perl->get_logger('Foo.Bar');
}

$logger->info("Starting $0 @param");

my $cfg;
my $db_string;
my $host_url;
my $db_password;
my $dbh;
my $user;

$cfg = new Config::Simple($dbi_configuration_file);
$user=$cfg->param('user');
$db_string=$cfg->param('database');
$host_url=$cfg->param('host_url');
$db_password=$cfg->param('db_pswd');
$dbh =  DBI->connect("DBI:mysql:database=${db_string};host=${host_url}",$user , $db_password, 
		     {'RaiseError' => 1, 'AutoCommit' => $autocommit});


my $sql = "SELECT META_PROP.DB_COL
from META_PROP_DATASET_PH, META_PROP 
where  DATASET like '$dataset_id'
and PH = '$phenotype' 
and META_PROP.PROP=META_PROP_DATASET_PH.PROP
and META_PROP.MEANING = 'P_VALUE'
order by META_PROP.MEANING, META_PROP.SORT limit 1";
$logger->debug($sql);
my $sth = $dbh->prepare($sql);
$sth->execute();
my $pvalue_column = $sth->fetchrow();
$sth->finish();


$sql = "SELECT META_PROP.DB_COL
from META_PROP_DATASET_PH, META_PROP 
where  DATASET like '$dataset_id'
and PH = '$phenotype' 
and META_PROP.PROP=META_PROP_DATASET_PH.PROP
and META_PROP.MEANING = 'BETA'
order by META_PROP.MEANING, META_PROP.SORT limit 1";
$logger->debug($sql);
$sth = $dbh->prepare($sql);
$sth->execute();
my $beta_column = $sth->fetchrow();
$sth->finish();

$sql = "SELECT META_PROP.DB_COL
from META_PROP_DATASET_PH, META_PROP 
where  DATASET like '$dataset_id'
and PH = '$phenotype' 
and META_PROP.PROP=META_PROP_DATASET_PH.PROP
and META_PROP.MEANING = 'ODDS_RATIO'
order by META_PROP.MEANING, META_PROP.SORT limit 1";
$logger->debug($sql);
$sth = $dbh->prepare($sql);
$sth->execute();
my $or_column = $sth->fetchrow();
$sth->finish();

$sql = "SELECT META_PROP.DB_COL
from META_PROP_DATASET_PH, META_PROP 
where  DATASET like '$dataset_id'
and PH = '$phenotype' 
and META_PROP.PROP=META_PROP_DATASET_PH.PROP
and META_PROP.MEANING = 'ZSCORE'
order by META_PROP.MEANING, META_PROP.SORT limit 1";
$logger->debug($sql);
$sth = $dbh->prepare($sql);
$sth->execute();
my $z_column = $sth->fetchrow();
$sth->finish();

print "ID\tChr\tPos\tEffect\tRef\tP\tBeta\n";

if (defined $beta_column || defined $or_column || defined $z_column)
{
		my $sql  = "select t.VAR_ID,t.$pvalue_column," . (defined $beta_column ? "t.$beta_column" : (defined $or_column ? "LOG(t.$or_column)" : "t.$z_column")) . " from $dataset t where $pvalue_column < $min_pvalue";
		print STDERR "$sql\n";
		$sth = $dbh->prepare($sql);
		$sth->execute();
		while(my ($var_id,$pvalue,$effect)=$sth->fetchrow()){
				my ($chrom,$pos,$ref,$effect_allele) = split(/_/,$var_id);
				next unless $effect_allele;
				$pvalue = "NA" if (!defined $pvalue || $pvalue == "");
				$effect = "NA" if (!defined $effect || $effect == "");
				print $var_id . "\t" . $chrom . "\t" . $pos . "\t" . $effect_allele . "\t" . $ref . "\t" . $pvalue. "\t" . $effect . "\n";
		}
		$sth->finish();
}
=head1 NAME

 get_significant_associations.pl 

=head1 SYNOPSIS

 extracts significant associations from tables

=head1 DESCRIPTION

 

=head1 OPTIONS

 debug                       set logger to DEBUG
 trace                       set logger to TRACE
 help                        help
 man                         manual

 dbi-configuration-file


=head1 EXAMPLE

 
 
=cut
